"""GCP Lustre import/export functionality for isabl_cli.

This module provides functionality to:
- Import input data from Google Cloud Storage (GCS) to Lustre scratch before pipeline runs
- Export analysis results from Lustre scratch to GCS after pipeline completion

This enables high-performance computing on GCP Slurm clusters using Lustre scratch storage.
"""

import fcntl
import json
import os
import shutil
import subprocess
import sys
import time

import click

from isabl_cli.settings import system_settings


class GCPLustreExportError(Exception):
    """Exception raised when GCP Lustre export fails."""

    pass


class GCPLustreImportError(Exception):
    """Exception raised when GCP Lustre import fails."""

    pass


def get_gcp_config():
    """Get GCP configuration from system settings.

    Returns:
        dict: GCP configuration dictionary, or empty dict if not configured.
    """
    return getattr(system_settings, "GCP_CONFIGURATION", None) or {}


def get_gcs_path_from_analysis(analysis_pk, gcs_base_uri, base_storage_directory):
    """Compute GCS target path from analysis storage_url.

    Fetches the analysis from the API to get its canonical storage_url,
    extracts the relative path, and combines with gcs_base_uri.

    Arguments:
        analysis_pk (int): Analysis primary key to fetch from API.
        gcs_base_uri (str): Base GCS bucket URI (e.g., "gs://my-bucket").
        base_storage_directory (str): Base storage directory (e.g., "/datalake").

    Returns:
        str: GCS path URI (e.g., "gs://bucket/analyses/00/01/123/").

    Raises:
        GCPLustreExportError: If analysis cannot be fetched or path cannot be computed.

    Example:
        >>> get_gcs_path_from_analysis(123, "gs://bucket", "/datalake")
        'gs://bucket/analyses/00/01/123/'
    """
    from isabl_cli import api

    try:
        analysis = api.get_instance("analyses", analysis_pk, fields="storage_url")
    except Exception as e:
        raise GCPLustreExportError(f"Failed to fetch analysis {analysis_pk}: {e}")

    storage_url = analysis.get("storage_url")
    if not storage_url:
        raise GCPLustreExportError(f"Analysis {analysis_pk} has no storage_url")

    # Extract relative path from storage_url by removing base_storage_directory prefix
    base = base_storage_directory.rstrip("/")
    if storage_url.startswith(base):
        relative_path = storage_url[len(base) :]
    else:
        # If storage_url doesn't start with base, use the full path
        relative_path = storage_url

    # Ensure relative_path starts with /
    if not relative_path.startswith("/"):
        relative_path = "/" + relative_path

    # Ensure path ends with / for directory export
    if not relative_path.endswith("/"):
        relative_path = relative_path + "/"

    # Build GCS target path
    gcs_base = gcs_base_uri.rstrip("/")
    return f"{gcs_base}{relative_path}"


def normalize_lustre_path(lustre_path, lustre_mount_path=None):
    """Normalize the Lustre path for the gcloud command.

    The gcloud command expects a path relative to the Lustre filesystem root,
    not an absolute system path. This function strips the mount path prefix
    if provided.

    Arguments:
        lustre_path (str): Path on Lustre filesystem (e.g., "/scratch/isabl/data/analyses/00/01/1").
        lustre_mount_path (str, optional): The Lustre mount point (e.g., "/scratch").
            If provided, this prefix will be stripped from lustre_path.

    Returns:
        str: Normalized Lustre path for gcloud command (e.g., "/isabl/data/analyses/00/01/1/").

    Example:
        >>> normalize_lustre_path("/scratch/isabl/data/analyses/00/01/1", "/scratch")
        '/isabl/data/analyses/00/01/1/'
    """
    path = lustre_path

    # Strip the mount path prefix if provided
    if lustre_mount_path:
        mount = lustre_mount_path.rstrip("/")
        if path.startswith(mount):
            path = path[len(mount):]

    # Ensure path starts with /
    if not path.startswith("/"):
        path = "/" + path

    # Ensure path ends with / for directory export
    if not path.endswith("/"):
        path = path + "/"

    return path


def _retry_with_fixed_interval(func, max_retries=180, retry_interval=60, operation_name="operation", error_class=None):
    """Retry a function with fixed interval.

    Used for waiting on resources that may be busy (e.g., only one Lustre export/import
    allowed at a time). Retries at a fixed interval until success or max retries.

    Arguments:
        func (callable): Function to retry. Should raise exception on failure.
        max_retries (int): Maximum number of attempts (default: 180 = 3 hours at 60s interval).
        retry_interval (int): Delay between retries in seconds (default: 60).
        operation_name (str): Name for logging.
        error_class (type): Exception class to raise on exhaustion (default: GCPLustreExportError).

    Returns:
        Any: Return value of func on success.

    Raises:
        GCPLustreExportError or error_class: If all retries exhausted.
    """
    if error_class is None:
        error_class = GCPLustreExportError
    for attempt in range(1, max_retries + 1):
        try:
            return func()
        except Exception as e:
            if attempt == max_retries:
                raise error_class(
                    f"{operation_name} failed after {max_retries} attempts: {e}"
                )

            click.echo(
                f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] "
                f"{operation_name} attempt {attempt}/{max_retries} failed: {e}. "
                f"Retrying in {retry_interval}s..."
            )
            time.sleep(retry_interval)


def initiate_export(gcp_config, lustre_path, gcs_path):
    """Initiate an async export from Lustre to GCS.

    Arguments:
        gcp_config (dict): GCP configuration dictionary.
        lustre_path (str): Lustre path to export (e.g., "/analyses/00/01/123/").
        gcs_path (str): GCS destination URI (e.g., "gs://bucket/analyses/00/01/123/").

    Returns:
        str: Operation name for tracking the export.

    Raises:
        GCPLustreExportError: If export initiation fails.
    """
    cmd = [
        "gcloud",
        "lustre",
        "instances",
        "export-data",
        gcp_config["lustre_instance"],
        f"--location={gcp_config['lustre_location']}",
        f"--gcs-path-uri={gcs_path}",
        f"--lustre-path={lustre_path}",
        "--async",
        f"--project={gcp_config['lustre_project']}",
        "--format=json",
    ]

    click.echo(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Initiating Lustre export...")
    click.echo(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Lustre path: {lustre_path}")
    click.echo(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] GCS target: {gcs_path}")

    def _run_export_command():
        """Inner function for retry logic."""
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True,
            )
        except subprocess.CalledProcessError as e:
            raise GCPLustreExportError(
                f"Failed to initiate Lustre export: {e.stderr or e.stdout}"
            )

        # Parse JSON output to extract operation name
        try:
            output = json.loads(result.stdout)
            operation_name = output.get("name")
            if not operation_name:
                raise GCPLustreExportError(
                    f"Could not extract operation name from output: {result.stdout}"
                )
        except json.JSONDecodeError as e:
            raise GCPLustreExportError(
                f"Failed to parse export output as JSON: {result.stdout}"
            )

        return operation_name

    # Retry with fixed interval (default: every 60s for 3 hours)
    # Only one Lustre export can run at a time, so we wait for any in-progress export
    max_retries = gcp_config.get("lustre_export_max_retries", 180)
    retry_interval = gcp_config.get("lustre_export_retry_interval", 60)

    operation_name = _retry_with_fixed_interval(
        _run_export_command,
        max_retries=max_retries,
        retry_interval=retry_interval,
        operation_name="Lustre export initiation"
    )

    click.echo(
        f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Export initiated. Operation: {operation_name}"
    )
    return operation_name


def check_export_status(gcp_config, operation_name):
    """Check the status of an export operation.

    Arguments:
        gcp_config (dict): GCP configuration dictionary.
        operation_name (str): Operation name to check.

    Returns:
        tuple: (done, error) where done is bool and error is str or None.

    Raises:
        GCPLustreExportError: If status check fails.
    """
    cmd = [
        "gcloud",
        "lustre",
        "operations",
        "describe",
        operation_name,
        f"--location={gcp_config['lustre_location']}",
        f"--project={gcp_config['lustre_project']}",
        "--format=json",
    ]

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True,
        )
    except subprocess.CalledProcessError as e:
        raise GCPLustreExportError(
            f"Failed to check operation status: {e.stderr or e.stdout}"
        )

    try:
        output = json.loads(result.stdout)
    except json.JSONDecodeError as e:
        raise GCPLustreExportError(
            f"Failed to parse status output as JSON: {result.stdout}"
        )

    done = output.get("done", False)
    error = output.get("error")

    if error:
        error_msg = error.get("message", str(error))
        return done, error_msg

    return done, None


def wait_for_export(gcp_config, operation_name):
    """Poll until export completes or times out.

    Arguments:
        gcp_config (dict): GCP configuration dictionary.
        operation_name (str): Operation name to wait for.

    Raises:
        GCPLustreExportError: If export fails or times out.
    """
    poll_interval = gcp_config.get("lustre_poll_interval", 30)
    max_attempts = gcp_config.get("lustre_max_poll_attempts", 360)

    for attempt in range(1, max_attempts + 1):
        click.echo(
            f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Sleeping {poll_interval}s before check..."
        )
        sys.stdout.flush()
        sys.stderr.flush()

        time.sleep(poll_interval)

        click.echo(
            f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Checking export status "
            f"(attempt {attempt}/{max_attempts})..."
        )
        sys.stdout.flush()
        sys.stderr.flush()

        try:
            done, error = check_export_status(gcp_config, operation_name)
        except GCPLustreExportError as e:
            click.echo(
                f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] WARNING: {e}, retrying..."
            )
            sys.stdout.flush()
            sys.stderr.flush()
            continue

        if done:
            if error:
                raise GCPLustreExportError(f"Export operation failed: {error}")
            click.echo(
                f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Export completed successfully!"
            )
            sys.stdout.flush()
            sys.stderr.flush()
            return

        click.echo(
            f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Export still in progress..."
        )
        sys.stdout.flush()
        sys.stderr.flush()

    raise GCPLustreExportError(f"Export timed out after {max_attempts} attempts")


def run_export(lustre_path, analysis_pk, delete_after=None):
    """Run the full export process from Lustre to GCS.

    This is the main entry point for the lustre-export CLI command.

    Arguments:
        lustre_path (str): Path on Lustre scratch where files are located.
        analysis_pk (int): Analysis primary key to look up GCS target path.
        delete_after (bool or None): Whether to delete scratch after export.
            If None, uses config value.

    Raises:
        GCPLustreExportError: If export fails.
        SystemExit: On error (for CLI usage).
    """
    gcp_config = get_gcp_config()

    if not gcp_config.get("lustre_export_enabled"):
        raise GCPLustreExportError("GCP Lustre export is not enabled")

    # Validate required settings
    required = [
        "lustre_instance",
        "lustre_location",
        "lustre_project",
        "gcs_base_uri",
    ]
    missing = [s for s in required if not gcp_config.get(s)]
    if missing:
        raise GCPLustreExportError(f"Missing required GCP settings: {missing}")

    # Compute GCS target path from analysis
    gcs_path = get_gcs_path_from_analysis(
        analysis_pk,
        gcp_config["gcs_base_uri"],
        system_settings.BASE_STORAGE_DIRECTORY,
    )

    # Normalize Lustre path for gcloud command (relative to Lustre filesystem root)
    normalized_lustre_path = normalize_lustre_path(
        lustre_path, gcp_config.get("lustre_mount_path")
    )

    # Sync filesystem to ensure all writes are flushed before export
    click.echo(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Syncing filesystem...")
    subprocess.run(["sync"], check=True)

    # Initiate export
    operation_name = initiate_export(gcp_config, normalized_lustre_path, gcs_path)

    # Wait for completion
    wait_for_export(gcp_config, operation_name)

    # Delete scratch if configured
    should_delete = (
        delete_after
        if delete_after is not None
        else gcp_config.get("lustre_delete_after_export", True)
    )

    # Compute full system path for deletion by prepending mount path
    full_lustre_path = lustre_path
    if gcp_config.get("lustre_mount_path"):
        mount = gcp_config["lustre_mount_path"].rstrip("/")
        # lustre_path is relative to mount, so construct full path
        lustre_relative = lustre_path.lstrip("/")
        full_lustre_path = f"{mount}/{lustre_relative}"

    if should_delete:
        click.echo(
            f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Deleting scratch data from Lustre..."
        )
        try:
            shutil.rmtree(full_lustre_path)
        except Exception as e:
            click.echo(
                f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] WARNING: Failed to delete scratch: {e}"
            )
        else:
            click.echo(
                f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Scratch data deleted successfully"
            )

    click.echo(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] GCP Lustre export complete")


def get_export_command_for_script(analysis, lustre_path):
    """Get the CLI command string to be embedded in the analysis script.

    This is called from AbstractApplication.write_command_script() when
    GCP Lustre export is enabled.

    Arguments:
        analysis (dict): Analysis instance.
        lustre_path (str): Path on Lustre where output files are written.
            This should be provided by the user via get_lustre_output_path().

    Returns:
        str: CLI command string, or empty string if export disabled.
    """
    gcp_config = get_gcp_config()

    if not gcp_config.get("lustre_export_enabled"):
        return ""

    # Validate required settings
    required = [
        "lustre_instance",
        "lustre_location",
        "lustre_project",
        "gcs_base_uri",
    ]

    missing = [s for s in required if not gcp_config.get(s)]
    if missing:
        click.secho(
            f"GCP Lustre export enabled but missing settings: {missing}",
            err=True,
            fg="yellow",
        )
        return ""

    analysis_pk = analysis["pk"]

    # Build command with explicit --delete-after flag based on config
    delete_flag = (
        "--delete-after"
        if gcp_config.get("lustre_delete_after_export", True)
        else "--no-delete-after"
    )

    # Redirect output to /tmp to avoid FILE_MODIFIED_FAILURE
    # The export command outputs progress updates which would modify head_job.log
    # inside the directory being exported, causing the export to fail.
    export_log = f"/tmp/lustre_export_{analysis_pk}.log"
    return f"isabl lustre-export --lustre-path {lustre_path} --analysis-pk {analysis_pk} {delete_flag} > {export_log} 2>&1"


# ============================================================================
# IMPORT FUNCTIONS
# ============================================================================


def initiate_import(gcp_config, gcs_path, lustre_path):
    """Initiate an async import from GCS to Lustre.

    Arguments:
        gcp_config (dict): GCP configuration dictionary.
        gcs_path (str): GCS source URI (e.g., "gs://bucket/data/file.fastq").
        lustre_path (str): Lustre destination path (e.g., "/123/inputs/file.fastq").

    Returns:
        str: Operation name for tracking the import.

    Raises:
        GCPLustreImportError: If import initiation fails.
    """
    cmd = [
        "gcloud",
        "lustre",
        "instances",
        "import-data",
        gcp_config["lustre_instance"],
        f"--location={gcp_config['lustre_location']}",
        f"--gcs-path-uri={gcs_path}",
        f"--lustre-path={lustre_path}",
        "--async",
        f"--project={gcp_config['lustre_project']}",
        "--format=json",
    ]

    click.echo(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Initiating Lustre import...")
    click.echo(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] GCS source: {gcs_path}")
    click.echo(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Lustre target: {lustre_path}")

    def _run_import_command():
        """Inner function for retry logic."""
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True,
            )
        except subprocess.CalledProcessError as e:
            raise GCPLustreImportError(
                f"Failed to initiate Lustre import: {e.stderr or e.stdout}"
            )

        # Parse JSON output to extract operation name
        try:
            output = json.loads(result.stdout)
            operation_name = output.get("name")
            if not operation_name:
                raise GCPLustreImportError(
                    f"Could not extract operation name from output: {result.stdout}"
                )
        except json.JSONDecodeError as e:
            raise GCPLustreImportError(
                f"Failed to parse import output as JSON: {result.stdout}"
            )

        return operation_name

    # Retry with fixed interval (only one Lustre import/export can run at a time)
    max_retries = gcp_config.get("lustre_import_max_retries", 180)
    retry_interval = gcp_config.get("lustre_import_retry_interval", 60)

    operation_name = _retry_with_fixed_interval(
        _run_import_command,
        max_retries=max_retries,
        retry_interval=retry_interval,
        operation_name="Lustre import initiation",
        error_class=GCPLustreImportError,
    )

    click.echo(
        f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Import initiated. Operation: {operation_name}"
    )
    return operation_name


def check_import_status(gcp_config, operation_name):
    """Check the status of an import operation.

    Arguments:
        gcp_config (dict): GCP configuration dictionary.
        operation_name (str): Operation name to check.

    Returns:
        tuple: (done, error) where done is bool and error is str or None.

    Raises:
        GCPLustreImportError: If status check fails.
    """
    cmd = [
        "gcloud",
        "lustre",
        "operations",
        "describe",
        operation_name,
        f"--location={gcp_config['lustre_location']}",
        f"--project={gcp_config['lustre_project']}",
        "--format=json",
    ]

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True,
        )
    except subprocess.CalledProcessError as e:
        raise GCPLustreImportError(
            f"Failed to check operation status: {e.stderr or e.stdout}"
        )

    try:
        output = json.loads(result.stdout)
    except json.JSONDecodeError as e:
        raise GCPLustreImportError(
            f"Failed to parse status output as JSON: {result.stdout}"
        )

    done = output.get("done", False)
    error = output.get("error")

    if error:
        error_msg = error.get("message", str(error))
        return done, error_msg

    return done, None


def wait_for_imports(gcp_config, operations):
    """Poll until all import operations complete or timeout.

    Arguments:
        gcp_config (dict): GCP configuration dictionary.
        operations (list): List of operation names to wait for.

    Raises:
        GCPLustreImportError: If any import fails or times out.
    """
    if not operations:
        return

    poll_interval = gcp_config.get("lustre_poll_interval", 30)
    max_attempts = gcp_config.get("lustre_max_poll_attempts", 360)

    pending = set(operations)
    failed = []

    for attempt in range(1, max_attempts + 1):
        time.sleep(poll_interval)

        click.echo(
            f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Checking import status "
            f"(attempt {attempt}/{max_attempts}, {len(pending)} pending)..."
        )

        still_pending = set()
        for operation_name in pending:
            try:
                done, error = check_import_status(gcp_config, operation_name)
            except GCPLustreImportError as e:
                click.echo(
                    f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] WARNING: {e}, retrying..."
                )
                still_pending.add(operation_name)
                continue

            if done:
                if error:
                    failed.append((operation_name, error))
                    click.echo(
                        f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Import FAILED: {operation_name}: {error}"
                    )
                else:
                    click.echo(
                        f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Import completed: {operation_name}"
                    )
            else:
                still_pending.add(operation_name)

        pending = still_pending

        if not pending:
            break

        click.echo(
            f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {len(pending)} imports still in progress..."
        )

    if failed:
        errors = "; ".join(f"{op}: {err}" for op, err in failed)
        raise GCPLustreImportError(f"Import operation(s) failed: {errors}")

    if pending:
        raise GCPLustreImportError(
            f"Import timed out after {max_attempts} attempts. "
            f"Still pending: {list(pending)}"
        )

    click.echo(
        f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] All imports completed successfully!"
    )


def run_batch_import(import_specs):
    """Import multiple directories from GCS to Lustre in parallel.

    This is the main entry point for the lustre-import CLI command.
    Note: gcloud lustre import-data only supports directory-level imports.

    Arguments:
        import_specs (list): List of [gcs_dir, lustre_dir] pairs (directories, not files).

    Raises:
        GCPLustreImportError: If any import fails.
    """
    gcp_config = get_gcp_config()

    if not gcp_config.get("lustre_import_enabled"):
        raise GCPLustreImportError("GCP Lustre import is not enabled")

    # Validate required settings
    required = [
        "lustre_instance",
        "lustre_location",
        "lustre_project",
    ]
    missing = [s for s in required if not gcp_config.get(s)]
    if missing:
        raise GCPLustreImportError(f"Missing required GCP settings: {missing}")

    if not import_specs:
        click.echo(
            f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] No directories to import, skipping."
        )
        return

    click.echo(
        f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Starting batch import of {len(import_specs)} directory(ies)..."
    )

    # Initiate all imports
    operations = []
    for gcs_path, lustre_path in import_specs:
        # Ensure GCS path ends with / for directory import
        if not gcs_path.endswith("/"):
            gcs_path = gcs_path + "/"
        # Normalize lustre path for gcloud command
        normalized_path = normalize_lustre_path(
            lustre_path, gcp_config.get("lustre_mount_path")
        )
        operation_name = initiate_import(gcp_config, gcs_path, normalized_path)
        operations.append(operation_name)

    # Wait for all to complete
    wait_for_imports(gcp_config, operations)

    click.echo(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] GCP Lustre batch import complete")


# ============================================================================
# SHARED IMPORT/CLEANUP FUNCTIONS (with reference counting and flock locking)
# ============================================================================


def run_shared_import(analysis_pk, import_specs, strategy="lustre-import"):
    """Import shared inputs with reference counting and flock locking.

    For each (gcs_dir, lustre_dir) in import_specs:
    1. Acquire exclusive flock on {lustre_dir}.lock
    2. Check if data already exists (.import_complete marker)
    3. If not, import data via the chosen strategy
    4. Add ref file: {lustre_dir}.refs/{analysis_pk}
    5. Release lock

    Arguments:
        analysis_pk (int): Analysis primary key (for ref file naming).
        import_specs (list): List of [gcs_dir, lustre_dir] pairs.
            lustre_dir is relative to the lustre mount point.
        strategy (str): Import strategy - "lustre-import" or "rsync".

    Raises:
        GCPLustreImportError: If import fails.
    """
    gcp_config = get_gcp_config()

    if not gcp_config.get("lustre_import_enabled"):
        raise GCPLustreImportError("GCP Lustre import is not enabled")

    lustre_mount = gcp_config.get("lustre_mount_path", "/scratch")

    click.echo(
        f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Starting shared import for "
        f"analysis {analysis_pk} ({len(import_specs)} directories, strategy={strategy})"
    )

    for gcs_path, lustre_rel_path in import_specs:
        full_data_dir = f"{lustre_mount}/{lustre_rel_path.strip('/')}"
        refs_dir = f"{full_data_dir}.refs"
        lock_file = f"{full_data_dir}.lock"
        marker_file = os.path.join(full_data_dir, ".import_complete")

        # Ensure directories exist (create data dir before import so we own it)
        os.makedirs(full_data_dir, exist_ok=True)
        os.makedirs(refs_dir, exist_ok=True)

        click.echo(
            f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] "
            f"Acquiring lock for {lustre_rel_path}..."
        )

        lock_fd = open(lock_file, "w")
        try:
            fcntl.flock(lock_fd, fcntl.LOCK_EX)
            click.echo(
                f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Lock acquired."
            )

            # Check if data already imported
            if os.path.exists(marker_file):
                click.echo(
                    f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] "
                    f"Data already imported at {full_data_dir}, skipping import."
                )
            else:
                click.echo(
                    f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] "
                    f"Importing {gcs_path} -> {full_data_dir} (strategy={strategy})"
                )

                if strategy == "rsync":
                    _rsync_import(gcs_path, full_data_dir)
                else:
                    _lustre_import_single(gcp_config, gcs_path, lustre_rel_path)

                # Write completion marker
                os.makedirs(full_data_dir, exist_ok=True)
                with open(marker_file, "w") as f:
                    f.write(gcs_path)

                click.echo(
                    f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Import complete."
                )

            # Add reference file (always, even if data already existed)
            ref_file = os.path.join(refs_dir, str(analysis_pk))
            with open(ref_file, "w") as f:
                f.write("")

            click.echo(
                f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] "
                f"Added ref for analysis {analysis_pk}."
            )
        finally:
            fcntl.flock(lock_fd, fcntl.LOCK_UN)
            lock_fd.close()

    click.echo(
        f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] "
        f"Shared import complete for analysis {analysis_pk}."
    )


def run_shared_cleanup(analysis_pk, import_specs):
    """Release reference and optionally delete shared import data.

    For each (gcs_dir, lustre_dir) in import_specs:
    1. Acquire exclusive flock on {lustre_dir}.lock
    2. Remove ref file: {lustre_dir}.refs/{analysis_pk}
    3. If refs dir is empty, delete data directory and refs directory
    4. Release lock

    Arguments:
        analysis_pk (int): Analysis primary key.
        import_specs (list): List of [gcs_dir, lustre_dir] pairs.

    Raises:
        GCPLustreImportError: If cleanup encounters a critical error.
    """
    gcp_config = get_gcp_config()
    lustre_mount = gcp_config.get("lustre_mount_path", "/scratch")

    click.echo(
        f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Starting shared cleanup for "
        f"analysis {analysis_pk} ({len(import_specs)} directories)"
    )

    for gcs_path, lustre_rel_path in import_specs:
        full_data_dir = f"{lustre_mount}/{lustre_rel_path.strip('/')}"
        refs_dir = f"{full_data_dir}.refs"
        lock_file = f"{full_data_dir}.lock"

        if not os.path.exists(refs_dir):
            click.echo(
                f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] "
                f"No refs dir found for {lustre_rel_path}, skipping."
            )
            continue

        lock_fd = open(lock_file, "w")
        try:
            fcntl.flock(lock_fd, fcntl.LOCK_EX)

            # Remove our ref
            ref_file = os.path.join(refs_dir, str(analysis_pk))
            if os.path.exists(ref_file):
                os.remove(ref_file)
                click.echo(
                    f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] "
                    f"Removed ref for analysis {analysis_pk} from {lustre_rel_path}."
                )

            # Check if any refs remain
            remaining = os.listdir(refs_dir)
            if not remaining:
                click.echo(
                    f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] "
                    f"No more refs, deleting shared data at {full_data_dir}."
                )
                shutil.rmtree(full_data_dir, ignore_errors=True)
                shutil.rmtree(refs_dir, ignore_errors=True)
            else:
                click.echo(
                    f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] "
                    f"{len(remaining)} ref(s) remaining for {lustre_rel_path}, keeping data."
                )
        finally:
            fcntl.flock(lock_fd, fcntl.LOCK_UN)
            lock_fd.close()

            # Remove lock file if data was deleted (no refs remain)
            if not os.path.exists(refs_dir):
                try:
                    os.remove(lock_file)
                except OSError:
                    pass

    click.echo(
        f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] "
        f"Shared cleanup complete for analysis {analysis_pk}."
    )


def _rsync_import(gcs_path, local_dest):
    """Import from GCS to local path using gsutil rsync.

    Arguments:
        gcs_path (str): GCS URI (e.g., gs://bucket/data/dir/).
        local_dest (str): Local destination directory.

    Raises:
        GCPLustreImportError: If gsutil rsync fails.
    """
    os.makedirs(local_dest, exist_ok=True)
    cmd = [
        "gsutil", "-m", "rsync", "-r",
        gcs_path.rstrip("/") + "/",
        local_dest,
    ]
    click.echo(
        f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Running: {' '.join(cmd)}"
    )
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        raise GCPLustreImportError(
            f"gsutil rsync failed: {e.stderr or e.stdout}"
        )


def _lustre_import_single(gcp_config, gcs_path, lustre_rel_path):
    """Import a single directory using gcloud lustre instances import-data.

    Wraps initiate_import + wait_for_imports for a single directory.

    Arguments:
        gcp_config (dict): GCP configuration dictionary.
        gcs_path (str): GCS source URI.
        lustre_rel_path (str): Lustre path relative to mount point.

    Raises:
        GCPLustreImportError: If import fails.
    """
    if not gcs_path.endswith("/"):
        gcs_path = gcs_path + "/"
    normalized_path = normalize_lustre_path(
        lustre_rel_path, gcp_config.get("lustre_mount_path")
    )
    operation_name = initiate_import(gcp_config, gcs_path, normalized_path)
    wait_for_imports(gcp_config, [operation_name])
