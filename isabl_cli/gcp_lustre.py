"""GCP Lustre export functionality for isabl_cli.

This module provides functionality to export analysis results from Google Cloud
Managed Lustre scratch storage to Google Cloud Storage (GCS) after pipeline completion.
"""

import json
import shutil
import subprocess
import time

import click

from isabl_cli.settings import system_settings


class GCPLustreExportError(Exception):
    """Exception raised when GCP Lustre export fails."""

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


def normalize_lustre_path(lustre_path):
    """Normalize the Lustre path for the gcloud command.

    The gcloud command expects a path relative to the Lustre filesystem root,
    not an absolute system path.

    Arguments:
        lustre_path (str): Path on Lustre filesystem (e.g., "/havasove/isabl/3").

    Returns:
        str: Normalized Lustre path for gcloud command (e.g., "/havasove/isabl/3/").
    """
    path = lustre_path

    # Ensure path starts with /
    if not path.startswith("/"):
        path = "/" + path

    # Ensure path ends with / for directory export
    if not path.endswith("/"):
        path = path + "/"

    return path


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
        time.sleep(poll_interval)

        click.echo(
            f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Checking export status "
            f"(attempt {attempt}/{max_attempts})..."
        )

        try:
            done, error = check_export_status(gcp_config, operation_name)
        except GCPLustreExportError as e:
            click.echo(
                f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] WARNING: {e}, retrying..."
            )
            continue

        if done:
            if error:
                raise GCPLustreExportError(f"Export operation failed: {error}")
            click.echo(
                f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Export completed successfully!"
            )
            return

        click.echo(
            f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Export still in progress..."
        )

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
    normalized_lustre_path = normalize_lustre_path(lustre_path)

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
    return f"isabl lustre-export --lustre-path {lustre_path} --analysis-pk {analysis_pk}"
