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


def compute_export_paths(storage_url, lustre_mount_path, gcs_base_uri):
    """Compute the Lustre and GCS paths for export.

    Arguments:
        storage_url (str): Full path to analysis output directory on Lustre.
        lustre_mount_path (str): Lustre mount path prefix (e.g., "/lustre").
        gcs_base_uri (str): Base GCS bucket URI (e.g., "gs://my-bucket").

    Returns:
        tuple: (lustre_path, gcs_path_uri) for the export command.

    Example:
        >>> compute_export_paths("/lustre/analyses/00/01/123", "/lustre", "gs://bucket")
        ('/analyses/00/01/123/', 'gs://bucket/analyses/00/01/123/')
    """
    # Strip the mount path prefix to get the relative Lustre path
    if storage_url.startswith(lustre_mount_path):
        lustre_path = storage_url[len(lustre_mount_path) :]
    else:
        lustre_path = storage_url

    # Ensure lustre_path starts with /
    if not lustre_path.startswith("/"):
        lustre_path = "/" + lustre_path

    # Ensure paths end with / for directory export
    if not lustre_path.endswith("/"):
        lustre_path = lustre_path + "/"

    # Build GCS target path
    gcs_base = gcs_base_uri.rstrip("/")
    gcs_path_uri = f"{gcs_base}{lustre_path}"

    return lustre_path, gcs_path_uri


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


def run_export(storage_url, delete_after=None):
    """Run the full export process from Lustre to GCS.

    This is the main entry point for the lustre-export CLI command.

    Arguments:
        storage_url (str): Full path to analysis output directory on Lustre.
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
        "lustre_mount_path",
        "gcs_base_uri",
    ]
    missing = [s for s in required if not gcp_config.get(s)]
    if missing:
        raise GCPLustreExportError(f"Missing required GCP settings: {missing}")

    # Compute paths
    lustre_path, gcs_path = compute_export_paths(
        storage_url,
        gcp_config["lustre_mount_path"],
        gcp_config["gcs_base_uri"],
    )

    # Initiate export
    operation_name = initiate_export(gcp_config, lustre_path, gcs_path)

    # Wait for completion
    wait_for_export(gcp_config, operation_name)

    # Delete scratch if configured
    should_delete = (
        delete_after
        if delete_after is not None
        else gcp_config.get("lustre_delete_after_export", True)
    )

    if should_delete:
        click.echo(
            f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Deleting scratch data from Lustre..."
        )
        try:
            shutil.rmtree(storage_url)
        except Exception as e:
            click.echo(
                f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] WARNING: Failed to delete scratch: {e}"
            )
        else:
            click.echo(
                f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Scratch data deleted successfully"
            )

    click.echo(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] GCP Lustre export complete")


def get_export_command_for_script(analysis):
    """Get the CLI command string to be embedded in the analysis script.

    This is called from AbstractApplication.write_command_script() when
    GCP Lustre export is enabled.

    Arguments:
        analysis (dict): Analysis instance.

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
        "lustre_mount_path",
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

    storage_url = analysis["storage_url"]
    return f"isabl lustre-export --storage-url {storage_url}"
