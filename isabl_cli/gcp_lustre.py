"""GCP Lustre export functionality for isabl_cli.

This module provides functionality to export analysis results from Google Cloud
Managed Lustre scratch storage to Google Cloud Storage (GCS) after pipeline completion.
"""

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


def build_export_script(analysis, gcp_config):
    """Build a bash script for exporting data from Lustre to GCS.

    This script:
    1. Initiates an async export from Lustre to GCS
    2. Polls for completion
    3. Optionally deletes scratch data on success
    4. Exits with error code on failure

    Arguments:
        analysis (dict): Analysis instance with 'pk' and 'storage_url' keys.
        gcp_config (dict): GCP configuration dictionary.

    Returns:
        str: Bash script for export.
    """
    storage_url = analysis["storage_url"]
    lustre_instance = gcp_config["lustre_instance"]
    location = gcp_config["lustre_location"]
    project = gcp_config["lustre_project"]
    lustre_mount_path = gcp_config["lustre_mount_path"]
    gcs_base_uri = gcp_config["gcs_base_uri"]
    poll_interval = gcp_config.get("lustre_poll_interval", 30)
    max_poll_attempts = gcp_config.get("lustre_max_poll_attempts", 360)
    delete_after_export = gcp_config.get("lustre_delete_after_export", True)

    lustre_path, gcs_path = compute_export_paths(
        storage_url, lustre_mount_path, gcs_base_uri
    )

    script = f'''
# GCP Lustre Export Script for Analysis {analysis["pk"]}
echo "[$(date)] Starting GCP Lustre export..."
echo "[$(date)] Lustre path: {lustre_path}"
echo "[$(date)] GCS target: {gcs_path}"

# Initiate async export
EXPORT_OUTPUT=$(gcloud lustre instances export-data {lustre_instance} \\
    --location={location} \\
    --gcs-path-uri="{gcs_path}" \\
    --lustre-path="{lustre_path}" \\
    --async \\
    --project {project} \\
    --format=json 2>&1)

EXPORT_EXIT_CODE=$?
if [ $EXPORT_EXIT_CODE -ne 0 ]; then
    echo "[$(date)] ERROR: Failed to initiate Lustre export"
    echo "$EXPORT_OUTPUT"
    exit 1
fi

# Extract operation name from output
OPERATION_NAME=$(echo "$EXPORT_OUTPUT" | grep -o '"name": "[^"]*"' | head -1 | cut -d'"' -f4)

if [ -z "$OPERATION_NAME" ]; then
    echo "[$(date)] ERROR: Could not extract operation name from export output"
    echo "$EXPORT_OUTPUT"
    exit 1
fi

echo "[$(date)] Export initiated. Operation: $OPERATION_NAME"

# Poll for completion
POLL_COUNT=0
MAX_POLLS={max_poll_attempts}
POLL_INTERVAL={poll_interval}

while [ $POLL_COUNT -lt $MAX_POLLS ]; do
    sleep $POLL_INTERVAL
    POLL_COUNT=$((POLL_COUNT + 1))

    echo "[$(date)] Checking export status (attempt $POLL_COUNT/$MAX_POLLS)..."

    STATUS_OUTPUT=$(gcloud lustre operations describe "$OPERATION_NAME" \\
        --location={location} \\
        --project {project} \\
        --format=json 2>&1)

    STATUS_EXIT_CODE=$?
    if [ $STATUS_EXIT_CODE -ne 0 ]; then
        echo "[$(date)] WARNING: Failed to get operation status, retrying..."
        echo "$STATUS_OUTPUT"
        continue
    fi

    # Check if operation is done
    IS_DONE=$(echo "$STATUS_OUTPUT" | grep -o '"done": true')

    if [ -n "$IS_DONE" ]; then
        # Check for errors
        HAS_ERROR=$(echo "$STATUS_OUTPUT" | grep '"error":')

        if [ -n "$HAS_ERROR" ]; then
            echo "[$(date)] ERROR: Export operation failed"
            echo "$STATUS_OUTPUT"
            exit 1
        fi

        echo "[$(date)] Export completed successfully!"
        break
    fi

    echo "[$(date)] Export still in progress..."
done

if [ $POLL_COUNT -ge $MAX_POLLS ]; then
    echo "[$(date)] ERROR: Export timed out after $MAX_POLLS attempts"
    exit 1
fi
'''

    if delete_after_export:
        script += f'''
# Delete scratch data after successful export
echo "[$(date)] Deleting scratch data from Lustre..."
rm -rf "{storage_url}"/*
echo "[$(date)] Scratch data deleted successfully"
'''

    script += '''
echo "[$(date)] GCP Lustre export complete"
'''

    return script


def get_export_command_for_script(analysis):
    """Get the export command string to be embedded in the analysis script.

    This is called from AbstractApplication.write_command_script() when
    GCP Lustre export is enabled.

    Arguments:
        analysis (dict): Analysis instance.

    Returns:
        str: Export command string, or empty string if export disabled.
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

    return build_export_script(analysis, gcp_config)
