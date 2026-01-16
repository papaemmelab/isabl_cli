"""LustreInputs helper class for apps to manage input imports from GCS to Lustre.

This module provides a convenient API for isabl applications to register input files
that need to be imported from GCS to Lustre scratch storage before pipeline execution.
"""

import hashlib
import json
import os

from isabl_cli.gcp_lustre import get_gcp_config, GCPLustreImportError


class LustreInputs:
    """Helper for apps to manage input imports from GCS to Lustre.

    Input paths in isabl are stored as local gcsfuse mount paths (e.g., /mnt/gcsfuse/...).
    This class converts them to GCS URIs for import and provides Lustre-local paths.

    Usage in get_command():
        lustre = LustreInputs(analysis, settings)
        lustre.add("/mnt/gcsfuse/data/file.bam")  # gcsfuse path
        local_path = lustre.get("/mnt/gcsfuse/data/file.bam")  # returns /scratch/123/inputs/file.bam

    Example:
        def get_command(self, analysis, inputs, settings):
            lustre = LustreInputs(analysis, settings)

            # Register raw data files
            for target in analysis["targets"]:
                for raw_file in target["raw_data"]:
                    lustre.add(raw_file["file_url"])

            # Register dependency results
            lustre.add(inputs["bam_file"])

            # Get Lustre-local paths
            bam_lustre = lustre.get(inputs["bam_file"])

            # Store for write_command_script to access
            self._lustre_inputs = lustre

            return f"my_pipeline --bam {bam_lustre}"
    """

    def __init__(self, analysis, settings=None):
        """Initialize LustreInputs helper.

        Arguments:
            analysis (dict): Analysis instance from API.
            settings (object, optional): Application settings (unused but kept for API consistency).
        """
        self.analysis = analysis
        self.settings = settings
        self._original_to_lustre = {}  # Maps original path -> Lustre relative path
        self._import_specs = []  # List of (gcs_uri, lustre_path) for batch import

        # Get config for path conversion
        gcp_config = get_gcp_config()
        self.lustre_mount = gcp_config.get("lustre_mount_path", "/scratch")
        self.gcsfuse_mount = gcp_config.get("gcsfuse_mount_path")  # e.g., "/mnt/gcsfuse"
        self.gcs_input_uri = gcp_config.get("gcs_input_uri")  # e.g., "gs://input-bucket"
        self.import_enabled = gcp_config.get("lustre_import_enabled", False)

        # Base directory for this analysis's inputs on Lustre (relative to mount)
        self.input_dir = f"/{analysis['pk']}/inputs"

    def _gcsfuse_to_gcs_uri(self, gcsfuse_path):
        """Convert gcsfuse mount path to GCS URI.

        Arguments:
            gcsfuse_path (str): Local gcsfuse path (e.g., /mnt/gcsfuse/data/file.fastq).

        Returns:
            str: GCS URI (e.g., gs://input-bucket/data/file.fastq).

        Raises:
            GCPLustreImportError: If path conversion fails.

        Example:
            /mnt/gcsfuse/data/file.fastq -> gs://input-bucket/data/file.fastq
        """
        # If already a GCS URI, return as-is
        if gcsfuse_path.startswith("gs://"):
            return gcsfuse_path

        if not self.gcsfuse_mount or not self.gcs_input_uri:
            raise GCPLustreImportError(
                "gcsfuse_mount_path and gcs_input_uri must be configured in GCP_CONFIGURATION"
            )

        if not gcsfuse_path.startswith(self.gcsfuse_mount):
            raise GCPLustreImportError(
                f"Path doesn't start with gcsfuse mount ({self.gcsfuse_mount}): {gcsfuse_path}"
            )

        # Extract relative path and build GCS URI
        relative = gcsfuse_path[len(self.gcsfuse_mount) :]
        if not relative.startswith("/"):
            relative = "/" + relative
        return f"{self.gcs_input_uri.rstrip('/')}{relative}"

    def add(self, path):
        """Register a file path to be imported to Lustre.

        Arguments:
            path (str): Local gcsfuse path (e.g., /mnt/gcsfuse/data/file.bam)
                        or GCS URI (e.g., gs://bucket/data/file.bam).

        Returns:
            str: The Lustre-local path where the file will be available.
        """
        if not self.import_enabled:
            # When import is disabled, return the original path
            return path

        if path in self._original_to_lustre:
            # Already registered, return existing Lustre path
            return f"{self.lustre_mount}{self._original_to_lustre[path]}"

        # Convert to GCS URI for import command
        gcs_uri = self._gcsfuse_to_gcs_uri(path)

        # Compute Lustre destination (preserve filename)
        filename = os.path.basename(path)
        lustre_relative = f"{self.input_dir}/{filename}"

        # Handle duplicates by adding hash suffix
        if lustre_relative in self._original_to_lustre.values():
            hash_suffix = hashlib.md5(path.encode()).hexdigest()[:8]
            base, ext = os.path.splitext(filename)
            lustre_relative = f"{self.input_dir}/{base}_{hash_suffix}{ext}"

        self._original_to_lustre[path] = lustre_relative
        self._import_specs.append((gcs_uri, lustre_relative))

        return f"{self.lustre_mount}{lustre_relative}"

    def get(self, path):
        """Get the Lustre-local path for a registered file.

        Arguments:
            path (str): The original path passed to add().

        Returns:
            str: Full Lustre path (e.g., /scratch/123/inputs/file.bam).

        Raises:
            ValueError: If path was not previously registered with add().
        """
        if not self.import_enabled:
            # When import is disabled, return the original path
            return path

        if path not in self._original_to_lustre:
            raise ValueError(f"Path not registered: {path}. Call add() first.")
        # Return full path including mount
        return f"{self.lustre_mount}{self._original_to_lustre[path]}"

    def get_import_specs(self):
        """Get the list of import specifications.

        Returns:
            list: List of (gcs_path, lustre_path) tuples.
        """
        return list(self._import_specs)

    def get_import_command(self):
        """Get CLI command to run imports (for embedding in script).

        Returns:
            str: CLI command string, or empty string if no imports needed.
        """
        if not self._import_specs or not self.import_enabled:
            return ""

        specs_json = json.dumps(self._import_specs)
        # Use single quotes around JSON and escape any single quotes within
        escaped_json = specs_json.replace("'", "'\"'\"'")
        return f"isabl lustre-import --specs '{escaped_json}'"

    def __len__(self):
        """Return the number of registered files."""
        return len(self._import_specs)

    def __repr__(self):
        """String representation of LustreInputs."""
        return (
            f"LustreInputs(analysis_pk={self.analysis['pk']}, "
            f"files={len(self._import_specs)}, enabled={self.import_enabled})"
        )
