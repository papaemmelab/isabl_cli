"""LustreInputs helper class for apps to manage shared input imports from GCS to Lustre.

This module provides a convenient API for isabl applications to register input files
that need to be imported from GCS to Lustre scratch storage before pipeline execution.

Inputs are placed in a shared location on Lustre (/shared_inputs/{md5}/) so that
multiple analyses needing the same data only import once. File-based reference counting
tracks which analyses are using each shared import, enabling safe cleanup.

Note: The gcloud lustre import-data command only supports directory-level imports,
so this class groups files by their parent directory and imports entire directories.
"""

import hashlib
import json
import os

from isabl_cli.gcp_lustre import get_gcp_config, GCPLustreImportError


class LustreInputs:
    """Helper for apps to manage shared input imports from GCS to Lustre.

    Input paths in isabl are stored as local gcsfuse mount paths (e.g., /mnt/gcsfuse/...).
    This class converts them to GCS URIs for import and provides shared Lustre-local paths.

    Inputs are stored in a shared location: {lustre_mount}{shared_inputs_path}/{md5_of_gcs_dir}/
    Multiple analyses referencing the same GCS directory share the same Lustre copy.
    Reference counting ensures data is only deleted when no analyses need it.

    Usage in get_command():
        lustre = LustreInputs(analysis, settings)
        lustre.add("/mnt/gcsfuse/data/sample1/file.bam")
        local_path = lustre.get("/mnt/gcsfuse/data/sample1/file.bam")
        # returns /scratch/shared_inputs/a1b2c3.../file.bam
        self._lustre_inputs = lustre  # store for write_command_script

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

        # Maps GCS parent directory -> md5 hash key
        self._gcs_dirs = {}

        # Maps original file path -> (gcs_parent_dir, filename)
        self._file_to_info = {}

        # Get config for path conversion
        gcp_config = get_gcp_config()
        self.lustre_mount = gcp_config.get("lustre_mount_path", "/scratch")
        self.gcsfuse_mount = gcp_config.get("gcsfuse_mount_path")  # e.g., "/mnt/gcsfuse"
        self.gcs_input_uri = gcp_config.get("gcs_input_uri")  # e.g., "gs://input-bucket"
        self.gcsfuse_output_mount = gcp_config.get("gcsfuse_output_mount_path")  # e.g., "/isabl/data"
        self.gcs_base_uri = gcp_config.get("gcs_base_uri")  # e.g., "gs://output-bucket"
        self.import_enabled = gcp_config.get("lustre_import_enabled", False)
        self.shared_base = gcp_config.get("shared_inputs_path", "/shared_inputs")

    def _gcsfuse_to_gcs_uri(self, gcsfuse_path):
        """Convert gcsfuse mount path to GCS URI.

        Supports both input data paths (from gcsfuse_mount_path) and analysis output
        paths (from gcsfuse_output_mount_path).

        Arguments:
            gcsfuse_path (str): Local gcsfuse path (e.g., /mnt/gcsfuse/data/file.fastq
                or /isabl/data/analyses/00/01/123/file.bam).

        Returns:
            str: GCS URI (e.g., gs://input-bucket/data/file.fastq or
                gs://output-bucket/analyses/00/01/123/file.bam).

        Raises:
            GCPLustreImportError: If path conversion fails.

        Example:
            /mnt/gcsfuse/data/file.fastq -> gs://input-bucket/data/file.fastq
            /isabl/data/analyses/00/01/123/file.bam -> gs://output-bucket/analyses/00/01/123/file.bam
        """
        # If already a GCS URI, return as-is
        if gcsfuse_path.startswith("gs://"):
            return gcsfuse_path

        # Try input mount path first
        if self.gcsfuse_mount and self.gcs_input_uri:
            if gcsfuse_path.startswith(self.gcsfuse_mount):
                relative = gcsfuse_path[len(self.gcsfuse_mount):]
                if not relative.startswith("/"):
                    relative = "/" + relative
                return f"{self.gcs_input_uri.rstrip('/')}{relative}"

        # Try output mount path (for analysis dependencies)
        if self.gcsfuse_output_mount and self.gcs_base_uri:
            if gcsfuse_path.startswith(self.gcsfuse_output_mount):
                relative = gcsfuse_path[len(self.gcsfuse_output_mount):]
                if not relative.startswith("/"):
                    relative = "/" + relative
                return f"{self.gcs_base_uri.rstrip('/')}{relative}"

        # Build helpful error message
        configured_mounts = []
        if self.gcsfuse_mount:
            configured_mounts.append(f"input: {self.gcsfuse_mount}")
        if self.gcsfuse_output_mount:
            configured_mounts.append(f"output: {self.gcsfuse_output_mount}")

        if not configured_mounts:
            raise GCPLustreImportError(
                "No gcsfuse mount paths configured. Set gcsfuse_mount_path and/or "
                "gcsfuse_output_mount_path in GCP_CONFIGURATION"
            )

        raise GCPLustreImportError(
            f"Path doesn't start with any configured gcsfuse mount "
            f"({', '.join(configured_mounts)}): {gcsfuse_path}"
        )

    def _get_shared_dir_key(self, gcs_dir):
        """Generate deterministic hash for a GCS directory.

        Uses full MD5 hex digest for collision resistance.

        Arguments:
            gcs_dir (str): GCS directory URI (e.g., gs://bucket/exp1/raw_data/).

        Returns:
            str: MD5 hex digest of the GCS directory URI.
        """
        return hashlib.md5(gcs_dir.encode()).hexdigest()

    def _get_shared_lustre_path(self, gcs_dir):
        """Get the shared Lustre path for a GCS directory.

        Arguments:
            gcs_dir (str): GCS directory URI.

        Returns:
            str: Full Lustre path (e.g., /scratch/shared_inputs/a1b2c3.../).
        """
        hash_key = self._get_shared_dir_key(gcs_dir)
        return f"{self.lustre_mount}{self.shared_base}/{hash_key}"

    def add(self, path):
        """Register a file path to be imported to Lustre.

        Files are grouped by their parent directory since gcloud lustre import-data
        only supports directory-level imports. Files are placed in a shared location
        so multiple analyses can share the same import.

        Arguments:
            path (str): Local gcsfuse path (e.g., /mnt/gcsfuse/data/sample1/file.bam)
                        or GCS URI (e.g., gs://bucket/data/sample1/file.bam).

        Returns:
            str: The Lustre-local path where the file will be available.
        """
        if not self.import_enabled:
            # When import is disabled, return the original path
            return path

        if path in self._file_to_info:
            # Already registered, compute and return existing Lustre path
            gcs_parent_dir, filename = self._file_to_info[path]
            shared_path = self._get_shared_lustre_path(gcs_parent_dir)
            return f"{shared_path}/{filename}"

        # Convert to GCS URI
        gcs_uri = self._gcsfuse_to_gcs_uri(path)

        # Extract parent directory and filename
        filename = os.path.basename(gcs_uri.rstrip("/"))
        gcs_parent_dir = os.path.dirname(gcs_uri)
        if not gcs_parent_dir.endswith("/"):
            gcs_parent_dir += "/"

        # Register the directory if not already tracked
        if gcs_parent_dir not in self._gcs_dirs:
            self._gcs_dirs[gcs_parent_dir] = self._get_shared_dir_key(gcs_parent_dir)

        # Track the file
        self._file_to_info[path] = (gcs_parent_dir, filename)

        shared_path = self._get_shared_lustre_path(gcs_parent_dir)
        return f"{shared_path}/{filename}"

    def get(self, path):
        """Get the Lustre-local path for a registered file.

        Arguments:
            path (str): The original path passed to add().

        Returns:
            str: Full Lustre path (e.g., /scratch/shared_inputs/a1b2c3.../file.bam).

        Raises:
            ValueError: If path was not previously registered with add().
        """
        if not self.import_enabled:
            # When import is disabled, return the original path
            return path

        if path not in self._file_to_info:
            raise ValueError(f"Path not registered: {path}. Call add() first.")

        gcs_parent_dir, filename = self._file_to_info[path]
        shared_path = self._get_shared_lustre_path(gcs_parent_dir)
        return f"{shared_path}/{filename}"

    def get_import_specs(self):
        """Get the list of directory import specifications.

        Returns:
            list: List of (gcs_dir, lustre_dir) tuples for directory imports.
                  lustre_dir is relative to the lustre mount point.
        """
        specs = []
        for gcs_dir, hash_key in self._gcs_dirs.items():
            lustre_path = f"{self.shared_base}/{hash_key}/"
            specs.append((gcs_dir, lustre_path))
        return specs

    def get_import_command(self, strategy="lustre-import"):
        """Get CLI command to run shared imports (for embedding in script).

        Arguments:
            strategy (str): Import strategy - "lustre-import" or "rsync".

        Returns:
            str: CLI command string, or empty string if no imports needed.
        """
        if not self._gcs_dirs or not self.import_enabled:
            return ""

        specs = self.get_import_specs()
        specs_json = json.dumps(specs)
        # Use single quotes around JSON and escape any single quotes within
        escaped_json = specs_json.replace("'", "'\"'\"'")
        analysis_pk = self.analysis["pk"]
        return (
            f"isabl lustre-shared-import "
            f"--analysis-pk {analysis_pk} "
            f"--strategy {strategy} "
            f"--specs '{escaped_json}'"
        )

    def get_cleanup_command(self):
        """Get CLI command to release shared input refs (for embedding in script).

        Returns:
            str: CLI command string, or empty string if no imports needed.
        """
        if not self._gcs_dirs or not self.import_enabled:
            return ""

        specs = self.get_import_specs()
        specs_json = json.dumps(specs)
        escaped_json = specs_json.replace("'", "'\"'\"'")
        analysis_pk = self.analysis["pk"]
        return (
            f"isabl lustre-shared-cleanup "
            f"--analysis-pk {analysis_pk} "
            f"--specs '{escaped_json}'"
        )

    def get_directories(self):
        """Get the list of unique GCS directories to import.

        Returns:
            list: List of GCS directory URIs.
        """
        return list(self._gcs_dirs.keys())

    def __len__(self):
        """Return the number of registered files."""
        return len(self._file_to_info)

    def __repr__(self):
        """String representation of LustreInputs."""
        return (
            f"LustreInputs(analysis_pk={self.analysis['pk']}, "
            f"files={len(self._file_to_info)}, "
            f"directories={len(self._gcs_dirs)}, "
            f"enabled={self.import_enabled})"
        )
