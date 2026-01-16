"""LustreInputs helper class for apps to manage input imports from GCS to Lustre.

This module provides a convenient API for isabl applications to register input files
that need to be imported from GCS to Lustre scratch storage before pipeline execution.

Note: The gcloud lustre import-data command only supports directory-level imports,
so this class groups files by their parent directory and imports entire directories.
"""

import hashlib
import json
import os

from isabl_cli.gcp_lustre import get_gcp_config, GCPLustreImportError


class LustreInputs:
    """Helper for apps to manage input imports from GCS to Lustre.

    Input paths in isabl are stored as local gcsfuse mount paths (e.g., /mnt/gcsfuse/...).
    This class converts them to GCS URIs for import and provides Lustre-local paths.

    Since gcloud lustre import-data only supports directory imports, files are grouped
    by their parent directory. Each unique GCS directory is imported once, and files
    are accessed at their original relative paths within the imported directory.

    Usage in get_command():
        lustre = LustreInputs(analysis, settings)
        lustre.add("/mnt/gcsfuse/data/sample1/file.bam")  # gcsfuse path
        local_path = lustre.get("/mnt/gcsfuse/data/sample1/file.bam")
        # returns /scratch/123/inputs/data_sample1_a1b2/file.bam

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

        # Maps GCS parent directory -> Lustre subdirectory name (relative to input_dir)
        self._gcs_dirs = {}

        # Maps original file path -> (gcs_parent_dir, filename)
        self._file_to_info = {}

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

    def _generate_lustre_subdir_name(self, gcs_dir):
        """Generate a unique Lustre subdirectory name for a GCS directory.

        Uses the last 2 path components plus a 4-character hash suffix for uniqueness.

        Arguments:
            gcs_dir (str): GCS directory URI (e.g., gs://bucket/exp1/raw_data/).

        Returns:
            str: Subdirectory name (e.g., exp1_raw_data_a1b2).
        """
        # Remove gs://bucket prefix and trailing slash
        path = gcs_dir.split("://", 1)[-1]  # Remove gs://
        if "/" in path:
            path = path.split("/", 1)[-1]  # Remove bucket name
        path = path.rstrip("/")

        # Get last 2 path components
        parts = path.split("/")
        if len(parts) >= 2:
            name_parts = parts[-2:]
        else:
            name_parts = parts

        # Create base name from path components
        base_name = "_".join(name_parts)

        # Add hash suffix for uniqueness (in case different paths have same last components)
        hash_suffix = hashlib.md5(gcs_dir.encode()).hexdigest()[:4]

        return f"{base_name}_{hash_suffix}"

    def add(self, path):
        """Register a file path to be imported to Lustre.

        Files are grouped by their parent directory since gcloud lustre import-data
        only supports directory-level imports.

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
            lustre_subdir = self._gcs_dirs[gcs_parent_dir]
            return f"{self.lustre_mount}{self.input_dir}/{lustre_subdir}/{filename}"

        # Convert to GCS URI
        gcs_uri = self._gcsfuse_to_gcs_uri(path)

        # Extract parent directory and filename
        filename = os.path.basename(gcs_uri.rstrip("/"))
        gcs_parent_dir = os.path.dirname(gcs_uri)
        if not gcs_parent_dir.endswith("/"):
            gcs_parent_dir += "/"

        # Register the directory if not already tracked
        if gcs_parent_dir not in self._gcs_dirs:
            lustre_subdir = self._generate_lustre_subdir_name(gcs_parent_dir)
            self._gcs_dirs[gcs_parent_dir] = lustre_subdir
        else:
            lustre_subdir = self._gcs_dirs[gcs_parent_dir]

        # Track the file
        self._file_to_info[path] = (gcs_parent_dir, filename)

        return f"{self.lustre_mount}{self.input_dir}/{lustre_subdir}/{filename}"

    def get(self, path):
        """Get the Lustre-local path for a registered file.

        Arguments:
            path (str): The original path passed to add().

        Returns:
            str: Full Lustre path (e.g., /scratch/123/inputs/exp1_raw_data_a1b2/file.bam).

        Raises:
            ValueError: If path was not previously registered with add().
        """
        if not self.import_enabled:
            # When import is disabled, return the original path
            return path

        if path not in self._file_to_info:
            raise ValueError(f"Path not registered: {path}. Call add() first.")

        gcs_parent_dir, filename = self._file_to_info[path]
        lustre_subdir = self._gcs_dirs[gcs_parent_dir]
        return f"{self.lustre_mount}{self.input_dir}/{lustre_subdir}/{filename}"

    def get_import_specs(self):
        """Get the list of directory import specifications.

        Returns:
            list: List of (gcs_dir, lustre_dir) tuples for directory imports.
        """
        specs = []
        for gcs_dir, lustre_subdir in self._gcs_dirs.items():
            lustre_path = f"{self.input_dir}/{lustre_subdir}/"
            specs.append((gcs_dir, lustre_path))
        return specs

    def get_import_command(self):
        """Get CLI command to run imports (for embedding in script).

        Returns:
            str: CLI command string, or empty string if no imports needed.
        """
        if not self._gcs_dirs or not self.import_enabled:
            return ""

        specs = self.get_import_specs()
        specs_json = json.dumps(specs)
        # Use single quotes around JSON and escape any single quotes within
        escaped_json = specs_json.replace("'", "'\"'\"'")
        return f"isabl lustre-import --specs '{escaped_json}'"

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
