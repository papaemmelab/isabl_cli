"""Tests for GCP Lustre import/export functionality."""

import json
import subprocess
from unittest.mock import MagicMock, patch

import pytest

from isabl_cli import gcp_lustre
from isabl_cli.lustre_inputs import LustreInputs
from isabl_cli.settings import _DEFAULTS


class TestGetGcsPathFromAnalysis:
    """Tests for get_gcs_path_from_analysis function."""

    def test_extracts_relative_path_from_storage_url(self):
        """Test that relative path is extracted from analysis storage_url."""
        mock_analysis = {"storage_url": "/datalake/analyses/00/01/123"}

        with patch("isabl_cli.api.get_instance", return_value=mock_analysis):
            result = gcp_lustre.get_gcs_path_from_analysis(
                analysis_pk=123,
                gcs_base_uri="gs://my-bucket",
                base_storage_directory="/datalake",
            )
            assert result == "gs://my-bucket/analyses/00/01/123/"

    def test_handles_trailing_slash_in_gcs_uri(self):
        """Test handling of trailing slashes in GCS URI."""
        mock_analysis = {"storage_url": "/datalake/analyses/00/01/123"}

        with patch("isabl_cli.api.get_instance", return_value=mock_analysis):
            result = gcp_lustre.get_gcs_path_from_analysis(
                analysis_pk=123,
                gcs_base_uri="gs://my-bucket/",
                base_storage_directory="/datalake",
            )
            assert result == "gs://my-bucket/analyses/00/01/123/"

    def test_handles_storage_url_with_trailing_slash(self):
        """Test handling of storage_url with trailing slash."""
        mock_analysis = {"storage_url": "/datalake/analyses/00/01/123/"}

        with patch("isabl_cli.api.get_instance", return_value=mock_analysis):
            result = gcp_lustre.get_gcs_path_from_analysis(
                analysis_pk=123,
                gcs_base_uri="gs://my-bucket",
                base_storage_directory="/datalake",
            )
            assert result == "gs://my-bucket/analyses/00/01/123/"

    def test_handles_base_storage_with_trailing_slash(self):
        """Test handling of base_storage_directory with trailing slash."""
        mock_analysis = {"storage_url": "/datalake/analyses/00/01/123"}

        with patch("isabl_cli.api.get_instance", return_value=mock_analysis):
            result = gcp_lustre.get_gcs_path_from_analysis(
                analysis_pk=123,
                gcs_base_uri="gs://my-bucket",
                base_storage_directory="/datalake/",
            )
            assert result == "gs://my-bucket/analyses/00/01/123/"

    def test_handles_storage_url_not_starting_with_base(self):
        """Test when storage_url doesn't start with base_storage_directory."""
        mock_analysis = {"storage_url": "/other/path/analyses/123"}

        with patch("isabl_cli.api.get_instance", return_value=mock_analysis):
            result = gcp_lustre.get_gcs_path_from_analysis(
                analysis_pk=123,
                gcs_base_uri="gs://my-bucket",
                base_storage_directory="/datalake",
            )
            assert result == "gs://my-bucket/other/path/analyses/123/"

    def test_raises_on_api_error(self):
        """Test that error is raised when API call fails."""
        with patch("isabl_cli.api.get_instance", side_effect=Exception("API error")):
            with pytest.raises(gcp_lustre.GCPLustreExportError) as exc_info:
                gcp_lustre.get_gcs_path_from_analysis(
                    analysis_pk=123,
                    gcs_base_uri="gs://my-bucket",
                    base_storage_directory="/datalake",
                )
            assert "Failed to fetch analysis" in str(exc_info.value)

    def test_raises_on_missing_storage_url(self):
        """Test that error is raised when analysis has no storage_url."""
        mock_analysis = {"storage_url": None}

        with patch("isabl_cli.api.get_instance", return_value=mock_analysis):
            with pytest.raises(gcp_lustre.GCPLustreExportError) as exc_info:
                gcp_lustre.get_gcs_path_from_analysis(
                    analysis_pk=123,
                    gcs_base_uri="gs://my-bucket",
                    base_storage_directory="/datalake",
                )
            assert "has no storage_url" in str(exc_info.value)


class TestNormalizeLustrePath:
    """Tests for normalize_lustre_path function."""

    def test_normalizes_path_with_leading_slash(self):
        """Test that path with leading slash is normalized."""
        result = gcp_lustre.normalize_lustre_path(
            lustre_path="/havasove/isabl/3",
        )
        assert result == "/havasove/isabl/3/"

    def test_handles_path_without_leading_slash(self):
        """Test handling of path without leading slash."""
        result = gcp_lustre.normalize_lustre_path(
            lustre_path="havasove/isabl/3",
        )
        assert result == "/havasove/isabl/3/"

    def test_ensures_trailing_slash(self):
        """Test that trailing slash is added."""
        result = gcp_lustre.normalize_lustre_path(
            lustre_path="/scratch/output",
        )
        assert result.endswith("/")

    def test_preserves_trailing_slash_if_present(self):
        """Test that existing trailing slash is preserved."""
        result = gcp_lustre.normalize_lustre_path(
            lustre_path="/havasove/isabl/3/",
        )
        assert result == "/havasove/isabl/3/"


class TestInitiateExport:
    """Tests for initiate_export function."""

    @pytest.fixture
    def gcp_config(self):
        """Sample GCP configuration."""
        return {
            "lustre_instance": "test-instance",
            "lustre_location": "us-east4-a",
            "lustre_project": "test-project",
        }

    def test_builds_correct_command(self, gcp_config):
        """Test that the correct gcloud command is built."""
        mock_result = MagicMock()
        mock_result.stdout = json.dumps({"name": "operations/test-op-123"})

        with patch("subprocess.run", return_value=mock_result) as mock_run:
            result = gcp_lustre.initiate_export(
                gcp_config, "/analyses/123/", "gs://bucket/analyses/123/"
            )

            mock_run.assert_called_once()
            cmd = mock_run.call_args[0][0]
            assert cmd[0] == "gcloud"
            assert cmd[1] == "lustre"
            assert cmd[2] == "instances"
            assert cmd[3] == "export-data"
            assert cmd[4] == "test-instance"
            assert "--async" in cmd
            assert "--format=json" in cmd
            assert "--location=us-east4-a" in cmd
            assert "--project=test-project" in cmd
            assert "--gcs-path-uri=gs://bucket/analyses/123/" in cmd
            assert "--lustre-path=/analyses/123/" in cmd

    def test_returns_operation_name(self, gcp_config):
        """Test that operation name is extracted from output."""
        mock_result = MagicMock()
        mock_result.stdout = json.dumps({"name": "operations/my-export-op"})

        with patch("subprocess.run", return_value=mock_result):
            result = gcp_lustre.initiate_export(
                gcp_config, "/analyses/123/", "gs://bucket/analyses/123/"
            )
            assert result == "operations/my-export-op"

    def test_raises_on_subprocess_error(self, gcp_config):
        """Test that GCPLustreExportError is raised on subprocess failure."""
        with patch(
            "subprocess.run",
            side_effect=subprocess.CalledProcessError(1, "gcloud", stderr="error"),
        ):
            with pytest.raises(gcp_lustre.GCPLustreExportError) as exc_info:
                gcp_lustre.initiate_export(
                    gcp_config, "/analyses/123/", "gs://bucket/analyses/123/"
                )
            assert "Failed to initiate" in str(exc_info.value)

    def test_raises_on_missing_operation_name(self, gcp_config):
        """Test that error is raised when operation name is missing."""
        mock_result = MagicMock()
        mock_result.stdout = json.dumps({"other": "data"})

        with patch("subprocess.run", return_value=mock_result):
            with pytest.raises(gcp_lustre.GCPLustreExportError) as exc_info:
                gcp_lustre.initiate_export(
                    gcp_config, "/analyses/123/", "gs://bucket/analyses/123/"
                )
            assert "Could not extract operation name" in str(exc_info.value)


class TestCheckExportStatus:
    """Tests for check_export_status function."""

    @pytest.fixture
    def gcp_config(self):
        """Sample GCP configuration."""
        return {
            "lustre_location": "us-east4-a",
            "lustre_project": "test-project",
        }

    def test_returns_done_true_on_completion(self, gcp_config):
        """Test that done=True is returned when operation completes."""
        mock_result = MagicMock()
        mock_result.stdout = json.dumps({"done": True})

        with patch("subprocess.run", return_value=mock_result):
            done, error = gcp_lustre.check_export_status(gcp_config, "op-123")
            assert done is True
            assert error is None

    def test_returns_done_false_when_in_progress(self, gcp_config):
        """Test that done=False is returned when still in progress."""
        mock_result = MagicMock()
        mock_result.stdout = json.dumps({"done": False})

        with patch("subprocess.run", return_value=mock_result):
            done, error = gcp_lustre.check_export_status(gcp_config, "op-123")
            assert done is False
            assert error is None

    def test_returns_error_on_failure(self, gcp_config):
        """Test that error message is returned when operation fails."""
        mock_result = MagicMock()
        mock_result.stdout = json.dumps({
            "done": True,
            "error": {"message": "Export failed: permission denied"}
        })

        with patch("subprocess.run", return_value=mock_result):
            done, error = gcp_lustre.check_export_status(gcp_config, "op-123")
            assert done is True
            assert error == "Export failed: permission denied"

    def test_builds_correct_command(self, gcp_config):
        """Test that the correct gcloud command is built."""
        mock_result = MagicMock()
        mock_result.stdout = json.dumps({"done": False})

        with patch("subprocess.run", return_value=mock_result) as mock_run:
            gcp_lustre.check_export_status(gcp_config, "operations/test-op")

            cmd = mock_run.call_args[0][0]
            assert cmd[0] == "gcloud"
            assert cmd[1] == "lustre"
            assert cmd[2] == "operations"
            assert cmd[3] == "describe"
            assert cmd[4] == "operations/test-op"
            assert "--format=json" in cmd


class TestWaitForExport:
    """Tests for wait_for_export function."""

    @pytest.fixture
    def gcp_config(self):
        """Sample GCP configuration with short intervals for testing."""
        return {
            "lustre_location": "us-east4-a",
            "lustre_project": "test-project",
            "lustre_poll_interval": 0.01,  # Very short for testing
            "lustre_max_poll_attempts": 3,
        }

    def test_returns_on_success(self, gcp_config):
        """Test that function returns when export succeeds."""
        with patch(
            "isabl_cli.gcp_lustre.check_export_status", return_value=(True, None)
        ):
            with patch("time.sleep"):
                # Should not raise
                gcp_lustre.wait_for_export(gcp_config, "op-123")

    def test_raises_on_export_error(self, gcp_config):
        """Test that error is raised when export fails."""
        with patch(
            "isabl_cli.gcp_lustre.check_export_status",
            return_value=(True, "Permission denied"),
        ):
            with patch("time.sleep"):
                with pytest.raises(gcp_lustre.GCPLustreExportError) as exc_info:
                    gcp_lustre.wait_for_export(gcp_config, "op-123")
                assert "Permission denied" in str(exc_info.value)

    def test_raises_on_timeout(self, gcp_config):
        """Test that error is raised when polling times out."""
        with patch(
            "isabl_cli.gcp_lustre.check_export_status", return_value=(False, None)
        ):
            with patch("time.sleep"):
                with pytest.raises(gcp_lustre.GCPLustreExportError) as exc_info:
                    gcp_lustre.wait_for_export(gcp_config, "op-123")
                assert "timed out" in str(exc_info.value)

    def test_retries_on_status_check_error(self, gcp_config):
        """Test that polling continues when status check fails."""
        call_count = [0]

        def mock_check(*args):
            call_count[0] += 1
            if call_count[0] < 2:
                raise gcp_lustre.GCPLustreExportError("Temporary error")
            return (True, None)

        with patch("isabl_cli.gcp_lustre.check_export_status", side_effect=mock_check):
            with patch("time.sleep"):
                gcp_lustre.wait_for_export(gcp_config, "op-123")
                assert call_count[0] == 2


class TestRunExport:
    """Tests for run_export function."""

    @pytest.fixture
    def full_gcp_config(self):
        """Complete GCP configuration."""
        return {
            "lustre_export_enabled": True,
            "lustre_instance": "test-instance",
            "lustre_location": "us-east4-a",
            "lustre_project": "test-project",
            "lustre_mount_path": "/lustre",
            "gcs_base_uri": "gs://my-bucket",
            "lustre_poll_interval": 0.01,
            "lustre_max_poll_attempts": 3,
            "lustre_delete_after_export": False,
        }

    def test_raises_when_not_enabled(self):
        """Test that error is raised when export is not enabled."""
        with patch(
            "isabl_cli.gcp_lustre.get_gcp_config",
            return_value={"lustre_export_enabled": False},
        ):
            with pytest.raises(gcp_lustre.GCPLustreExportError) as exc_info:
                gcp_lustre.run_export("/scratch/output", 123)
            assert "not enabled" in str(exc_info.value)

    def test_raises_when_missing_settings(self):
        """Test that error is raised when required settings are missing."""
        with patch(
            "isabl_cli.gcp_lustre.get_gcp_config",
            return_value={"lustre_export_enabled": True, "lustre_instance": "test"},
        ):
            with pytest.raises(gcp_lustre.GCPLustreExportError) as exc_info:
                gcp_lustre.run_export("/scratch/output", 123)
            assert "Missing required" in str(exc_info.value)

    def test_full_export_workflow(self, full_gcp_config):
        """Test the complete export workflow."""
        with patch("isabl_cli.gcp_lustre.get_gcp_config", return_value=full_gcp_config):
            with patch(
                "isabl_cli.gcp_lustre.get_gcs_path_from_analysis",
                return_value="gs://my-bucket/analyses/00/01/123/",
            ):
                with patch(
                    "isabl_cli.gcp_lustre.initiate_export", return_value="op-123"
                ) as mock_init:
                    with patch("isabl_cli.gcp_lustre.wait_for_export") as mock_wait:
                        # Pass path relative to Lustre filesystem root
                        gcp_lustre.run_export("/havasove/isabl/3", 123)

                        mock_init.assert_called_once()
                        # Check that lustre path is normalized (relative path, not absolute)
                        call_args = mock_init.call_args
                        assert call_args[0][1] == "/havasove/isabl/3/"
                        assert call_args[0][2] == "gs://my-bucket/analyses/00/01/123/"
                        mock_wait.assert_called_once_with(full_gcp_config, "op-123")

    def test_deletes_scratch_when_configured(self, full_gcp_config, tmp_path):
        """Test that scratch is deleted when configured."""
        full_gcp_config["lustre_delete_after_export"] = True
        full_gcp_config["lustre_mount_path"] = str(tmp_path)

        # Create a temp directory to delete
        scratch_dir = tmp_path / "havasove" / "isabl" / "3"
        scratch_dir.mkdir(parents=True)
        (scratch_dir / "test.txt").write_text("test")

        with patch("isabl_cli.gcp_lustre.get_gcp_config", return_value=full_gcp_config):
            with patch(
                "isabl_cli.gcp_lustre.get_gcs_path_from_analysis",
                return_value="gs://my-bucket/analyses/00/01/123/",
            ):
                with patch("isabl_cli.gcp_lustre.initiate_export", return_value="op-123"):
                    with patch("isabl_cli.gcp_lustre.wait_for_export"):
                        # Pass path relative to Lustre filesystem (without mount prefix)
                        gcp_lustre.run_export("/havasove/isabl/3", 123)

        assert not scratch_dir.exists()

    def test_respects_delete_after_override(self, full_gcp_config, tmp_path):
        """Test that delete_after parameter overrides config."""
        full_gcp_config["lustre_delete_after_export"] = True
        full_gcp_config["lustre_mount_path"] = str(tmp_path)

        scratch_dir = tmp_path / "havasove" / "isabl" / "3"
        scratch_dir.mkdir(parents=True)
        (scratch_dir / "test.txt").write_text("test")

        with patch("isabl_cli.gcp_lustre.get_gcp_config", return_value=full_gcp_config):
            with patch(
                "isabl_cli.gcp_lustre.get_gcs_path_from_analysis",
                return_value="gs://my-bucket/analyses/00/01/123/",
            ):
                with patch("isabl_cli.gcp_lustre.initiate_export", return_value="op-123"):
                    with patch("isabl_cli.gcp_lustre.wait_for_export"):
                        # Override to not delete
                        gcp_lustre.run_export("/havasove/isabl/3", 123, delete_after=False)

        assert scratch_dir.exists()


class TestGetExportCommandForScript:
    """Tests for get_export_command_for_script function."""

    @pytest.fixture
    def analysis(self):
        """Sample analysis dict."""
        return {"pk": 123, "storage_url": "/datalake/analyses/00/01/123"}

    def test_returns_empty_when_disabled(self, analysis, monkeypatch):
        """Test returns empty string when export is disabled."""
        monkeypatch.setattr(
            "isabl_cli.gcp_lustre.get_gcp_config",
            lambda: {"lustre_export_enabled": False},
        )

        result = gcp_lustre.get_export_command_for_script(analysis, "/scratch/output")
        assert result == ""

    def test_returns_empty_when_config_missing(self, analysis, monkeypatch):
        """Test returns empty string when GCP config is missing."""
        monkeypatch.setattr("isabl_cli.gcp_lustre.get_gcp_config", lambda: {})

        result = gcp_lustre.get_export_command_for_script(analysis, "/scratch/output")
        assert result == ""

    def test_returns_empty_when_required_settings_missing(self, analysis, monkeypatch):
        """Test returns empty string when required settings are missing."""
        monkeypatch.setattr(
            "isabl_cli.gcp_lustre.get_gcp_config",
            lambda: {
                "lustre_export_enabled": True,
                "lustre_instance": "test",
                # Missing other required settings
            },
        )

        result = gcp_lustre.get_export_command_for_script(analysis, "/scratch/output")
        assert result == ""

    def test_returns_cli_command_when_properly_configured(self, analysis, monkeypatch):
        """Test returns CLI command when properly configured."""
        monkeypatch.setattr(
            "isabl_cli.gcp_lustre.get_gcp_config",
            lambda: {
                "lustre_export_enabled": True,
                "lustre_instance": "test-instance",
                "lustre_location": "us-east4-a",
                "lustre_project": "test-project",
                "gcs_base_uri": "gs://my-bucket",
                "lustre_delete_after_export": True,
            },
        )

        result = gcp_lustre.get_export_command_for_script(analysis, "/scratch/output")
        assert result == "isabl lustre-export --lustre-path /scratch/output --analysis-pk 123 --delete-after"

    def test_includes_delete_after_flag_when_true(self, analysis, monkeypatch):
        """Test that --delete-after flag is included when config is True."""
        monkeypatch.setattr(
            "isabl_cli.gcp_lustre.get_gcp_config",
            lambda: {
                "lustre_export_enabled": True,
                "lustre_instance": "test-instance",
                "lustre_location": "us-east4-a",
                "lustre_project": "test-project",
                "gcs_base_uri": "gs://my-bucket",
                "lustre_delete_after_export": True,
            },
        )

        result = gcp_lustre.get_export_command_for_script(analysis, "/scratch/output")
        assert "--delete-after" in result
        assert "--no-delete-after" not in result

    def test_includes_no_delete_after_flag_when_false(self, analysis, monkeypatch):
        """Test that --no-delete-after flag is included when config is False."""
        monkeypatch.setattr(
            "isabl_cli.gcp_lustre.get_gcp_config",
            lambda: {
                "lustre_export_enabled": True,
                "lustre_instance": "test-instance",
                "lustre_location": "us-east4-a",
                "lustre_project": "test-project",
                "gcs_base_uri": "gs://my-bucket",
                "lustre_delete_after_export": False,
            },
        )

        result = gcp_lustre.get_export_command_for_script(analysis, "/scratch/output")
        assert "--no-delete-after" in result

    def test_defaults_to_delete_after_when_not_in_config(self, analysis, monkeypatch):
        """Test that --delete-after is used when config key is missing (defaults to True)."""
        monkeypatch.setattr(
            "isabl_cli.gcp_lustre.get_gcp_config",
            lambda: {
                "lustre_export_enabled": True,
                "lustre_instance": "test-instance",
                "lustre_location": "us-east4-a",
                "lustre_project": "test-project",
                "gcs_base_uri": "gs://my-bucket",
                # lustre_delete_after_export not set, should default to True
            },
        )

        result = gcp_lustre.get_export_command_for_script(analysis, "/scratch/output")
        assert "--delete-after" in result
        assert "--no-delete-after" not in result


class TestGCPConfiguration:
    """Tests for GCP configuration in settings."""

    def test_default_config_has_export_disabled(self):
        """Test that default GCP config has export disabled."""
        gcp_config = _DEFAULTS.get("GCP_CONFIGURATION", {})
        assert gcp_config.get("lustre_export_enabled") is False

    def test_default_config_has_all_required_keys(self):
        """Test that default GCP config has all expected keys."""
        gcp_config = _DEFAULTS.get("GCP_CONFIGURATION", {})
        expected_keys = [
            "lustre_export_enabled",
            "lustre_instance",
            "lustre_location",
            "lustre_project",
            "lustre_mount_path",
            "gcs_base_uri",
            "lustre_poll_interval",
            "lustre_max_poll_attempts",
            "lustre_delete_after_export",
        ]
        for key in expected_keys:
            assert key in gcp_config, f"Missing key: {key}"

    def test_default_poll_interval(self):
        """Test default poll interval value."""
        gcp_config = _DEFAULTS.get("GCP_CONFIGURATION", {})
        assert gcp_config.get("lustre_poll_interval") == 30

    def test_default_max_poll_attempts(self):
        """Test default max poll attempts value."""
        gcp_config = _DEFAULTS.get("GCP_CONFIGURATION", {})
        assert gcp_config.get("lustre_max_poll_attempts") == 360

    def test_default_delete_after_export(self):
        """Test default delete after export value."""
        gcp_config = _DEFAULTS.get("GCP_CONFIGURATION", {})
        assert gcp_config.get("lustre_delete_after_export") is True

    def test_default_import_disabled(self):
        """Test that default GCP config has import disabled."""
        gcp_config = _DEFAULTS.get("GCP_CONFIGURATION", {})
        assert gcp_config.get("lustre_import_enabled") is False

    def test_default_config_has_import_keys(self):
        """Test that default GCP config has import-related keys."""
        gcp_config = _DEFAULTS.get("GCP_CONFIGURATION", {})
        assert "lustre_import_enabled" in gcp_config
        assert "gcs_input_uri" in gcp_config
        assert "gcsfuse_mount_path" in gcp_config


# ============================================================================
# IMPORT TESTS
# ============================================================================


class TestInitiateImport:
    """Tests for initiate_import function."""

    @pytest.fixture
    def gcp_config(self):
        """Sample GCP configuration."""
        return {
            "lustre_instance": "test-instance",
            "lustre_location": "us-east4-a",
            "lustre_project": "test-project",
        }

    def test_builds_correct_command(self, gcp_config):
        """Test that the correct gcloud command is built."""
        mock_result = MagicMock()
        mock_result.stdout = json.dumps({"name": "operations/test-op-123"})

        with patch("subprocess.run", return_value=mock_result) as mock_run:
            result = gcp_lustre.initiate_import(
                gcp_config, "gs://bucket/file.fastq", "/123/inputs/file.fastq"
            )

            mock_run.assert_called_once()
            cmd = mock_run.call_args[0][0]
            assert cmd[0] == "gcloud"
            assert cmd[1] == "lustre"
            assert cmd[2] == "instances"
            assert cmd[3] == "import-data"
            assert cmd[4] == "test-instance"
            assert "--async" in cmd
            assert "--format=json" in cmd
            assert "--location=us-east4-a" in cmd
            assert "--project=test-project" in cmd
            assert "--gcs-path-uri=gs://bucket/file.fastq" in cmd
            assert "--lustre-path=/123/inputs/file.fastq" in cmd

    def test_returns_operation_name(self, gcp_config):
        """Test that operation name is extracted from output."""
        mock_result = MagicMock()
        mock_result.stdout = json.dumps({"name": "operations/my-import-op"})

        with patch("subprocess.run", return_value=mock_result):
            result = gcp_lustre.initiate_import(
                gcp_config, "gs://bucket/file.fastq", "/123/inputs/file.fastq"
            )
            assert result == "operations/my-import-op"

    def test_raises_on_subprocess_error(self, gcp_config):
        """Test that GCPLustreImportError is raised on subprocess failure."""
        with patch(
            "subprocess.run",
            side_effect=subprocess.CalledProcessError(1, "gcloud", stderr="error"),
        ):
            with pytest.raises(gcp_lustre.GCPLustreImportError) as exc_info:
                gcp_lustre.initiate_import(
                    gcp_config, "gs://bucket/file.fastq", "/123/inputs/file.fastq"
                )
            assert "Failed to initiate" in str(exc_info.value)

    def test_raises_on_missing_operation_name(self, gcp_config):
        """Test that error is raised when operation name is missing."""
        mock_result = MagicMock()
        mock_result.stdout = json.dumps({"other": "data"})

        with patch("subprocess.run", return_value=mock_result):
            with pytest.raises(gcp_lustre.GCPLustreImportError) as exc_info:
                gcp_lustre.initiate_import(
                    gcp_config, "gs://bucket/file.fastq", "/123/inputs/file.fastq"
                )
            assert "Could not extract operation name" in str(exc_info.value)


class TestCheckImportStatus:
    """Tests for check_import_status function."""

    @pytest.fixture
    def gcp_config(self):
        """Sample GCP configuration."""
        return {
            "lustre_location": "us-east4-a",
            "lustre_project": "test-project",
        }

    def test_returns_done_true_on_completion(self, gcp_config):
        """Test that done=True is returned when operation completes."""
        mock_result = MagicMock()
        mock_result.stdout = json.dumps({"done": True})

        with patch("subprocess.run", return_value=mock_result):
            done, error = gcp_lustre.check_import_status(gcp_config, "op-123")
            assert done is True
            assert error is None

    def test_returns_done_false_when_in_progress(self, gcp_config):
        """Test that done=False is returned when still in progress."""
        mock_result = MagicMock()
        mock_result.stdout = json.dumps({"done": False})

        with patch("subprocess.run", return_value=mock_result):
            done, error = gcp_lustre.check_import_status(gcp_config, "op-123")
            assert done is False
            assert error is None

    def test_returns_error_on_failure(self, gcp_config):
        """Test that error message is returned when operation fails."""
        mock_result = MagicMock()
        mock_result.stdout = json.dumps({
            "done": True,
            "error": {"message": "Import failed: permission denied"}
        })

        with patch("subprocess.run", return_value=mock_result):
            done, error = gcp_lustre.check_import_status(gcp_config, "op-123")
            assert done is True
            assert error == "Import failed: permission denied"


class TestWaitForImports:
    """Tests for wait_for_imports function."""

    @pytest.fixture
    def gcp_config(self):
        """Sample GCP configuration with short intervals for testing."""
        return {
            "lustre_location": "us-east4-a",
            "lustre_project": "test-project",
            "lustre_poll_interval": 0.01,  # Very short for testing
            "lustre_max_poll_attempts": 3,
        }

    def test_returns_on_all_success(self, gcp_config):
        """Test that function returns when all imports succeed."""
        with patch(
            "isabl_cli.gcp_lustre.check_import_status", return_value=(True, None)
        ):
            with patch("time.sleep"):
                # Should not raise
                gcp_lustre.wait_for_imports(gcp_config, ["op-1", "op-2"])

    def test_raises_on_import_error(self, gcp_config):
        """Test that error is raised when any import fails."""
        with patch(
            "isabl_cli.gcp_lustre.check_import_status",
            return_value=(True, "Permission denied"),
        ):
            with patch("time.sleep"):
                with pytest.raises(gcp_lustre.GCPLustreImportError) as exc_info:
                    gcp_lustre.wait_for_imports(gcp_config, ["op-123"])
                assert "Permission denied" in str(exc_info.value)

    def test_raises_on_timeout(self, gcp_config):
        """Test that error is raised when polling times out."""
        with patch(
            "isabl_cli.gcp_lustre.check_import_status", return_value=(False, None)
        ):
            with patch("time.sleep"):
                with pytest.raises(gcp_lustre.GCPLustreImportError) as exc_info:
                    gcp_lustre.wait_for_imports(gcp_config, ["op-123"])
                assert "timed out" in str(exc_info.value)

    def test_handles_empty_operations_list(self, gcp_config):
        """Test that empty operations list returns immediately."""
        gcp_lustre.wait_for_imports(gcp_config, [])  # Should not raise


class TestRunBatchImport:
    """Tests for run_batch_import function."""

    @pytest.fixture
    def full_gcp_config(self):
        """Complete GCP configuration."""
        return {
            "lustre_import_enabled": True,
            "lustre_instance": "test-instance",
            "lustre_location": "us-east4-a",
            "lustre_project": "test-project",
            "lustre_poll_interval": 0.01,
            "lustre_max_poll_attempts": 3,
        }

    def test_raises_when_not_enabled(self):
        """Test that error is raised when import is not enabled."""
        with patch(
            "isabl_cli.gcp_lustre.get_gcp_config",
            return_value={"lustre_import_enabled": False},
        ):
            with pytest.raises(gcp_lustre.GCPLustreImportError) as exc_info:
                gcp_lustre.run_batch_import([["gs://bucket/file", "/path/file"]])
            assert "not enabled" in str(exc_info.value)

    def test_raises_when_missing_settings(self):
        """Test that error is raised when required settings are missing."""
        with patch(
            "isabl_cli.gcp_lustre.get_gcp_config",
            return_value={"lustre_import_enabled": True, "lustre_instance": "test"},
        ):
            with pytest.raises(gcp_lustre.GCPLustreImportError) as exc_info:
                gcp_lustre.run_batch_import([["gs://bucket/file", "/path/file"]])
            assert "Missing required" in str(exc_info.value)

    def test_handles_empty_specs(self, full_gcp_config):
        """Test that empty import specs returns without error."""
        with patch("isabl_cli.gcp_lustre.get_gcp_config", return_value=full_gcp_config):
            # Should not raise
            gcp_lustre.run_batch_import([])

    def test_full_import_workflow(self, full_gcp_config):
        """Test the complete import workflow."""
        with patch("isabl_cli.gcp_lustre.get_gcp_config", return_value=full_gcp_config):
            with patch(
                "isabl_cli.gcp_lustre.initiate_import", return_value="op-123"
            ) as mock_init:
                with patch("isabl_cli.gcp_lustre.wait_for_imports") as mock_wait:
                    gcp_lustre.run_batch_import([
                        ["gs://bucket/file1.fastq", "/123/inputs/file1.fastq"],
                        ["gs://bucket/file2.fastq", "/123/inputs/file2.fastq"],
                    ])

                    assert mock_init.call_count == 2
                    mock_wait.assert_called_once()
                    # Check operations were passed to wait_for_imports
                    wait_call_args = mock_wait.call_args
                    assert len(wait_call_args[0][1]) == 2  # Two operations


class TestLustreInputs:
    """Tests for LustreInputs helper class."""

    @pytest.fixture
    def gcp_config(self):
        """Sample GCP configuration."""
        return {
            "lustre_import_enabled": True,
            "lustre_mount_path": "/scratch",
            "gcsfuse_mount_path": "/mnt/gcsfuse",
            "gcs_input_uri": "gs://input-bucket",
        }

    @pytest.fixture
    def analysis(self):
        """Sample analysis dict."""
        return {"pk": 123, "storage_url": "/datalake/analyses/00/01/123"}

    def test_add_and_get_path(self, gcp_config, analysis):
        """Test adding and getting a path."""
        with patch("isabl_cli.lustre_inputs.get_gcp_config", return_value=gcp_config):
            lustre = LustreInputs(analysis)
            lustre.add("/mnt/gcsfuse/data/file.bam")

            result = lustre.get("/mnt/gcsfuse/data/file.bam")
            assert result == "/scratch/123/inputs/file.bam"

    def test_add_returns_lustre_path(self, gcp_config, analysis):
        """Test that add() returns the Lustre path."""
        with patch("isabl_cli.lustre_inputs.get_gcp_config", return_value=gcp_config):
            lustre = LustreInputs(analysis)
            result = lustre.add("/mnt/gcsfuse/data/file.bam")

            assert result == "/scratch/123/inputs/file.bam"

    def test_gcsfuse_to_gcs_uri_conversion(self, gcp_config, analysis):
        """Test conversion of gcsfuse path to GCS URI."""
        with patch("isabl_cli.lustre_inputs.get_gcp_config", return_value=gcp_config):
            lustre = LustreInputs(analysis)
            gcs_uri = lustre._gcsfuse_to_gcs_uri("/mnt/gcsfuse/data/file.fastq")

            assert gcs_uri == "gs://input-bucket/data/file.fastq"

    def test_gcs_uri_passthrough(self, gcp_config, analysis):
        """Test that GCS URIs are passed through unchanged."""
        with patch("isabl_cli.lustre_inputs.get_gcp_config", return_value=gcp_config):
            lustre = LustreInputs(analysis)
            gcs_uri = lustre._gcsfuse_to_gcs_uri("gs://bucket/data/file.fastq")

            assert gcs_uri == "gs://bucket/data/file.fastq"

    def test_raises_on_invalid_path(self, gcp_config, analysis):
        """Test that error is raised for invalid paths."""
        with patch("isabl_cli.lustre_inputs.get_gcp_config", return_value=gcp_config):
            lustre = LustreInputs(analysis)

            with pytest.raises(gcp_lustre.GCPLustreImportError) as exc_info:
                lustre._gcsfuse_to_gcs_uri("/other/path/file.bam")
            assert "doesn't start with gcsfuse mount" in str(exc_info.value)

    def test_raises_on_missing_config(self, analysis):
        """Test that error is raised when gcsfuse config is missing."""
        with patch("isabl_cli.lustre_inputs.get_gcp_config", return_value={
            "lustre_import_enabled": True,
        }):
            lustre = LustreInputs(analysis)

            with pytest.raises(gcp_lustre.GCPLustreImportError) as exc_info:
                lustre._gcsfuse_to_gcs_uri("/mnt/gcsfuse/file.bam")
            assert "must be configured" in str(exc_info.value)

    def test_get_raises_on_unregistered_path(self, gcp_config, analysis):
        """Test that get() raises for unregistered paths."""
        with patch("isabl_cli.lustre_inputs.get_gcp_config", return_value=gcp_config):
            lustre = LustreInputs(analysis)

            with pytest.raises(ValueError) as exc_info:
                lustre.get("/mnt/gcsfuse/data/unknown.bam")
            assert "not registered" in str(exc_info.value)

    def test_handles_duplicate_filenames(self, gcp_config, analysis):
        """Test that duplicate filenames get unique paths."""
        with patch("isabl_cli.lustre_inputs.get_gcp_config", return_value=gcp_config):
            lustre = LustreInputs(analysis)

            path1 = lustre.add("/mnt/gcsfuse/dir1/file.bam")
            path2 = lustre.add("/mnt/gcsfuse/dir2/file.bam")

            # Both should have unique paths
            assert path1 != path2
            assert "file.bam" in path1
            assert "file" in path2  # Has hash suffix

    def test_get_import_specs(self, gcp_config, analysis):
        """Test get_import_specs returns correct format."""
        with patch("isabl_cli.lustre_inputs.get_gcp_config", return_value=gcp_config):
            lustre = LustreInputs(analysis)
            lustre.add("/mnt/gcsfuse/data/file.bam")

            specs = lustre.get_import_specs()
            assert len(specs) == 1
            assert specs[0][0] == "gs://input-bucket/data/file.bam"
            assert specs[0][1] == "/123/inputs/file.bam"

    def test_get_import_command(self, gcp_config, analysis):
        """Test get_import_command generates correct CLI command."""
        with patch("isabl_cli.lustre_inputs.get_gcp_config", return_value=gcp_config):
            lustre = LustreInputs(analysis)
            lustre.add("/mnt/gcsfuse/data/file.bam")

            cmd = lustre.get_import_command()
            assert "isabl lustre-import" in cmd
            assert "--specs" in cmd
            assert "gs://input-bucket/data/file.bam" in cmd

    def test_get_import_command_empty_when_no_files(self, gcp_config, analysis):
        """Test get_import_command returns empty string when no files."""
        with patch("isabl_cli.lustre_inputs.get_gcp_config", return_value=gcp_config):
            lustre = LustreInputs(analysis)

            cmd = lustre.get_import_command()
            assert cmd == ""

    def test_disabled_returns_original_paths(self, analysis):
        """Test that when import is disabled, original paths are returned."""
        with patch("isabl_cli.lustre_inputs.get_gcp_config", return_value={
            "lustre_import_enabled": False,
            "lustre_mount_path": "/scratch",
            "gcsfuse_mount_path": "/mnt/gcsfuse",
            "gcs_input_uri": "gs://input-bucket",
        }):
            lustre = LustreInputs(analysis)

            original_path = "/mnt/gcsfuse/data/file.bam"
            result = lustre.add(original_path)

            # When disabled, should return original path
            assert result == original_path
            assert lustre.get(original_path) == original_path

    def test_len(self, gcp_config, analysis):
        """Test __len__ returns number of registered files."""
        with patch("isabl_cli.lustre_inputs.get_gcp_config", return_value=gcp_config):
            lustre = LustreInputs(analysis)

            assert len(lustre) == 0
            lustre.add("/mnt/gcsfuse/file1.bam")
            assert len(lustre) == 1
            lustre.add("/mnt/gcsfuse/file2.bam")
            assert len(lustre) == 2

    def test_repr(self, gcp_config, analysis):
        """Test __repr__ returns useful information."""
        with patch("isabl_cli.lustre_inputs.get_gcp_config", return_value=gcp_config):
            lustre = LustreInputs(analysis)
            lustre.add("/mnt/gcsfuse/file.bam")

            repr_str = repr(lustre)
            assert "123" in repr_str  # analysis pk
            assert "1" in repr_str  # number of files
