"""Tests for GCP Lustre export functionality."""

import json
import subprocess
from unittest.mock import MagicMock, patch

import pytest

from isabl_cli import gcp_lustre
from isabl_cli.settings import _DEFAULTS


class TestComputeExportPaths:
    """Tests for compute_export_paths function."""

    def test_strips_mount_path(self):
        """Test that mount path is properly stripped."""
        lustre_path, gcs_path = gcp_lustre.compute_export_paths(
            storage_url="/lustre/analyses/00/01/123",
            lustre_mount_path="/lustre",
            gcs_base_uri="gs://my-bucket",
        )
        assert lustre_path == "/analyses/00/01/123/"
        assert gcs_path == "gs://my-bucket/analyses/00/01/123/"

    def test_handles_trailing_slashes_in_gcs_uri(self):
        """Test handling of trailing slashes in GCS URI."""
        lustre_path, gcs_path = gcp_lustre.compute_export_paths(
            storage_url="/lustre/analyses/00/01/123",
            lustre_mount_path="/lustre",
            gcs_base_uri="gs://my-bucket/",
        )
        assert gcs_path == "gs://my-bucket/analyses/00/01/123/"

    def test_handles_storage_url_with_trailing_slash(self):
        """Test handling of storage_url with trailing slash."""
        lustre_path, gcs_path = gcp_lustre.compute_export_paths(
            storage_url="/lustre/analyses/00/01/123/",
            lustre_mount_path="/lustre",
            gcs_base_uri="gs://my-bucket",
        )
        assert lustre_path == "/analyses/00/01/123/"
        assert gcs_path == "gs://my-bucket/analyses/00/01/123/"

    def test_handles_storage_url_without_mount_prefix(self):
        """Test when storage_url doesn't start with mount path."""
        lustre_path, gcs_path = gcp_lustre.compute_export_paths(
            storage_url="/other/path/analyses/123",
            lustre_mount_path="/lustre",
            gcs_base_uri="gs://my-bucket",
        )
        assert lustre_path == "/other/path/analyses/123/"
        assert gcs_path == "gs://my-bucket/other/path/analyses/123/"

    def test_handles_deep_paths(self):
        """Test with deeply nested paths."""
        lustre_path, gcs_path = gcp_lustre.compute_export_paths(
            storage_url="/mnt/lustre/data/project/analyses/ab/cd/12345",
            lustre_mount_path="/mnt/lustre",
            gcs_base_uri="gs://gcs-output-bucket",
        )
        assert lustre_path == "/data/project/analyses/ab/cd/12345/"
        assert gcs_path == "gs://gcs-output-bucket/data/project/analyses/ab/cd/12345/"


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
                gcp_lustre.run_export("/lustre/analyses/123")
            assert "not enabled" in str(exc_info.value)

    def test_raises_when_missing_settings(self):
        """Test that error is raised when required settings are missing."""
        with patch(
            "isabl_cli.gcp_lustre.get_gcp_config",
            return_value={"lustre_export_enabled": True, "lustre_instance": "test"},
        ):
            with pytest.raises(gcp_lustre.GCPLustreExportError) as exc_info:
                gcp_lustre.run_export("/lustre/analyses/123")
            assert "Missing required" in str(exc_info.value)

    def test_full_export_workflow(self, full_gcp_config):
        """Test the complete export workflow."""
        with patch("isabl_cli.gcp_lustre.get_gcp_config", return_value=full_gcp_config):
            with patch(
                "isabl_cli.gcp_lustre.initiate_export", return_value="op-123"
            ) as mock_init:
                with patch("isabl_cli.gcp_lustre.wait_for_export") as mock_wait:
                    gcp_lustre.run_export("/lustre/analyses/123")

                    mock_init.assert_called_once()
                    mock_wait.assert_called_once_with(full_gcp_config, "op-123")

    def test_deletes_scratch_when_configured(self, full_gcp_config, tmp_path):
        """Test that scratch is deleted when configured."""
        full_gcp_config["lustre_delete_after_export"] = True
        full_gcp_config["lustre_mount_path"] = str(tmp_path)

        # Create a temp directory to delete
        scratch_dir = tmp_path / "analyses" / "123"
        scratch_dir.mkdir(parents=True)
        (scratch_dir / "test.txt").write_text("test")

        with patch("isabl_cli.gcp_lustre.get_gcp_config", return_value=full_gcp_config):
            with patch("isabl_cli.gcp_lustre.initiate_export", return_value="op-123"):
                with patch("isabl_cli.gcp_lustre.wait_for_export"):
                    gcp_lustre.run_export(str(scratch_dir))

        assert not scratch_dir.exists()

    def test_respects_delete_after_override(self, full_gcp_config, tmp_path):
        """Test that delete_after parameter overrides config."""
        full_gcp_config["lustre_delete_after_export"] = True
        full_gcp_config["lustre_mount_path"] = str(tmp_path)

        scratch_dir = tmp_path / "analyses" / "123"
        scratch_dir.mkdir(parents=True)
        (scratch_dir / "test.txt").write_text("test")

        with patch("isabl_cli.gcp_lustre.get_gcp_config", return_value=full_gcp_config):
            with patch("isabl_cli.gcp_lustre.initiate_export", return_value="op-123"):
                with patch("isabl_cli.gcp_lustre.wait_for_export"):
                    # Override to not delete
                    gcp_lustre.run_export(str(scratch_dir), delete_after=False)

        assert scratch_dir.exists()


class TestGetExportCommandForScript:
    """Tests for get_export_command_for_script function."""

    @pytest.fixture
    def analysis(self):
        """Sample analysis dict."""
        return {"pk": 123, "storage_url": "/lustre/analyses/00/01/123"}

    def test_returns_empty_when_disabled(self, analysis, monkeypatch):
        """Test returns empty string when export is disabled."""
        monkeypatch.setattr(
            "isabl_cli.gcp_lustre.get_gcp_config",
            lambda: {"lustre_export_enabled": False},
        )

        result = gcp_lustre.get_export_command_for_script(analysis)
        assert result == ""

    def test_returns_empty_when_config_missing(self, analysis, monkeypatch):
        """Test returns empty string when GCP config is missing."""
        monkeypatch.setattr("isabl_cli.gcp_lustre.get_gcp_config", lambda: {})

        result = gcp_lustre.get_export_command_for_script(analysis)
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

        result = gcp_lustre.get_export_command_for_script(analysis)
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
                "lustre_mount_path": "/lustre",
                "gcs_base_uri": "gs://my-bucket",
            },
        )

        result = gcp_lustre.get_export_command_for_script(analysis)
        assert result == "isabl lustre-export --storage-url /lustre/analyses/00/01/123"


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
