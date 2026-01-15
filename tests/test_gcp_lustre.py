"""Tests for GCP Lustre export functionality."""

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


class TestBuildExportScript:
    """Tests for build_export_script function."""

    @pytest.fixture
    def gcp_config(self):
        """Sample GCP configuration."""
        return {
            "lustre_export_enabled": True,
            "lustre_instance": "test-instance",
            "lustre_location": "us-east4-a",
            "lustre_project": "test-project",
            "lustre_mount_path": "/lustre",
            "gcs_base_uri": "gs://my-bucket",
            "lustre_poll_interval": 30,
            "lustre_max_poll_attempts": 10,
            "lustre_delete_after_export": False,
        }

    @pytest.fixture
    def analysis(self):
        """Sample analysis dict."""
        return {"pk": 123, "storage_url": "/lustre/analyses/00/01/123"}

    def test_script_contains_export_command(self, analysis, gcp_config):
        """Test that generated script contains gcloud export command."""
        script = gcp_lustre.build_export_script(analysis, gcp_config)

        assert "gcloud lustre instances export-data" in script
        assert "test-instance" in script
        assert "us-east4-a" in script
        assert "test-project" in script
        assert "--async" in script

    def test_script_contains_correct_paths(self, analysis, gcp_config):
        """Test that script uses correct lustre and GCS paths."""
        script = gcp_lustre.build_export_script(analysis, gcp_config)

        assert '--lustre-path="/analyses/00/01/123/"' in script
        assert '--gcs-path-uri="gs://my-bucket/analyses/00/01/123/"' in script

    def test_script_contains_polling_logic(self, analysis, gcp_config):
        """Test that generated script contains polling logic."""
        script = gcp_lustre.build_export_script(analysis, gcp_config)

        assert "gcloud lustre operations describe" in script
        assert "MAX_POLLS=10" in script
        assert "POLL_INTERVAL=30" in script
        assert "while" in script

    def test_script_includes_cleanup_when_enabled(self, analysis, gcp_config):
        """Test that script includes cleanup when lustre_delete_after_export is True."""
        gcp_config["lustre_delete_after_export"] = True
        script = gcp_lustre.build_export_script(analysis, gcp_config)

        assert "rm -rf" in script
        assert "Deleting scratch data" in script

    def test_script_excludes_cleanup_when_disabled(self, analysis, gcp_config):
        """Test that script excludes cleanup when lustre_delete_after_export is False."""
        gcp_config["lustre_delete_after_export"] = False
        script = gcp_lustre.build_export_script(analysis, gcp_config)

        assert "rm -rf" not in script

    def test_script_contains_error_handling(self, analysis, gcp_config):
        """Test that script has proper error handling."""
        script = gcp_lustre.build_export_script(analysis, gcp_config)

        assert "exit 1" in script
        assert "ERROR" in script
        assert "EXPORT_EXIT_CODE" in script


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

    def test_returns_script_when_properly_configured(self, analysis, monkeypatch):
        """Test returns export script when properly configured."""
        monkeypatch.setattr(
            "isabl_cli.gcp_lustre.get_gcp_config",
            lambda: {
                "lustre_export_enabled": True,
                "lustre_instance": "test-instance",
                "lustre_location": "us-east4-a",
                "lustre_project": "test-project",
                "lustre_mount_path": "/lustre",
                "gcs_base_uri": "gs://my-bucket",
                "lustre_poll_interval": 30,
                "lustre_max_poll_attempts": 10,
                "lustre_delete_after_export": True,
            },
        )

        result = gcp_lustre.get_export_command_for_script(analysis)
        assert result != ""
        assert "gcloud lustre instances export-data" in result


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
