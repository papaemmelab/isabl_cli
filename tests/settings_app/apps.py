"""Fixture apps that read their settings from the `settings.yaml` here."""

from typing import ClassVar
import uuid

from isabl_cli import AbstractApplication


class SettingsFileApp(AbstractApplication):
    """Echo the `greeting` setting from the settings manifest."""

    NAME = str(uuid.uuid4())
    VERSION = "STILL_TESTING"
    ASSEMBLY = "GRCh4000"
    SPECIES = "HUMAN"

    application_settings: ClassVar[dict[str, type]] = {
        "greeting": str,
        "reference": str,
    }

    def validate_experiments(self, targets, references):
        """Accept every tuple."""
        return True

    def get_command(self, analysis, inputs, settings):
        """Echo the greeting."""
        return f"echo {settings.greeting}"


class SettingsFileAppMissingPath(SettingsFileApp):
    """Name a reference that does not exist."""

    NAME = str(uuid.uuid4())
