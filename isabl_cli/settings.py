"""System and user settings."""

from getpass import getuser
from importlib import import_module
from os.path import abspath
from os.path import dirname
from os.path import expanduser
from os.path import join
import getpass
import json
import os

from cached_property import cached_property
from munch import Munch
import click
import pytz
import six

from isabl_cli import exceptions

_DEFAULTS = {
    "DEFAULT_LINUX_GROUP": None,
    "MAKE_STORAGE_DIRECTORY": "isabl_cli.data._make_storage_directory",
    "TRASH_ANALYSIS_STORAGE": "isabl_cli.data.trash_analysis_storage",
    "REFERENCE_DATA_IMPORTER": "isabl_cli.data.LocalReferenceDataImporter",
    "REFERENCE_GENOME_IMPORTER": "isabl_cli.data.LocalReferenceGenomeImporter",
    "DATA_IMPORTER": "isabl_cli.data.LocalDataImporter",
    "BED_IMPORTER": "isabl_cli.data.LocalBedImporter",
    "BASE_STORAGE_DIRECTORY": join(expanduser("~"), "isabl_storage"),
    "FASTQ_READ_SUFFIX": "",
    "ADMIN_USER": getpass.getuser(),
    "TIME_ZONE": "America/New_York",
    "INSTALLED_APPLICATIONS": [],
    "CUSTOM_COMMANDS": [],
    "ON_DATA_IMPORT": ["isabl_cli.data.symlink_experiment_to_projects"],
    "ON_STATUS_CHANGE": [
        "isabl_cli.data.symlink_analysis_to_targets",
        "isabl_cli.data.trigger_analyses_merge",
    ],
    "SUBMIT_MERGE_ANALYSIS": None,
    "ON_SIGNAL_FAILURE": None,
    "ADMIN_COMMANDS": [
        "isabl_cli.commands.process_finished",
        "isabl_cli.commands.patch_results",
    ],
    "SYSTEM_COMMANDS": [
        "isabl_cli.commands.get_bams",
        "isabl_cli.commands.get_bed",
        "isabl_cli.commands.get_count",
        "isabl_cli.commands.get_data",
        "isabl_cli.commands.get_metadata",
        "isabl_cli.commands.get_outdirs",
        "isabl_cli.commands.get_paths",
        "isabl_cli.commands.get_reference",
        "isabl_cli.commands.get_results",
        "isabl_cli.commands.login",
        "isabl_cli.commands.merge_project_analyses",
        "isabl_cli.commands.merge_individual_analyses",
        "isabl_cli.commands.patch_status",
        "isabl_cli.commands.rerun_signals",
        "isabl_cli.commands.run_web_signals",
    ],
}

_IMPORT_STRINGS = {
    "ADMIN_COMMANDS",
    "BED_IMPORTER",
    "CUSTOM_COMMANDS",
    "DATA_IMPORTER",
    "INSTALLED_APPLICATIONS",
    "MAKE_STORAGE_DIRECTORY",
    "ON_DATA_IMPORT",
    "ON_SIGNAL_FAILURE",
    "ON_STATUS_CHANGE",
    "REFERENCE_DATA_IMPORTER",
    "REFERENCE_GENOME_IMPORTER",
    "SUBMIT_MERGE_ANALYSIS",
    "SYSTEM_COMMANDS",
    "TRASH_ANALYSIS_STORAGE",
}

_PATH_STRINGS = {"BASE_STORAGE_DIRECTORY"}


def perform_import(val, setting_name):
    """Perform the necessary import or imports."""
    if isinstance(val, six.string_types):
        return import_from_string(val, setting_name)
    elif isinstance(val, (list, tuple)):
        return [import_from_string(item, setting_name) for item in val]
    return val  # pragma: no cover


def import_from_string(val, setting_name=None):
    """Attempt to import a class from a string representation."""
    try:
        module_path, class_name = val.rsplit(".", 1)
        module = import_module(module_path)
        return getattr(module, class_name)
    except (ImportError, AttributeError) as error:  # pragma: no cover
        raise ImportError(
            f"Could not import '{val}' for setting '{setting_name}'. "
            f"{error.__class__.__name__}: {error}."
        )


class UserSettings:

    """A class used to manage user specific configurations."""

    settings_path = join(expanduser("~"), ".isabl", "settings.json")

    def __repr__(self):  # pragma: no cover
        """Show all configurations in `settings_path`."""
        return f"UserSettings({self._read()})"

    def __setattr__(self, name, value):
        """Write attribute to `settings_path`."""
        self._write(name, value)

    def __getattr__(self, name):
        """Get attribute from `settings_path`."""
        return self._read().get(name, None)

    def _read(self):
        """Read settings.json."""
        try:
            with open(self.settings_path, "r") as f:
                return json.load(f) or {}
        except:  # pylint: disable=W0702
            return {}

    def _write(self, name, value):
        """Write key, value to settings.json."""
        os.makedirs(dirname(self.settings_path), exist_ok=True)
        data = self._read()

        with open(self.settings_path, "w") as f:
            data[name] = value
            json.dump(data, f, indent=4, sort_keys=True)

        os.chmod(self.settings_path, 0o600)
        return data


class BaseSettings:

    """Obtained from drf-yasg/rest_framework."""

    def __init__(self, defaults, import_strings=None, path_strings=None):
        """See: https://github.com/axnsan12/drf-yasg/."""
        self.defaults = defaults
        self.path_strings = path_strings or []
        self.import_strings = import_strings or []

    @property
    def _settings(self):
        """Return dictionary system with settings."""
        raise NotImplementedError

    def __getattr__(self, attr):
        """Check if present in user settings or fall back to defaults."""
        if attr not in self.defaults:  # pragma: no cover
            raise AttributeError("Invalid setting: '%s'" % attr)

        try:
            val = self._settings[attr]
        except KeyError:
            val = self.defaults[attr]

        if attr in self.import_strings:  # coerce import strings into object
            val = perform_import(val, attr)
        elif val and attr in self.path_strings:
            val = abspath(val)
        elif attr == "TIME_ZONE":
            val = pytz.timezone(val)

        return val


class SystemSettings(BaseSettings):

    client_id = os.environ.get("ISABL_CLIENT_ID")

    @cached_property
    def is_admin_user(self):
        """Return True if current user is ADMIN_USER."""
        return getuser() == self.ADMIN_USER

    @cached_property
    def api_username(self):
        """Get current username from database."""
        from isabl_cli.api import api_request

        return api_request("get", url="/rest-auth/user/").json()["username"]

    @cached_property
    def client(self):
        """Get client configuration from database."""
        from isabl_cli.api import get_instance

        if not self.client_id:
            click.secho(
                "Set environment variable ISABL_CLIENT_ID "
                "to your databased client primary key or slug "
                "to configure Isabl CLI directly from the API.",
                fg="yellow",
            )

            return {}
        return get_instance("clients", self.client_id)

    @cached_property
    def _settings(self):
        """Return dictionary system with settings."""
        return self.client.get("settings", {})


def get_application_settings(defaults, settings, reference_data, import_strings):
    """Build settings dot dict given defaults."""
    reference_data = reference_data or {}
    import_strings = import_strings or set()
    errors = []

    def _validate(default, setting, attr):
        setting = (
            NotImplemented
            if (setting is None and default == NotImplemented)
            else setting or default
        )

        if isinstance(setting, str) and "reference_data_id:" in setting:
            setting = reference_data.get(setting.split(":", 1)[1])
            setting = setting["url"] if setting else NotImplemented
        elif attr in import_strings:  # coerce import strings into object
            setting = perform_import(setting, attr)

        if setting == NotImplemented:
            errors.append(f"Missing required setting: '{attr}'")

        return setting

    def _settingfy(default, setting=None, attr=None, skip_check=False):
        if isinstance(default, dict):
            skip_check = skip_check or default.get("skip_check", False)
            setting = setting or {}
            tuples = {}

            if not isinstance(setting, dict):
                errors.append(f"Invalid setting expected dict, got: {setting}")
                setting = {}

            for i in [] if skip_check else setting.keys():
                if i not in default:
                    errors.append(f"Got unexpected setting '{i}' for '{attr}'...")

            for key, value in (
                setting.items() if skip_check and setting else default.items()
            ):
                tuples[key] = _settingfy(
                    default=value,
                    setting=setting.get(key),
                    attr=key,
                    skip_check=skip_check,  # propagate skip_check
                )

            return Munch(**tuples)

        # if the default setting is a list, it's not possible to validate inner keys
        elif isinstance(default, (list, tuple)):
            return type(default)(
                _settingfy(default=i, attr=attr, skip_check=True)
                for i in setting or default
            )

        return _validate(default, setting, attr)

    settings = _settingfy(defaults, settings, f"App Settings")

    if errors:
        raise exceptions.ConfigurationError(
            "Invalid app configuration:\n\t- " + "\n\t- ".join(errors)
        )

    return settings


# pylint: disable=C0103
system_settings = SystemSettings(_DEFAULTS, _IMPORT_STRINGS, _PATH_STRINGS)
user_settings = UserSettings()
