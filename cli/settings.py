"""System and user settings."""

from getpass import getuser
from importlib import import_module
from os import environ
from os.path import abspath
from os.path import expanduser
from os.path import join
import getpass
import json
import os

from cached_property import cached_property
import pytz
import six
import yaml

from cli import exceptions

_DEFAULTS = {
    "API_BASE_URL": "http://0.0.0.0:8000/api/v1",
    "MAKE_STORAGE_DIRECTORY": "cli.data.make_storage_directory",
    "TRASH_ANALYSIS_STORAGE": "cli.data.trash_analysis_storage",
    "REFERENCE_DATA_IMPORTER": "cli.data.ReferenceDataImporter",
    "DATA_IMPORTER": "cli.data.DataImporter",
    "BED_IMPORTER": "cli.data.BedImporter",
    "BASE_STORAGE_DIRECTORY": join(expanduser("~"), "bee_storage"),
    "FASTQ_READ_PREFIX": "",
    "ADMIN_USER": getpass.getuser(),
    "TIME_ZONE": "America/New_York",
    "APPLICATIONS_SETTINGS": {},
    "INSTALLED_APPLICATIONS": [],
    "CUSTOM_COMMANDS": [],
    "ON_DATA_IMPORT": ["cli.data.symlink_experiment_to_projects"],
    "ON_STATUS_CHANGE": [
        "cli.data.symlink_analysis_to_targets",
        "cli.data.trigger_analyses_merge",
    ],
    "ON_SIGNAL_FAILURE": None,
    "ADMIN_COMMANDS": ["cli.commands.processed_finished"],
    "SYSTEM_COMMANDS": [
        "cli.commands.get_attributes",
        "cli.commands.get_bams",
        "cli.commands.get_count",
        "cli.commands.get_paths",
        "cli.commands.get_sequencing_data",
        "cli.commands.merge_project_analyses",
        "cli.commands.patch_status",
    ],
}

_IMPORT_STRINGS = {
    "ON_SIGNAL_FAILURE",
    "ON_DATA_IMPORT",
    "ON_STATUS_CHANGE",
    "MAKE_STORAGE_DIRECTORY",
    "TRASH_ANALYSIS_STORAGE",
    "ADMIN_COMMANDS",
    "SYSTEM_COMMANDS",
    "CUSTOM_COMMANDS",
    "INSTALLED_APPLICATIONS",
    "DATA_IMPORTER",
    "BED_IMPORTER",
    "REFERENCE_DATA_IMPORTER",
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


class UserSettings(object):

    """A class used to manage user specific configurations."""

    settings_dir = join(expanduser("~"), ".bee")
    settings_path = join(settings_dir, "settings.json")

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
        os.makedirs(self.settings_dir, exist_ok=True)
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
            try:
                val = environ[f"BEE_{attr}"]
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
    @cached_property
    def is_admin_user(self):
        """Return True if current user is ADMIN_USER."""
        return getuser() == self.ADMIN_USER

    @cached_property
    def api_username(self):
        """Get current username from database."""
        from cli.api import api_request

        response = api_request("get", url=f"{self.API_BASE_URL}/rest-auth/user/")
        return response.json()["username"]

    @property
    def _settings(self):
        """Return dictionary system with settings."""
        settings = {}

        if "CLI_SETTINGS_PATH" in environ:
            with open(environ["CLI_SETTINGS_PATH"], "r") as f:
                settings = yaml.load(f.read())

        return settings


class ApplicationSettings:
    def __init__(self, application, defaults, import_strings=None):
        """Get application settings from system settings."""
        self._key = f"{application.NAME} {application.VERSION} {application.ASSEMBLY}"
        self.defaults = defaults
        self.application = application.application
        self.reference_data = application.assembly["reference_data"] or {}
        self.import_strings = import_strings or {}

    @property
    def system_settings(self):
        """Return dictionary with settings."""
        return system_settings

    def __getattr__(self, attr):
        """Check if present in user settings or fall back to defaults."""
        if attr not in self.defaults:  # pragma: no cover
            raise AttributeError("Invalid setting: '%s'" % attr)

        required = self.defaults[attr] is NotImplemented

        try:
            val = self.application["settings"][attr]
        except KeyError:
            val = self.defaults[attr]

        if isinstance(val, str) and "reference_data_id:" in val:
            val = self.reference_data.get(val.split(":", 1)[1])
            val = val["url"] if val else NotImplemented
        elif attr in self.import_strings:  # coerce import strings into object
            val = perform_import(val, attr)

        # raise error if not implemented
        if val is NotImplemented or (val is None and required):
            msg = f"Setting '{attr}' is required, contact an engineer"
            raise exceptions.MissingRequirementError(msg)

        return val


# pylint: disable=C0103
system_settings = SystemSettings(_DEFAULTS, _IMPORT_STRINGS, _PATH_STRINGS)
user_settings = UserSettings()
