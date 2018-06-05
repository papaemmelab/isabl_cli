"""cli module."""

from importlib import import_module
from os import environ
from os.path import abspath
from os.path import dirname
from os.path import expanduser
from os.path import join
from getpass import getuser
import importlib
import json
import os
import getpass
import six

from cached_property import cached_property
import pytz
import yaml

ROOT = abspath(dirname(__file__))  # make sure we use absolute paths

with open(join(ROOT, "VERSION"), "r") as version_file:
    VERSION = version_file.read().strip()

__version__ = VERSION

_DEFAULTS = {
    'API_BASE_URL': 'http://0.0.0.0:8000/api/v1',
    'GET_STORAGE_DIRECTORY': 'cli.data.get_storage_directory',
    'TRASH_ANALYSIS_STORAGE_FUNCTION': 'cli.data.trash_analysis_storage',
    'BASE_STORAGE_DIRECTORY': None,
    'FASTQ_READ_PREFIX': '',
    'ADMIN_USER': getpass.getuser(),
    'TIME_ZONE': 'America/New_York',
    'COMMANDS_LIST': [],
    'PIPELINES_SETTINGS': {},
    'INSTALLED_PIPELINES': [],
    }

_IMPORT_STRINGS = {
    'GET_STORAGE_DIRECTORY',
    'TRASH_ANALYSIS_STORAGE_FUNCTION',
    'COMMANDS_LIST',
    'INSTALLED_PIPELINES',
    }

_PATH_STRINGS = {'BASE_STORAGE_DIRECTORY'}


def perform_import(val, setting_name):
    """Perform the necessary import or imports."""
    if isinstance(val, six.string_types):
        return import_from_string(val, setting_name)
    elif isinstance(val, (list, tuple)):
        return [import_from_string(item, setting_name) for item in val]
    return val  # pragma: no cover


def import_from_string(val, setting_name):
    """Attempt to import a class from a string representation."""
    try:
        module_path, class_name = val.rsplit('.', 1)
        module = import_module(module_path)
        return getattr(module, class_name)
    except (ImportError, AttributeError) as error:  # pragma: no cover
        raise ImportError(
            f"Could not import '{val}' for setting '{setting_name}'. "
            f"{error.__class__.__name__}: {error}.")


class UserSettings(object):

    """A class used to manage user specific configurations."""

    settings_dir = join(expanduser('~'), '.bee',)
    settings_path = join(settings_dir, 'settings.json')

    def __repr__(self):  # pragma: no cover
        """Show all configurations in `settings_path`."""
        return f'UserSettings({self._read()})'

    def __setattr__(self, name, value):
        """Write attribute to `settings_path`."""
        self._write(name, value)

    def __getattr__(self, name):
        """Get attribute from `settings_path`."""
        return self._read().get(name, None)

    def _read(self):
        """Read settings.json."""
        try:
            with open(self.settings_path, 'r') as f:
                return json.load(f) or {}
        except:  # pylint: disable=W0702
            return {}

    def _write(self, name, value):
        """Write key, value to settings.json."""
        os.makedirs(self.settings_dir, exist_ok=True)

        with open(self.settings_path, 'w') as f:
            data = self._read()
            data[name] = value
            json.dump(data, f, indent=4, sort_keys=True)

        os.chmod(self.settings_path, 0o600)
        return data


class SystemSettings(object):

    """Obtained from drf-yasg."""

    def __init__(self, defaults, import_strings=None, path_strings=None):
        """See: https://github.com/axnsan12/drf-yasg/."""
        self.defaults = defaults
        self.path_strings = path_strings or []
        self.import_strings = import_strings or []

    @cached_property
    def is_admin_user(self):
        """Return True if current user is ADMIN_USER."""
        return getuser() == self.ADMIN_USER

    @cached_property
    def api_username(self):
        """Get current username from database."""
        from cli.api import api_request
        response = api_request('get', url=f'{self.API_BASE_URL}/preferences')
        return response.json()['user']

    @property
    def _settings(self):
        """Return dictionary system with settings."""
        settings = {}

        if 'CLI_SETTINGS_PATH' in environ:
            with open(environ['CLI_SETTINGS_PATH'], 'r') as f:
                settings = yaml.load(f.read())

        return settings

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
        elif attr == 'TIME_ZONE':
            val = pytz.timezone(val)

        return val


# pylint: disable=C0103
system_settings = SystemSettings(_DEFAULTS, _IMPORT_STRINGS, _PATH_STRINGS)
user_settings = UserSettings()
