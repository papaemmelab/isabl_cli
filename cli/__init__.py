"""cli module."""

import json
import os
from os.path import abspath
from os.path import dirname
from os.path import join
from os import environ
from os.path import expanduser

ROOT = abspath(dirname(__file__))  # make sure we use absolute paths

with open(join(ROOT, "VERSION"), "r") as version_file:
    VERSION = version_file.read().strip()

__version__ = VERSION

_DEFAULTS = {
    'API_V1_URL': 'http://0.0.0.0:8000/api/v1',
    'API_TOKEN_URL': 'http://0.0.0.0:8000/api/v1/auth/token/',
    }


class UserSettings(object):

    """A class used to manage user specific configurations."""

    settings_dir = join(expanduser('~'), '.bee',)
    settings_path = join(settings_dir, 'settings.json')

    def __repr__(self):
        return f'UserSettings({self._read_data()})'

    def __setattr__(self, name, value):
        self._write_data(name, value)

    def __getattr__(self, name):
        return self._read_data().get(name, None)

    def _read_data(self):
        try:
            with open(self.settings_path, 'r') as f:
                return json.load(f) or {}
        except:  # pylint: disable=W0702
            return {}

    def _write_data(self, name, value):
        os.makedirs(self.settings_dir, exist_ok=True)

        with open(self.settings_path, 'w') as f:
            data = self._read_data()
            data[name] = value
            json.dump(data, f, indent=4, sort_keys=True)

        os.chmod(self.settings_path, 0o600)
        return data


class SystemSettings(object):

    """Obtained from drf-yasg."""

    def __init__(self, defaults, import_strings=None):
        """See: https://github.com/axnsan12/drf-yasg/."""
        self.defaults = defaults
        self.import_strings = import_strings or []

    @property
    def _settings(self):
        """Return dictionary system with settings."""
        return environ

    def __getattr__(self, attr):
        """Check if present in user settings or fall back to defaults."""
        if attr not in self.defaults:  # pragma: no cover
            raise AttributeError("Invalid setting: '%s'" % attr)

        try:
            val = self._settings[attr]
        except KeyError:
            val = self.defaults[attr]

        if attr in self.import_strings:  # coerce import strings into object
            pass
            # val = perform_import(val, attr)

        return val


# pylint: disable=C0103
system_settings = SystemSettings(defaults=_DEFAULTS, import_strings={})
user_settings = UserSettings()
