"""System and user settings."""

from __future__ import annotations

from getpass import getuser
from importlib import import_module
from os.path import abspath
from os.path import dirname
from os.path import expanduser
from os.path import join
from os.path import realpath
from os.path import relpath
from typing import TypeAlias
import getpass
import inspect
import json
import os
import re
import subprocess

from cached_property import cached_property
from munch import Munch
from munch import munchify
import analytics
import click
import pytz
import six
import yaml

from isabl_cli import exceptions

_DEFAULTS = {
    "SUBMIT_ANALYSES": "isabl_cli.batch_systems.submit_local",
    "SUBMIT_CONFIGURATION": {},
    "DEFAULT_LINUX_GROUP": None,
    "MAKE_STORAGE_DIRECTORY": "isabl_cli.data._make_storage_directory",
    "TRASH_ANALYSIS_STORAGE": "isabl_cli.data.trash_analysis_storage",
    "REFERENCE_DATA_IMPORTER": "isabl_cli.data.LocalReferenceDataImporter",
    "REFERENCE_GENOME_IMPORTER": "isabl_cli.data.LocalReferenceGenomeImporter",
    "DATA_IMPORTER": "isabl_cli.data.LocalDataImporter",
    "YAML_DATA_IMPORTER": "isabl_cli.data.LocalYamlDataImporter",
    "BED_IMPORTER": "isabl_cli.data.LocalBedImporter",
    "BASE_STORAGE_DIRECTORY": join(expanduser("~"), "isabl_storage"),
    "FASTQ_READ_SUFFIX": "",
    "ADMIN_USER": getpass.getuser(),
    "CLIENT_NAME": None,
    "TIME_ZONE": "America/New_York",
    "INSTALLED_APPLICATIONS": [],
    "CUSTOM_COMMANDS": [],
    "ON_DATA_IMPORT": [],
    "ON_STATUS_CHANGE": [],
    "CREATE_SYMLINKS": True,
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
        "isabl_cli.commands.merge_individual_analyses",
        "isabl_cli.commands.merge_project_analyses",
        "isabl_cli.commands.patch_status",
        "isabl_cli.commands.reject_analysis",
        "isabl_cli.commands.rerun_signals",
        "isabl_cli.commands.run_signals",
        "isabl_cli.commands.run_web_signals",
    ],
    "EXTRA_RAW_DATA_FORMATS": [],
}

_IMPORT_STRINGS = {
    "SUBMIT_ANALYSES",
    "ADMIN_COMMANDS",
    "BED_IMPORTER",
    "CUSTOM_COMMANDS",
    "DATA_IMPORTER",
    "YAML_DATA_IMPORTER",
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


def import_valid_applications(val, setting_name):
    """Import all valid applications strings specified on the client settings."""
    if setting_name == "INSTALLED_APPLICATIONS":
        apps = []
        for item in val:
            try:
                apps.append(import_from_string(item, setting_name))
            except ImportError as error:
                click.secho(f"Failed to import applications: {error}", err=True, fg="red")
        return apps
    return [import_from_string(item, setting_name) for item in val]


def perform_import(val, setting_name):
    """Perform the necessary import or imports."""
    if isinstance(val, six.string_types):
        return import_from_string(val, setting_name)
    elif isinstance(val, (list, tuple)):
        return import_valid_applications(val, setting_name)
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
    api_url = os.environ.get("ISABL_API_URL", "http://localhost:8000/api/v1/")

    def __repr__(self):  # pragma: no cover
        """Show all configurations in `settings_path`."""
        return f"UserSettings({self._read()})"

    def __setattr__(self, name, value):
        """Write attribute to `settings_path`."""
        self._write(name, value)

    def __getattr__(self, name):
        """Get attribute from `settings_path`."""
        return self._read().get(self.api_url, {}).get(name, None)

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
            if self.api_url not in data:
                data[self.api_url] = {}
            data[self.api_url][name] = value
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
                err=True,
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


SETTINGS_SCHEMA = "isabl-app-settings/1"

_Declared: TypeAlias = "type | tuple[type, ...] | dict[str, _Declared]"

_HEAD_JOB: _Declared = {
    "partition": str,
    "cores": int,
    "memory_gb": int,
    "minutes": int,
}


class _StrictLoader(yaml.SafeLoader):
    """Collect anchors, aliases and duplicate keys as faults."""

    def __init__(self, stream: str) -> None:
        """Start with no faults."""
        super().__init__(stream)
        self.faults: list[str] = []

    def compose_node(self, parent: yaml.Node | None, index: object) -> yaml.Node | None:
        """Record a fault for each anchor or alias, then compose the node."""
        event = self.peek_event()
        line = event.start_mark.line + 1
        if isinstance(event, yaml.AliasEvent):
            self.faults.append(f"YAML alias at line {line}: copy the value instead")
        elif getattr(event, "anchor", None):
            self.faults.append(f"YAML anchor at line {line}: copy the value instead")
        return super().compose_node(parent, index)

    def construct_mapping(
        self, node: yaml.MappingNode, deep: bool = False
    ) -> dict[object, object]:
        """Record a fault for each duplicate key, then construct the map."""
        # PyYAML keeps the last duplicate without a warning.
        seen = set()
        for key_node, _ in node.value:
            if isinstance(key_node, yaml.ScalarNode):
                if key_node.value in seen:
                    line = key_node.start_mark.line + 1
                    self.faults.append(
                        f"duplicate key '{key_node.value}' at line {line}"
                    )
                seen.add(key_node.value)
        return super().construct_mapping(node, deep)


def _type_name(declared: type | tuple[type, ...]) -> str:
    if isinstance(declared, tuple):
        return " or ".join(i.__name__ for i in declared)
    return declared.__name__


def _check(value: object, declared: _Declared, key: str, faults: list[str]) -> None:
    # Append the faults of one value against its declared type or nested map.
    if isinstance(declared, dict):
        if not isinstance(value, dict):
            _check(value, dict, key, faults)
            return
        prefix = f"{key}." if key else ""
        for name in declared:
            if name not in value:
                faults.append(f"missing setting '{prefix}{name}'")
        for name, item in value.items():
            if name in declared:
                _check(item, declared[name], f"{prefix}{name}", faults)
            else:
                faults.append(f"unknown setting '{prefix}{name}'")
        return

    allowed = declared if isinstance(declared, tuple) else (declared,)
    if not all(isinstance(i, type) for i in allowed):
        faults.append(
            f"declaration of '{key}' must be a type, got {declared!r}: "
            "fix application_settings"
        )
        return
    # A YAML boolean is a Python int, so isinstance alone accepts it as int.
    if not isinstance(value, allowed) or (
        isinstance(value, bool) and bool not in allowed
    ):
        faults.append(
            f"setting '{key}' must be {_type_name(declared)}, "
            f"got {type(value).__name__} {value!r}"
        )


def settings_file(app_class: type) -> str | None:
    """Return the `settings.yaml` next to the module of `app_class`, if any."""
    path = join(dirname(inspect.getfile(app_class)), "settings.yaml")
    return path if os.path.isfile(path) else None


def load_app_settings(app_class: type, api_url: str | None = None) -> Munch:
    """Load the section of `app_class` from its settings manifest.

    :param api_url: The Isabl instance whose block to read. Defaults to
        `ISABL_API_URL`. A trailing slash does not matter.
    :raises ConfigurationError: If the manifest has one or more faults.
        The message lists every fault.
    """
    api_url = api_url or os.environ.get(
        "ISABL_API_URL", "http://localhost:8000/api/v1/"
    )
    name = app_class.__name__
    path = settings_file(app_class)
    faults: list[str] = []

    def _error() -> exceptions.ConfigurationError:
        return exceptions.ConfigurationError(
            f"Invalid settings manifest {path} [{api_url} {name}]:\n"
            + "\n".join(f"  - {i}" for i in faults)
        )

    if path is None:
        faults.append(
            f"no settings.yaml next to {inspect.getfile(app_class)}: create it"
        )
        raise _error()

    with open(path, encoding="utf-8") as stream:
        loader = _StrictLoader(stream.read())
    try:
        data = loader.get_single_data()
    except yaml.YAMLError as error:
        faults.append(f"YAML syntax error: {error}")
        raise _error() from error
    finally:
        loader.dispose()
    faults += loader.faults

    if not isinstance(data, dict):
        faults.append(f"the file must be a map, got {type(data).__name__}")
        raise _error()
    if data.get("schema") != SETTINGS_SCHEMA:
        faults.append(f"schema must be '{SETTINGS_SCHEMA}', got {data.get('schema')!r}")

    blocks = [i for i in data if i != "schema"]
    block = next(
        (data[i] for i in blocks if str(i).rstrip("/") == api_url.rstrip("/")), None
    )
    if not isinstance(block, dict):
        faults.append(
            f"no block for {api_url} (blocks: {', '.join(map(str, blocks)) or 'none'})"
        )
        raise _error()

    module = inspect.getmodule(app_class)
    for key in block:
        defined = getattr(module, str(key), None)
        if not (
            inspect.isclass(defined) and defined.__module__ == app_class.__module__
        ):
            faults.append(
                f"unknown class '{key}': {app_class.__module__} defines no such class"
            )

    section = block.get(name)
    if not isinstance(section, dict):
        faults.append(f"no section for {name}")
    else:
        declared = getattr(app_class, "application_settings", {})
        values = {k: v for k, v in section.items() if k != "head_job"}
        _check(values, declared, "", faults)
        head_job = section.get("head_job", {})
        _check(head_job, dict, "head_job", faults)
        if isinstance(head_job, dict):
            for method, resources in head_job.items():
                _check(resources, _HEAD_JOB, f"head_job.{method}", faults)

    if faults:
        raise _error()
    return munchify(section)


def setting_paths(values: object, key: str = "") -> list[tuple[str, str]]:
    """Return `(dotted_key, path)` for every path in the setting values.

    A path is a piece of a string value that starts with `/` but not
    with `//`. The pieces come from a split on whitespace, commas,
    colons, quotes and `=`. The split finds binds and quoted paths in
    Nextflow text. It does not find an unquoted path in Groovy code.

    :param key: The dotted key of `values`. Recursive calls set it.
    """
    if isinstance(values, dict):
        items = values.items()
    elif isinstance(values, list):
        items = enumerate(values)
    elif isinstance(values, str):
        pieces = re.split(r"[\s,:'\"=]+", values)
        return [
            (key, i) for i in pieces if i.startswith("/") and not i.startswith("//")
        ]
    else:
        return []
    prefix = f"{key}." if key else ""
    return [j for name, item in items for j in setting_paths(item, f"{prefix}{name}")]


def check_paths(values: object, production: list[str] | None = None) -> list[str]:
    """Return a fault for each path in the values that is not valid.

    :param production: Folder prefixes, such as `/data1/papaemme/isabl/`.
        If given, a path outside all of the prefixes is a fault.
    """
    faults = []
    for key, path in setting_paths(values):
        if not os.path.exists(path):
            faults.append(f"path does not exist: {key} = {path}")
        # The added slash lets a production folder itself pass.
        if production is not None and not f"{path}/".startswith(tuple(production)):
            faults.append(f"path outside production locations: {key} = {path}")
    return faults


def _jsonable(value: object) -> object:
    # Keep JSON types. Write callables as import strings, the rest with str.
    if isinstance(value, dict):
        return {str(k): _jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, set)):
        return [_jsonable(i) for i in value]
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    if callable(value) and hasattr(value, "__qualname__"):
        return f"{value.__module__}.{value.__qualname__}"
    return str(value)


def settings_record(app_class: type, values: object, api_url: str) -> dict:
    """Return the record of the settings that a run of `app_class` uses.

    The record holds the manifest path relative to the checkout, the
    instance, the section, the schema, the values, the commit, and
    whether the checkout has local edits. `file` and `schema` are None
    for an app without a settings manifest.

    :raises ConfigurationError: If the module of `app_class` is not in
        a git checkout.
    """
    folder = dirname(inspect.getfile(app_class))

    def git(*args: str) -> str:
        # The bot account does not own every checkout, hence safe.directory.
        command = ["git", "-c", "safe.directory=*", "-C", folder, *args]
        try:
            return subprocess.run(
                command, check=True, capture_output=True, text=True
            ).stdout
        except (OSError, subprocess.CalledProcessError) as error:
            detail = getattr(error, "stderr", None) or error
            raise exceptions.ConfigurationError(
                f"Can't read the git commit of {app_class.__name__} in {folder}: "
                f"{str(detail).strip()}. Run apps from a git checkout."
            ) from error

    top, commit = git("rev-parse", "--show-toplevel", "HEAD").splitlines()
    path = settings_file(app_class)
    return {
        "file": relpath(realpath(path), top) if path else None,
        "instance": api_url,
        "section": app_class.__name__,
        "schema": SETTINGS_SCHEMA if path else None,
        "values": _jsonable(values),
        "commit": commit,
        "dirty": bool(git("status", "--porcelain", "--untracked-files=no")),
    }


# load analytics token
analytics.write_key = "y1bj5g2hqjZJXEA2berZ7WIick8IDSdc"

# pylint: disable=C0103
system_settings = SystemSettings(_DEFAULTS, _IMPORT_STRINGS, _PATH_STRINGS)
user_settings = UserSettings()
