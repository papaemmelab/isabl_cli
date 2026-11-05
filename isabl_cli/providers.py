"""
Applications provided by installed packages through entry points.

A package lists its applications by exposing a provider function in the
`isabl_cli.applications` entry point group, for example in its `setup.py`:

    entry_points={
        "isabl_cli.applications": ["isabl_apps = isabl_apps.isabl:get_isabl_info"]
    }

The function returns `{"applications": [...]}`, where each item is either:

- an application class or import string, imported when the CLI starts just
  like the ones in INSTALLED_APPLICATIONS; or
- a dictionary from `isabl_cli.discovery.find_application_specs`, which lists the
  command from metadata and only imports the application when its command
  runs, keeping `isabl --help` fast as the number of applications grows.

Providers are only used if DISCOVER_APPLICATIONS is True. Use
EXCLUDED_APPLICATIONS to hide some of their applications.
"""

import click

from isabl_cli.discovery import cli_command_name
from isabl_cli.discovery import get_command_name
from isabl_cli.discovery import get_entry_points
from isabl_cli.discovery import import_string
from isabl_cli.discovery import is_excluded
from isabl_cli.settings import import_from_string

ENTRY_POINT_GROUP = "isabl_cli.applications"


class LazyApplication:
    """An application listed from metadata and imported when its command runs."""

    def __init__(
        self,
        import_string,  # pylint: disable=W0621
        name,
        version,
        *,
        assembly=None,
        command_name=None,
        help="",  # pylint: disable=W0622
        provider=None,
        **_,  # ignore keys added by newer versions of application_spec
    ):
        """Store the application metadata."""
        self.import_string = import_string
        self.name = name
        self.version = version
        self.assembly = assembly
        self.help = help or ""
        self.provider = provider
        self.command_name = command_name or cli_command_name(name, version)

    def __repr__(self):  # pragma: no cover
        """Show the application path."""
        return f"LazyApplication({self.import_string})"

    def load_command(self):
        """Import the application and return its click command."""
        try:
            app = import_from_string(self.import_string, ENTRY_POINT_GROUP)
            command = app.as_cli_command()
        except Exception as error:  # pylint: disable=W0703
            raise click.ClickException(
                f"Failed to load {self.import_string} "
                f"(provided by '{self.provider}'): {error}"
            )

        if command.name != self.command_name or app.ASSEMBLY != self.assembly:
            _warn(
                f"The metadata of {self.import_string} provided by '{self.provider}' "
                "is out of date, the package should regenerate it."
            )

        return command


class ApplicationGroup(click.Group):
    """A command group that imports lazy applications only when invoked."""

    def __init__(self, name=None, commands=None, lazy_applications=None, **attrs):
        """Add lazy applications, a dict of command names to LazyApplication."""
        super().__init__(name=name, commands=commands, **attrs)
        self.lazy_applications = dict(lazy_applications or {})

    def list_commands(self, ctx):
        """List loaded and lazy commands."""
        return sorted(set(self.commands) | set(self.lazy_applications))

    def get_command(self, ctx, cmd_name):
        """Import a lazy application the first time its command is requested."""
        if cmd_name not in self.commands and cmd_name in self.lazy_applications:
            command = self.lazy_applications[cmd_name].load_command()
            self.add_command(command, cmd_name)

        return self.commands.get(cmd_name)

    def format_commands(self, ctx, formatter):
        """List commands in --help, using metadata for lazy ones."""
        rows = []

        for name in self.list_commands(ctx):
            command = self.commands.get(name)

            if command is None or not command.hidden:
                rows.append((name, command))

        if not rows:
            return

        limit = formatter.width - 6 - max(len(name) for name, _ in rows)

        with formatter.section("Commands"):
            formatter.write_dl(
                [
                    (
                        name,
                        (
                            command.get_short_help_str(limit)
                            if command is not None
                            # a bare command shortens help the same way, no import needed
                            else click.Command(
                                name, help=self.lazy_applications[name].help
                            ).get_short_help_str(limit)
                        ),
                    )
                    for name, command in rows
                ]
            )


def get_providers():
    """
    Call the provider functions installed in the entry point group.

    Returns:
        list: tuples of provider name and its list of applications. Providers
        that fail are reported and skipped.
    """
    providers = []

    for entry_point in get_entry_points(ENTRY_POINT_GROUP):
        try:
            info = entry_point.load()()
            apps = info.get("applications", []) if isinstance(info, dict) else info
            providers.append((entry_point.name, list(apps)))
        except Exception as error:  # pylint: disable=W0703
            _warn(
                f"Failed to load applications from '{entry_point.name}' "
                f"({entry_point.value}): {error}",
                color="red",
            )

    return providers


def collect_applications(installed=(), excluded=(), discover=True):
    """
    Merge INSTALLED_APPLICATIONS with the applications of installed providers.

    Arguments:
        installed (list): imported INSTALLED_APPLICATIONS, they take precedence.
        excluded (list): EXCLUDED_APPLICATIONS patterns, see `is_excluded`.
        discover (bool): if False, providers are ignored.

    Returns:
        tuple: list of application classes to register now and list of
        LazyApplication to list from metadata. Applications are deduplicated
        by import string.
    """
    eager, lazy, seen = [], [], set()
    items = [(None, i) for i in installed]

    if discover:
        items += [(name, i) for name, apps in get_providers() for i in apps]

    for provider, item in items:
        try:
            if isinstance(item, dict):
                app = LazyApplication(provider=provider, **item)
                path, command_name = app.import_string, app.command_name
            else:
                if isinstance(item, str):
                    item = import_from_string(item, ENTRY_POINT_GROUP)
                app, path = item, import_string(item)
                # building command names instantiates apps, only do it if needed
                command_name = get_command_name(item) if excluded else None
        except Exception as error:  # pylint: disable=W0703
            source = provider or "INSTALLED_APPLICATIONS"
            _warn(f"Invalid application from {source}: {item} ({error})", color="red")
            continue

        if path in seen or is_excluded(path, command_name, excluded):
            continue

        seen.add(path)
        (lazy if isinstance(app, LazyApplication) else eager).append(app)

    return eager, lazy


def _warn(msg, color="yellow"):
    click.secho(msg, err=True, fg=color)
