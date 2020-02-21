"""
Module that contains the command line app.

Why does this file exist, and why not put this in __main__?
You might be tempted to import things from __main__ later, but that will
cause problems, the code will get executed twice:

    - When you run `python -m isabl_cli` python will execute
      `__main__.py` as a script. That means there won't be any
      `isabl_cli.__main__` in `sys.modules`.

    - When you import __main__ it will get executed again (as a module) because
      there's no `isabl_cli.__main__` in `sys.modules`.

Also see (1) from http://click.pocoo.org/5/setuptools/#setuptools-integration
"""

from collections import defaultdict
import time

from slugify import slugify
import click

from isabl_cli import __version__
from isabl_cli import exceptions
from isabl_cli.settings import system_settings
from isabl_cli.utils import send_analytics


@click.group()
@click.version_option(version=__version__)
def main():  # pragma: no cover
    """Run Isabl command line tools."""
    pass


def add_apps_groups(apps):
    """Group apps by assembly."""
    apps_dict = defaultdict(dict)

    for i in apps:  # pylint: disable=W0621
        start_time = time.time()

        try:
            command = i.as_cli_command()
            apps_dict[i.ASSEMBLY][command.name] = command
        except (exceptions.ConfigurationError, AttributeError) as error:
            click.secho(f"Invalid configuration, failed to register {i}: {error}")

        if (time.time() - start_time) > 0.1:  # pragma: no cover
            click.secho(f"{i.__name__} is loading slowly...", err=True, fg="yellow")

    for i, j in apps_dict.items():
        name = f"apps-{slugify(i)}" if i else "apps"
        help_text = f"{i} applications." if i else "data processing applications."
        main.add_command(click.Group(name=name, commands=j, help=help_text))


for i in system_settings.SYSTEM_COMMANDS:
    # Add analytics to every system command
    i.callback = send_analytics(i.callback)

    main.add_command(i)

for i in system_settings.CUSTOM_COMMANDS:  # pragma: no cover
    main.add_command(i)

try:
    if system_settings.INSTALLED_APPLICATIONS:  # pragma: no cover
        add_apps_groups(system_settings.INSTALLED_APPLICATIONS)
except ImportError as error:
    click.secho(f"Failed to import applications: {error}", fg="red")


if system_settings.is_admin_user:
    for i in system_settings.ADMIN_COMMANDS:
        main.add_command(i)

    if system_settings.REFERENCE_DATA_IMPORTER:
        main.add_command(system_settings.REFERENCE_DATA_IMPORTER.as_cli_command())

    if system_settings.REFERENCE_GENOME_IMPORTER:
        main.add_command(system_settings.REFERENCE_GENOME_IMPORTER.as_cli_command())

    if system_settings.DATA_IMPORTER:
        main.add_command(system_settings.DATA_IMPORTER.as_cli_command())

    if system_settings.BED_IMPORTER:
        main.add_command(system_settings.BED_IMPORTER.as_cli_command())
