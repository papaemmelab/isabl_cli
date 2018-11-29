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
from slugify import slugify
import click

from isabl_cli import __version__
from isabl_cli.settings import system_settings


@click.group()
@click.version_option(version=__version__)
def main():  # pragma: no cover
    """run isabl command line tools."""
    pass


def add_apps_groups(apps):
    """Group apps by assembly."""
    apps_dict = defaultdict(dict)

    for i in apps:  # pylint: disable=W0621
        command = i.as_cli_command()
        apps_dict[i.ASSEMBLY][command.name] = command

    for i, j in apps_dict.items():
        name, _help = f"apps-{slugify(i)}", f"{i} applications."
        main.add_command(click.Group(name=name, commands=j, help=_help))


for i in system_settings.SYSTEM_COMMANDS:
    main.add_command(i)

for i in system_settings.CUSTOM_COMMANDS:  # pragma: no cover
    main.add_command(i)

if system_settings.INSTALLED_APPLICATIONS:
    add_apps_groups(system_settings.INSTALLED_APPLICATIONS)

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
