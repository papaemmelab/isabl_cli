"""
Module that contains the command line app.

Why does this file exist, and why not put this in __main__?
You might be tempted to import things from __main__ later, but that will
cause problems, the code will get executed twice:

    - When you run `python -m cli` python will execute
      `__main__.py` as a script. That means there won't be any
      `cli.__main__` in `sys.modules`.

    - When you import __main__ it will get executed again (as a module) because
      there's no `cli.__main__` in `sys.modules`.

Also see (1) from http://click.pocoo.org/5/setuptools/#setuptools-integration
"""

import click

from cli import __version__
from cli.settings import system_settings


@click.group()
@click.version_option(version=__version__)
def main():  # pragma: no cover
    """CLI command line tools."""
    pass


@click.group()
def applications():  # pragma: no cover
    """Run registered applications."""
    pass


main.add_command(applications)

for i in system_settings.SYSTEM_COMMANDS:
    main.add_command(i)

for i in system_settings.CUSTOM_COMMANDS:  # pragma: no cover
    main.add_command(i)

for i in system_settings.INSTALLED_APPLICATIONS:  # pragma: no cover
    applications.add_command(i.as_cli_command())

if system_settings.is_admin_user:
    for i in system_settings.ADMIN_COMMANDS:
        main.add_command(i)

    if system_settings.REFERENCE_DATA_IMPORTER:
        main.add_command(
            system_settings.REFERENCE_DATA_IMPORTER.as_cli_command())

    if system_settings.DATA_IMPORTER:
        main.add_command(system_settings.DATA_IMPORTER.as_cli_command())

    if system_settings.BED_IMPORTER:
        main.add_command(system_settings.BED_IMPORTER.as_cli_command())
