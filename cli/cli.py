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
from cli import commands
from cli import system_settings


@click.group()
@click.version_option(version=__version__)
def main():  # pragma: no cover
    """CLI command line tools."""
    pass


for i in system_settings.COMMANDS_LIST:
    main.add_command(i)

main.add_command(commands.patch_status)
