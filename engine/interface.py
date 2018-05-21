"""Leuktools base validator."""

import click

from leuktools import exceptions
from leuktools import options
from leuktools import utils
from leuktools.engine.runner import Runner


class Interface(Runner):

    """An abstract class to run pipelines."""

    @staticmethod
    def get_cli_options():
        """Must return CLI options, see pipeline.py for documentation."""
        raise exceptions.NotImplementedError()

    @staticmethod
    def get_cli_help():
        """Must return CLI help, see pipeline.py for documentation."""
        raise exceptions.NotImplementedError()

    @staticmethod
    def process_cli_options(**cli_options):
        """Must return list of tuples, see pipeline.py for documentation."""
        raise exceptions.NotImplementedError()

    def _get_cli_command(self):
        """Get the click command line command."""
        @click.command(help=self.get_cli_help(), name=self._get_cmd_name())
        @utils.apply_decorators(self._get_cli_options())
        def command(commit, force, quiet, **cli_options):
            """Click command to be used in the CLI."""
            self._force = force
            self._commit = True if force else commit
            self._quiet = quiet
            self._run_tuples(self.process_cli_options(**cli_options))

        return command

    def _get_cmd_name(self):
        return f"run_{self.pipeline.name.lower()}_{self.pk}"

    def _get_cli_options(self):
        """Add default options."""
        return self.get_cli_options() + [
            options.COMMIT, options.FORCE, options.QUIET
            ]
