"""Pipeline engine command line interface logic."""

import click

from cli import options
from cli import utils
from .runner import Runner


class Interface(Runner):

    """An abstract class to run pipelines."""

    @staticmethod
    def get_cli_options():
        """Must return CLI options, see pipeline.py for documentation."""
        raise NotImplementedError

    @staticmethod
    def get_cli_help():
        """Must return CLI help, see pipeline.py for documentation."""
        raise NotImplementedError

    @staticmethod
    def process_cli_options(**cli_options):
        """Must return list of tuples, see pipeline.py for documentation."""
        raise NotImplementedError

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
        return f"run_{self.pipeline['name']}_{self.pipeline['pk']}".lower()

    def _get_cli_options(self):
        """Add default options."""
        return self.get_cli_options() + [
            options.COMMIT, options.FORCE, options.QUIET]
