"""Pipeline engine command line interface logic."""

from slugify import slugify
import click

from cli import utils
from .runner import Runner


COMMIT = click.option(
    "--commit",
    show_default=True,
    is_flag=True,
    help="Commit results.")


VERBOSE = click.option(
    "--verbose",
    show_default=True,
    is_flag=True,
    help="Verbose output.")


FORCE = click.option(
    "--force",
    help="Wipe all analyses and start from scratch.",
    is_flag=True)


class Interface(Runner):

    """An abstract class to run pipelines."""

    cli_help = ""
    cli_options = []

    @staticmethod
    def get_tuples(**cli_options):
        """Must return list of tuples, see pipeline.py for documentation."""
        raise NotImplementedError

    @classmethod
    def as_cli_command(cls):
        """Get pipeline as click command line interface."""
        pipe = cls()

        @click.command(help=cls.cli_help, name=pipe._get_cmd_name())
        @utils.apply_decorators(pipe._get_cli_options())
        def command(commit, force, verbose, **cli_options):
            """Click command to be used in the CLI."""
            if commit and force:
                raise click.UsageError(
                    '--commit is redundant with --force, simply use --force')

            pipe.run_tuples(
                tuples=pipe.get_tuples(**cli_options),
                commit=commit,
                force=force,
                verbose=verbose)

            if not (commit or force):
                utils.echo_add_commit_message()

        return command

    def _get_cmd_name(self):
        name = f"{self.NAME} {self.VERSION} {self.ASSEMBLY}"
        return slugify(name, separator='_')

    def _get_cli_options(self):
        """Add default options."""
        return self.cli_options + [COMMIT, FORCE, VERBOSE]
