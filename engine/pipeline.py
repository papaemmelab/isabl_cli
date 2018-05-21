"""Leuktools AbstractPipelineRunner."""

# pylint: disable=W9008
from leuktools import exceptions
from leuktools.engine.interface import Interface


class AbstractPipeline(Interface):

    """An abstract pipeline."""

    @property
    def pk(self):
        """
        Must return a unique pipeline primary key as an integer.

        Returns:
            int: the database primary key.
        """
        raise exceptions.NotImplementedError()

    @staticmethod
    def get_cli_help():
        """
        Get CLI command help message.

        Returns:
            str: a string with the CLI command help.
        """
        raise exceptions.NotImplementedError()

    @staticmethod
    def get_cli_options():
        """
        Get list of CLI options.

        Returns:
            list: of option decoratos (click.option).
        """
        raise exceptions.NotImplementedError()

    @staticmethod
    def process_cli_options(**cli_options):
        """
        Must return list of tuples given the parsed options.

        Returns:
            list: of (targets, references, analyses) tuples.
        """
        raise exceptions.NotImplementedError()

    @staticmethod
    def validate_tuple(targets, references, analyses):
        """
        Must raise UsageError if tuple combination isnt valid else return True.

        Raises:
            exceptions.UsageError: if tuple is invalid.

        Returns:
            bool: True if (targets, references, analyses) combination is ok.
        """
        raise exceptions.NotImplementedError()

    @staticmethod
    def build_command(analysis):
        """
        Must return command, final analysis status and env dictionary.

        The env dictionary should have the environment variables that will be
        available to the command.

        The allowed final status are SUCCEEDED and IN_PROGRESS.

        Returns:
            tuple: command (str), final status (str), and env (dict).
        """
        raise exceptions.NotImplementedError()

    @staticmethod
    def get_requirements(sequencing_methods):
        """
        Get submission requirements given a set of targets' methods.

        Arguments:
            sequencing_methods (set): targets sequencing methods.

        Returns:
            str: lsf requirements (e.g. -q test).
        """
        raise exceptions.NotImplementedError()
