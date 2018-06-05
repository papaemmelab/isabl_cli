"""Leuktools AbstractPipelineRunner."""

# pylint: disable=W9008
from .interface import Interface


class AbstractPipeline(Interface):

    """An abstract pipeline."""

    # pipeline's name and version must be set
    NAME = None
    VERSION = None

    @staticmethod
    def get_cli_help():
        """
        Get CLI command help message.

        Returns:
            str: a string with the CLI command help.
        """
        raise NotImplementedError()

    @staticmethod
    def get_cli_options():
        """
        Get list of CLI options decorators.

        Returns:
            list: of option decoratos (click.option).
        """
        raise NotImplementedError()

    @staticmethod
    def get_tuples(**cli_options):
        """
        Must return list of tuples given the parsed options.

        Arguments:
            cli_options (dict): parsed command line options.

        Returns:
            list: of (targets, references, analyses) tuples.
        """
        raise NotImplementedError()

    def validate_tuple(self, targets, references, analyses):
        """
        Must raise UsageError if tuple combination isnt valid else return True.

        Arguments:
            targets (list): list of targets dictionaries.
            references (list): list of references dictionaries.
            analyses (list): list of analyses dictionaries.

        Raises:
            click.UsageError: if tuple is invalid.

        Returns:
            bool: True if (targets, references, analyses) combination is ok.
        """
        raise NotImplementedError()

    @staticmethod
    def get_command(analysis, settings):
        """
        Must return command and final analysis status.

        Allowed final status are SUCCEEDED and IN_PROGRESS.

        Arguments:
            analysis (dict): an analysis object as retrieved from API.
            settings (dict): pipelines settings.

        Returns:
            tuple: command (str), final status (str)
        """
        raise NotImplementedError()

    @staticmethod
    def get_status(analysis):  # pylint: disable=W0613
        """Possible values are FINISHED and IN_PROGRESS."""
        return 'FINISHED'
