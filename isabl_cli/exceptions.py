"""Package exceptions."""

from click import ClickException
from click import UsageError


class PackageBaseException(Exception):

    """A base exception for isabl_cli."""


class ValidationError(PackageBaseException, UsageError, AssertionError):

    """A class to raise when a validation error occurs."""


class MissingRequirementError(PackageBaseException):

    """A class to raise when a requirement is missing."""


class MissingOutputError(PackageBaseException):

    """A class to raise when a file that should exist is missing."""


class ConfigurationError(PackageBaseException):

    """A class to raise when is not properly configured."""


class ApplicationNotRegisteredError(ConfigurationError, ClickException):

    """A class to raise when an application has no database object yet."""

    def __init__(self, message):
        """Show the message without a traceback when raised in a command."""
        ClickException.__init__(self, message)


class ImplementationError(PackageBaseException):

    """A class to raise when is not properly implemented."""


class CantBeRunError(PackageBaseException):

    """A class to raise when a application just cannot be run."""


class MissingDataError(PackageBaseException):

    """A class to raise when data is missing."""


class AutomationError(PackageBaseException):

    """A class to raise when automations fail."""


class TestSessionWriteError(PackageBaseException):

    """A class to raise when a test session would write to an API that isn't local."""

    __test__ = False  # not a pytest test class
