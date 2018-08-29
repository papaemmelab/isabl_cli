"""cli module."""

from os.path import abspath
from os.path import dirname
from os.path import join

from cli.app import AbstractApplication
from cli.exceptions import MissingDataError
from cli.exceptions import ValidationError

ROOT = abspath(dirname(__file__))  # make sure we use absolute paths

with open(join(ROOT, "VERSION"), "r") as version_file:
    VERSION = version_file.read().strip()

__version__ = VERSION
