"""cli module."""

from os.path import abspath
from os.path import dirname
from os.path import join

from cli.app import AbstractApplication
from cli.exceptions import MissingDataError
from cli.exceptions import ValidationError

with open(join(abspath(dirname(__file__)), "VERSION"), "r") as version_file:
    __version__ = version_file.read().strip()
