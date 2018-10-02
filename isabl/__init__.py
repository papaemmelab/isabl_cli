"""isabl module."""

from os.path import abspath
from os.path import dirname
from os.path import join

from isabl.app import AbstractApplication
from isabl.exceptions import MissingDataError
from isabl.exceptions import ValidationError

with open(join(abspath(dirname(__file__)), "VERSION"), "r") as version_file:
    __version__ = version_file.read().strip()
