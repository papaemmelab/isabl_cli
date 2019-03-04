"""isabl_cli module."""

from os.path import abspath
from os.path import dirname
from os.path import join

from isabl_cli.api import create_instance
from isabl_cli.api import get_analyses
from isabl_cli.api import get_experiments
from isabl_cli.api import get_instance
from isabl_cli.api import get_instances
from isabl_cli.api import get_instances_count
from isabl_cli.api import get_projects
from isabl_cli.api import get_tree
from isabl_cli.api import get_trees
from isabl_cli.api import patch_instance
from isabl_cli.app import AbstractApplication
from isabl_cli.exceptions import MissingDataError
from isabl_cli.exceptions import ValidationError

with open(join(abspath(dirname(__file__)), "VERSION"), "r") as version_file:
    __version__ = version_file.read().strip()
