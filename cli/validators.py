"""cli validators."""
import click

from glob import glob
import os

from cli import exceptions
from cli.api import get_instances


def validate_patterns_are_files(patterns, check_size=True):
    """
    Check that a list of `patterns` are valid files.

    Arguments:
        patterns (list): a list of patterns to be check.
        check_size (bool): check size is not zero for all files matched.

    Returns:
        bool: True if all patterns match existing files.
    """
    for pattern in patterns:
        files = list(glob(pattern))

        if not files:
            msg = "{} pattern matched no files.".format(pattern)
            raise exceptions.ValidationError(msg)

        for i in files:
            if not os.path.isfile(i):
                msg = "{} is not a file.".format(i)
                raise exceptions.ValidationError(msg)

            if check_size and not os.path.getsize(i) > 0:
                msg = "{} is an empty file.".format(i)
                raise exceptions.ValidationError(msg)

    return True


def validate_patterns_are_dirs(patterns):
    """
    Check that a list of `patterns` are valid dirs.

    Arguments:
        patterns (list): a list of directory patterns.

    Returns:
        bool: True if all patterns match existing directories.
    """
    for pattern in patterns:
        dirs = list(glob(pattern))

        if not dirs:
            msg = "{} pattern matched no dirs.".format(pattern)
            raise exceptions.ValidationError(msg)

        for i in dirs:
            if not os.path.isdir(i):
                msg = "{} is not a directory.".format(i)
                raise exceptions.ValidationError(msg)

    return True

def validate_pairs(pairs):
    """Get workflows for pairs."""
    workflows = {}
    ids = {w for pair in pairs for w in pair}

    for w in get_instances('workflows', ids):
        workflows[w['system_id']] = [w]

    for target, reference in pairs:
        if not target in workflows.keys():
            msg = 'Workflow {} does not exist.'.format(target)
            raise exceptions.ValidationError(msg)
        if not reference in workflows.keys():
            msg = 'Workflow {} does not exist.'.format(reference)
            raise exceptions.ValidationError(msg)
        yield workflows[str(target)], workflows[str(reference)], []

def validate_pairs_from_tuples(ctx, _, pairs):
    """Determine if workflow with ``pairs`` ID exists."""
    if not pairs: return
    for as_target, as_reference, analyses in validate_pairs(pairs):
        yield as_target, as_reference, analyses


def validate_pairs_from_file(ctx, _, path):
    """Return pairs from tsv file."""
    if not path: return
    
    pairs = []
    # Make sure is a tab delimited file of two columns.
    with open(path, "r") as f:
        line = f.readline()

        while line.startswith("#"):
            line = f.readline()

        ids = line.split("\t")
        if len(ids) != 2:
            msg = ("{} must be a tab delimited file"
                   " with two columns.").format(f.name)
            raise click.UsageError(msg)

        pairs.append((ids[0],ids[1]))

    for as_target, as_reference, analyses in validate_pairs(pairs):
        yield as_target, as_reference, analyses
