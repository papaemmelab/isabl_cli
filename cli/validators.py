"""cli validators."""

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
    if not pairs:
        return []

    ids = {i for pair in pairs for i in pair}
    workflows = {i['system_id']: i for i in get_instances('workflows', ids)}
    ret = []

    for target, reference in pairs:
        if not target in workflows.keys():
            raise exceptions.ValidationError(f'Workflow {target} not found.')
        if not reference in workflows.keys():
            raise exceptions.ValidationError(f'Workflow {reference} not found.')
        ret.append((workflows[str(target)], workflows[str(reference)], []))

    return ret


def validate_pairs_from_file(ctx, _, path):
    """Return pairs from tsv file."""
    pairs = []

    with open(path, "r") as f:
        for i in f:
            if i.startswith("#"):
                continue

            try:
                ids = i.strip().split('\t')
                pairs.append((ids[0], ids[1]))
            except ValueError:
                raise exceptions.ValidationError(f'two columns required: {i}')

    return validate_pairs(pairs)
