"""cli utils."""

import getpass
import os
import sys
import tarfile
import json

import click

from cli import system_settings


def traverse_dict(dictionary, keys, serialize=False):
    """
    Traverse a `dictionary` using a list of `keys`.

    Arguments:
        dictionary (dict): dict to be traversed.
        keys (list): keys to be explored.
        serialize (bool): force to string, if value is dict use json.dumps.

    Returns:
        str:  if `serialize` is True.
        object: if `serialize` is false.
    """
    value = dictionary

    for i in keys:
        try:
            value = value.get(i, f'INVALID KEY ({i})')
        except AttributeError:
            value = f'INVALID KEY ({i})'

    if serialize:
        if isinstance(value, dict):
            value = json.dumps(value)
        else:
            value = str(value)

    return value


def apply_decorators(decorators):
    """
    Apply a list of decorators to callable.

    See: http://stackoverflow.com/questions/4122815
    """
    def decorator(f):
        for i in reversed(decorators):
            f = i(f)
        return f

    return decorator


def get_rsync_command(src, dst, chmod='a-w'):
    """Get str for commant to move `src` directory to `dst`."""
    return (
        f'rsync -va --append-verify --chmod={chmod} '
        f'--remove-source-files {src}/ {dst}/ && '
        f'find {src}/ -depth -type d -empty '
        r'-exec rmdir "{}" \;')


def get_tree_size(path, follow_symlinks=False):
    """Return total size of directory in bytes."""
    total = 0
    for entry in os.scandir(path):
        if entry.is_dir(follow_symlinks=follow_symlinks):
            total += get_tree_size(entry.path)
        else:
            total += entry.stat(follow_symlinks=follow_symlinks).st_size
    return total


def force_link(src, dst):
    """Force a link between src and dst."""
    try:
        os.unlink(dst)
        os.link(src, dst)
    except OSError:
        os.link(src, dst)


def force_symlink(src, dst):
    """Force a symlink between src and dst."""
    try:
        os.unlink(dst)
        os.symlink(src, dst)
    except OSError:
        os.symlink(src, dst)


def tar_dir(output_path, source_dir):
    """Compress a `source_dir` in `output_path`."""
    with tarfile.open(output_path, "w:gz") as tar:
        tar.add(source_dir, arcname=os.path.basename(source_dir))


def check_admin(msg=None):
    """Raise `PermissionError` if user is not `system_settings.ADMIN_USER`."""
    admin = system_settings.ADMIN_USER
    msg = msg or f'Operation can only be performed by {admin}'

    if getpass.getuser() != admin:
        raise PermissionError(msg)


def echo_add_commit_message():
    """Echo add `--commit` flag message."""
    click.secho('\nAdd --commit to proceed.\n', fg='green', blink=True)


def echo_title(title, color="cyan", blink=False):
    """Echo a title."""
    title = "\n" + title.strip().upper() + "\n"
    title += "".join("-" for i in title.strip()) + "\n"
    click.secho(title, fg=color, file=sys.stderr, blink=blink)
