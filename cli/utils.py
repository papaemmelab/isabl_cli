"""cli utils."""

import os
import tarfile


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


def _check_admin(msg=None):
    """Raise error if user is not settings.ADMIN unless settings.TEST_MODE."""
    if settings.TEST_MODE:
        return

    if msg is None:
        msg = "Operation can only be performed by %s" % settings.ADMIN_USER

    if getpass.getuser() != settings.ADMIN_USER:
        raise PermissionError(msg)
