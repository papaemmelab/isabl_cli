"""Isabl CLI utils."""

from functools import update_wrapper
from getpass import getuser
from os import stat
from pathlib import Path
from pwd import getpwuid
import getpass
import json
import os
import re
import sys
import tarfile

import analytics
import click

from packaging import version
from isabl_cli.settings import system_settings


def makedirs(path, exist_ok=True, mode=0o777):
    """Make dirs ignoring umask."""
    original_umask = os.umask(0)
    os.makedirs(path, exist_ok=exist_ok, mode=mode)
    os.umask(original_umask)


def get_results(
    experiment,
    result_key,
    application_key=None,
    application_name=None,
    application_version=None,
    targets=None,
    references=None,
    analyses=None,
    status="SUCCEEDED",
    raise_error=True,
):
    """
    Match results from a experiment object.

    If targets, references or analyses are provided the analysis result must
    match these list of samples and dependencies. Match can be done by app pk, by app
    name, or by app name and cversion.

    Pass `result_key='storage_url'` to get the output directory.

    Arguments:
        experiment (dict): experiment object for which result will be retrieved.
        result_key (dict): name of the result.
        application_key (int): key of the application that generated the result.
        application_name (str): name of the application that generated the result.
        application_version (str): version of the application that generated the result.
        targets (list): target experiments dicts that must match.
        references (dict): reference experiments dicts that must match.
        analyses (dict): analyses dicts that must match.
        status (str): expected analysis status. For multiple, string with `\`.

    Returns:
        list: of tuples (result_value, analysis primary key).
    """
    targets = {i.pk for i in targets or []}
    references = {i.pk for i in references or []}
    analyses = {i.pk for i in analyses or []}

    # Filter candidates by same pk or name and/ version
    experiment_results = []
    for i in experiment.results:
        if application_key:
            if i.application.pk == application_key:
                experiment_results.append(i)
        elif application_name:
            if (
                application_version
                and application_version != "latest"
            ):
                if (
                    i.application.name == application_name
                    and i.application.version == application_version
                ):
                    experiment_results.append(i)
            else:
                if i.application.name == application_name:
                    experiment_results.append(i)

    # Filter candidates by same targets/references/analyses
    candidates = []
    for i in experiment_results:
        if targets and {j.pk for j in i.targets}.difference(
            targets
        ):  # pragma: no cover
            continue

        if references and {j.pk for j in i.references}.difference(
            references
        ):  # pragma: no cover
            continue

        if analyses and not analyses.issubset(
            {j.pk for j in i.analyses}
        ):  # pragma: no cover
            continue

        candidates.append(i)

    # Return latest if more than 1 result and version is `latest`.
    if len(candidates) > 1 and application_version == "latest":
        candidates = sorted(candidates, key=lambda x: x.pk, reverse=True)[:1]

    results = []
    for i in candidates:
        results_dict = i if result_key == "storage_url" else i.results
        result = results_dict.get(result_key)
        results.append((result, i.pk))
        if raise_error:
            assert result_key in results_dict, (
                f"Result '{result_key}' not found for analysis {i.pk}"
                f"({i.application.name} {i.application.version}) "
                f"with status: {i.status}"
            )

        assert i.status in status.split("/") if status else True, (
            f"Expected status '{status}' for result '{result_key}' did not match: "
            f"{i.pk}({i.application.name} {i.application.version}) is {i.status}"
        )
    return results


def get_result(*args, **kwargs):
    """
    See get_results for full signature.

    Arguments:
        args (list): see get_results.
        kwargs (dict): see get_results.
        application_name (str): app name to display a more explicit error.

    Returns:
        tuple: result value, analysis pk that produced the result
    """
    app_name = kwargs.get("application_name") or kwargs.get("application_key")
    results = get_results(*args, **kwargs)
    assert results, f"No results found for application: {app_name}"
    if kwargs.get("application_version") != "latest": 
        assert len(results) == 1, f"Multiple results returned {results}"
    result, key = results[0]
    return result, key


def traverse_dict(dictionary, keys, serialize=False):
    """
    Traverse a `dictionary` using a list of `keys`.

    Arguments:
        dictionary (dict): dict to be traversed.
        keys (list): keys to be explored.
        serialize (bool): force to string, if value is dict use json.dumps.

    Returns:
        str: if `serialize` is True.
        object: if `serialize` is false.
    """
    value = dictionary

    for i in keys:
        try:
            if isinstance(value, list):
                value = [j.get(i, f"INVALID KEY ({i})") for j in value]
            else:
                value = value.get(i, f"INVALID KEY ({i})")
        except AttributeError:
            value = f"INVALID KEY ({i})"

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
    decorators = reversed(decorators)

    def decorator(f):
        for i in decorators:
            f = i(f)
        return f

    return decorator


def get_rsync_command(src, dst, chmod="a-w"):
    """Get str for commant to move `src` directory to `dst`."""
    return (
        f"(chmod -R u+w {dst} || true) && "
        f"rsync -va --append-verify --remove-source-files {src}/ {dst}/ && "
        f"chmod -R {chmod} {dst} && "
        f"find {src}/ -depth -type d -empty "
        r'-exec rmdir "{}" \;'
    )

def check_rsync_version(version_stdout):
    """Check for outdated rsync lacking the --append-verify."""
    version_string =  re.search(r'.*(?P<ver>\d\.\d\.\d).*', version_stdout).group("ver")
    major_version = version.parse(version_string).major
    if major_version != 3:
        raise ValueError("Please upgrade your version of rsync!")


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
    msg = msg or f"Operation can only be performed by {admin}"

    if getpass.getuser() != admin:
        raise PermissionError(msg)


def echo_add_commit_message():
    """Echo add `--commit` flag message."""
    click.secho("\nAdd --commit to proceed.\n", fg="green", blink=True)


def echo_title(title, color="cyan", blink=False):
    """Echo a title."""
    title = "\n" + title.strip().upper() + "\n"
    title += "".join("-" for i in title.strip()) + "\n"
    click.secho(title, fg=color, file=sys.stderr, blink=blink)


def find_owner(filename):
    """Find directory owner."""
    return getpwuid(stat(filename).st_uid).pw_name


def assert_same_owner(path):
    """Validate that a path is owned by the same user."""
    try:
        assert find_owner(path) == getuser(), f"{path} must be owned by {getuser()}"
    except AssertionError as error:  # pragma: no cover
        raise click.UsageError(str(error))
    except FileNotFoundError:  # pragma: no cover
        pass


def called_from(depth=1, verbose=True):
    """Print where the current function was called from."""
    ret = sys._getframe(1).f_code.co_name + f" -> was called from: (depth={depth})\n"
    ret += "\n".join(
        "\t" + sys._getframe(i).f_code.co_name + str(i) for i in range(3, 3 + depth)
    )

    if verbose:  # pragma: no cover
        print(ret)

    return ret


def send_analytics(command):  # noqa
    """Can be used as method or decorator of click group commands to send analytics."""

    @click.pass_context
    def wrapper(ctx, *args, **kwargs):
        """Send track event after the command is executed."""
        ctx.invoke(command, *args, **kwargs)
        analytics.track(
            system_settings.api_username,
            "Ran cli command",
            {"command": command.__name__, "args": args, "kwargs": kwargs},
        )

    return update_wrapper(wrapper, command)


def first_matching_file(directory, pattern, exclude=None):
    """
    Recursively search within a dicrectory for the first file that matches the pattern.

    Args:
        directory (str): the directory to search within.
        pattern (str): a glob pattern (e.g., '*.txt') to match filenames.
        exclude (str, optional): files containing this substring will be skipped.

    Returns:
        str: the path to the first matching file as a string.

    Raises:
        NotADirectoryError: If `directory` is not a valid directory.
        FileNotFoundError: If no matching file is found.
    """
    root_path = Path(directory)
    if not root_path.is_dir():
        raise NotADirectoryError(f"{directory} is not a valid directory.")

    matching_files = (
        file for file in root_path.rglob(pattern)
        if not (exclude and exclude in file.name)
    )

    try:
        return str(next(matching_files))
    except StopIteration:
        raise FileNotFoundError(f"No file matching pattern '{pattern}' found in '{directory}'")