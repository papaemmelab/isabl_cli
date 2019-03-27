"""isabl_cli utils."""

import getpass
import json
import os
import sys
import tarfile

import click

from isabl_cli.settings import system_settings


def get_results(
    experiment,
    application_key,
    result_key,
    targets=None,
    references=None,
    analyses=None,
):
    """
    Match results from a experiment object.

    If targets, references or analyses are provided the analysis result must
    match these list of samples and dependencies.

    Pass `result_key='storage_url'` to get the output directory.

    Arguments:
        experiment (dict): experiment object for which result will be retrieved.
        application_key (int): key of the application that generated the result.
        result_key (dict): name of the result.
        targets (list): target experiments dicts that must match.
        references (dict): reference experiments dicts that must match.
        analyses (dict): analyses dicts that must match.

    Returns:
        list: of tuples (result_value, analysis primary key).
    """
    results = []
    targets = {i["pk"] for i in targets or []}
    references = {i["pk"] for i in references or []}
    analyses = {i["pk"] for i in analyses or []}

    for i in experiment["results"]:
        if i["application"]["pk"] == application_key:
            i_targets = {j["pk"] for j in i["targets"]}
            i_references = {j["pk"] for j in i["references"]}
            i_analyses = {j["pk"] for j in i["analyses"]}

            if targets and i_targets.difference(targets):
                continue

            if references and i_references.difference(references):
                continue

            if analyses and not analyses.issubset(i_analyses):
                continue

            if result_key == "storage_url":
                result = i["storage_url"]
            else:
                result = i["results"][result_key]

            results.append((result, i["pk"]))

    return results


def get_result(*args, **kwargs):
    """
    See get_experiments_results for full signature.

    Returns:
        tuple: result value, analysis pk that produced the result
    """
    results = get_results(*args, **kwargs)
    assert results, (
        f"No results were found for analysis: "
        f"{kwargs.get('application_pk') or kwargs.get('application_name')}"
    )
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
