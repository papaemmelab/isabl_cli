"""isabl_cli utils tests."""

from getpass import getuser
from os.path import join
import os
import tarfile

import pytest

from isabl_cli import utils
from isabl_cli.settings import _DEFAULTS


def test_find_user(tmpdir):
    assert getuser() == utils.find_owner(str(tmpdir))
    assert utils.assert_same_owner(str(tmpdir)) is None


def test_force_symlink(tmpdir):
    src = join(str(tmpdir), "src")
    dst = join(str(tmpdir), "dst")

    with open(src, "w") as f:
        f.write("Not empty.")

    utils.force_symlink(src, dst)
    assert os.path.islink(dst)


def test_force_symlink_overwrite(tmpdir):
    src = join(str(tmpdir), "src")
    dst = join(str(tmpdir), "dst")

    with open(src, "w") as f:
        f.write("Correct.")

    with open(dst, "w") as f:
        f.write("Wrong.")

    utils.force_symlink(src, dst)
    assert os.path.islink(dst)

    with open(dst, "r") as f:
        assert "Correct" in f.read()


def test_force_link(tmpdir):
    src = join(str(tmpdir), "src")
    dst = join(str(tmpdir), "dst")

    with open(src, "w") as f:
        f.write("Not empty.")

    utils.force_link(src, dst)
    assert os.path.isfile(dst)
    assert not os.path.islink(dst)


def test_force_link_overwrite(tmpdir):
    src = join(str(tmpdir), "src")
    dst = join(str(tmpdir), "dst")

    with open(src, "w") as f:
        f.write("Correct.")

    with open(dst, "w") as f:
        f.write("Wrong.")

    utils.force_link(src, dst)
    assert os.path.isfile(dst)
    assert not os.path.islink(dst)

    with open(dst, "r") as f:
        assert "Correct" in f.read()


def test_tar_dir(tmpdir):
    dst_dir = join(str(tmpdir), "dst_dir")
    source_dir = join(str(tmpdir), "source_dir")
    output_path = join(str(tmpdir), "source_dir.tar.gz")

    files = (
        (join(source_dir, "1"), "first file."),
        (join(source_dir, "2"), "second file."),
        (join(source_dir, "3"), "third file."),
    )

    os.makedirs(source_dir)
    os.makedirs(dst_dir)

    for i, j in files:
        with open(i, "w") as f:
            f.write(j)

    utils.tar_dir(output_path=output_path, source_dir=source_dir)
    assert tarfile.is_tarfile(output_path)

    # Remove files and dir
    for i in files:
        os.unlink(i[0])

    os.rmdir(source_dir)
    tar = tarfile.open(output_path)
    tar.extractall(path=dst_dir)
    tar.close()

    for i, j in files:
        i = i.replace(source_dir, join(dst_dir, "source_dir"))
        with open(i, "r") as f:
            assert j in f.read()


def test_get_tree_size(tmpdir):
    tmpdir.mkdir("l1").mkdir("l2").join("test").write("foo")
    tmpdir.mkdir("l3").mkdir("l2").join("test").write("foo")
    sym_false = utils.get_tree_size(tmpdir.join("l1"), follow_symlinks=False)
    os.symlink(tmpdir.join("l3"), tmpdir.join("l1").join("l3"))
    sym_true = utils.get_tree_size(tmpdir.join("l1"), follow_symlinks=True)

    assert sym_false
    assert sym_true
    assert sym_true > sym_false


def test_check_admin():
    admin = _DEFAULTS["ADMIN_USER"]
    _DEFAULTS["ADMIN_USER"] = "not the admin"

    with pytest.raises(PermissionError) as error:
        utils.check_admin()

    assert "not the admin" in str(error.value)
    _DEFAULTS["ADMIN_USER"] = admin


def test_traverse_dict():
    assert utils.traverse_dict({"1": {}}, ["1"], serialize=True) == "{}"
    assert utils.traverse_dict({"1": [{"2": 1}]}, ["1", "2"], serialize=True) == "[1]"


def test_called_from():
    def a_function():
        return utils.called_from(3, verbose=False)

    assert "a_function" in a_function()
    assert "test_called_from" in a_function()
    assert "pytest_pyfunc_call" in a_function()

def test_rsync_version_checking():
    incompatible_version_string = "rsync  version 2.6.9  protocol version 29"
    with pytest.raises(ValueError):
        utils.check_rsync_version(incompatible_version_string)
    compatible_version_string = "rsync  version 3.2.7  protocol version 31"
    assert utils.check_rsync_version(compatible_version_string) is None
