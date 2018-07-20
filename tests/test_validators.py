"""cli validators tests."""

from os.path import join

import click
import pytest

from cli import api
from cli import exceptions
from cli import factories
from cli import validators


def test_validate_patterns_are_files(tmpdir):
    """Create multiple files and test test_validate_patterns_are_files."""
    tmpdir_path = str(tmpdir)

    for i in range(5):
        with open(join(tmpdir_path, "empty" + str(i)), "w"):
            pass

    for i in range(11, 15):
        with open(join(tmpdir_path, "not_empty" + str(i)), "w") as f:
            f.write("I'm not empty.")

    empty = [join(tmpdir_path, "empty*")]
    not_empty = [join(tmpdir_path, "not_empty*")]
    none_existing = ["florentino", "ariza*"]
    none_file = [tmpdir_path]

    # check empty files exist
    assert validators.validate_patterns_are_files(empty, check_size=False)

    # check files exist amd are not empty
    assert validators.validate_patterns_are_files(not_empty, check_size=True)

    # check that empty files raise error with default setting
    with pytest.raises(exceptions.ValidationError):
        validators.validate_patterns_are_files(empty)

    # check that empty files raise error with flag
    with pytest.raises(exceptions.ValidationError):
        validators.validate_patterns_are_files(empty, check_size=True)

    # check that pattern is not file
    with pytest.raises(exceptions.ValidationError):
        validators.validate_patterns_are_files(none_file, check_size=True)

    # check that empty files raise error with flag
    with pytest.raises(exceptions.ValidationError):
        validators.validate_patterns_are_files(empty, check_size=True)

    # check that invalid patterns raise validationerror error
    with pytest.raises(exceptions.ValidationError):
        validators.validate_patterns_are_files(none_existing, check_size=True)


def test_validate_patterns_are_dirs(tmpdir):
    """test_validate_patterns_are_dirs."""
    tmpdir_path = str(tmpdir)
    file_patterns = [join(tmpdir_path, "a_file")]
    existing_patterns = [tmpdir_path]
    none_existing_patterns = ["florentino", "ariza*"]

    # touch the file
    with open(file_patterns[0], "w"):
        pass

    # check empty files exist
    assert validators.validate_patterns_are_dirs(existing_patterns)

    # check that empty files raise error with default setting
    with pytest.raises(exceptions.ValidationError):
        validators.validate_patterns_are_dirs(none_existing_patterns)

    # check that empty files raise error with flag
    with pytest.raises(exceptions.ValidationError):
        validators.validate_patterns_are_dirs(file_patterns)

def test_validate_pairs(tmpdir):
    """test_validate_pairs."""
    data = [factories.WorkflowFactory() for i in range(3)]
    instance_a = api.create_instance('workflows', **data[0])
    instance_b = api.create_instance('workflows', **data[1])
    instance_c = api.create_instance('workflows', **data[2])

    # pairs that exist
    real_pairs = [
        (instance_a['system_id'],instance_b['system_id']),
        (instance_c['system_id'],instance_b['system_id'])
    ]
    parsed_pairs = validators.validate_pairs(real_pairs)
    assert parsed_pairs[0] == ([instance_a], [instance_b])
    assert parsed_pairs[1] == ([instance_c], [instance_b])

    # pairs that dont exist
    fake_pair = [("not_real", 0), (2, "fake")]
    with pytest.raises(exceptions.ValidationError):
        validators.validate_pairs(fake_pair)

def test_validate_pairs_from_file(tmpdir):
    """test_validate_pairs_from_file."""
    tmpdir_path = str(tmpdir)
    data = [factories.WorkflowFactory() for i in range(3)]
    instance_a = api.create_instance('workflows', **data[0])
    instance_b = api.create_instance('workflows', **data[1])
    instance_c = api.create_instance('workflows', **data[2])

    realpairdf = join(tmpdir_path, 'realpairdf.tsv')
    with open(realpairdf, 'w') as f:
        f.write('# I am the header\n')
        f.write(instance_a['system_id'] + '\t' + instance_b['system_id'] + '\n')
        f.write(instance_c['system_id'] + '\t' + instance_b['system_id'] + '\n')

    badpairdf = join(tmpdir_path, 'badpairdf.tsv')
    with open(badpairdf, 'w') as f:
        f.write('\t'.join([instance_a['system_id']]))

    # valid pair file
    parsed_pairs = validators.validate_pairs_from_file(None, None, realpairdf)
    assert parsed_pairs[0] == ([instance_a], [instance_b])
    assert parsed_pairs[1] == ([instance_c], [instance_b])

    # invalid pair file
    with pytest.raises(click.UsageError):
        validators.validate_pairs_from_file(None, None, badpairdf)
