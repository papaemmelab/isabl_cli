from os.path import join
from os.path import isfile
from click.testing import CliRunner
import click
import pytest

from cli import api
from cli import factories
from cli import exceptions
from cli import options
from cli.engine import AbstractPipeline
from cli.settings import _DEFAULTS


class TestPipeline(AbstractPipeline):

    NAME = 'HELLO_WORLD'
    VERSION = 'STILL_TESTING'
    ASSEMBLY = 'GRCh4000'
    SPECIES = 'HUMAN'

    cli_help = "This is a test pipeline"
    cli_options = [options.TARGETS]

    def merge_analyses(self, storage_url, analyses):
        assert len(analyses) == 2, f'CARLOSSS {len(analyses)}'

    def get_tuples(self, targets):
        return [([i], [], []) for i in targets]

    def validate_tuple(self, targets, references, analyses):
        self.validate_one_target_no_references(targets, references)

        if targets[0]['center_id'] == '0':
            raise exceptions.ValidationError('Invalid Center ID')

        return True

    def get_command(self, analysis, settings):
        if analysis['targets'][0]['center_id'] == '1':
            return 'exit 1'
        return f"echo {analysis['targets'][0]['system_id']}"


def test_engine(tmpdir):
    data_storage_directory = tmpdir.mkdir('data_storage_directory')
    _DEFAULTS['BASE_STORAGE_DIRECTORY'] = data_storage_directory.strpath

    individual = factories.IndividualFactory(species='HUMAN')
    specimen = factories.SpecimenFactory(individual=individual)
    project = api.create_instance('projects', **factories.ProjectFactory())

    workflows = [
        factories.WorkflowFactory(
            center_id=str(i), specimen=specimen, projects=[project])
        for i in range(4)
        ]

    workflows = [api.create_instance('workflows', **i) for i in workflows]
    tuples = [([i], [], []) for i in workflows]
    command = TestPipeline.as_cli_command()
    pipeline = TestPipeline()
    pipeline.run_tuples(tuples, commit=True)

    runner = CliRunner()
    result = runner.invoke(command, ["--help"])

    assert "This is a test pipeline" in result.output
    assert "--commit" in result.output
    assert "--force" in result.output
    assert "--verbose" in result.output

    pks = ','.join(str(i['pk']) for i in workflows)
    args = ['-fi', 'pk__in', pks, '--verbose']
    result = runner.invoke(command, args, catch_exceptions=False)
    analysis = pipeline.get_project_level_analysis(project)

    assert 'FAILED' in result.output
    assert 'SUCCEEDED' in result.output
    assert 'SKIPPED 3' in result.output
    assert 'INVALID 1' in result.output

    args = ['-fi', 'pk__in', pks, '--commit', '--force']
    result = runner.invoke(command, args)
    assert '--commit is redundant with --force' in result.output

    args = ['-fi', 'pk__in', pks, '--force']
    result = runner.invoke(command, args)
    assert 'trashing:' in result.output


def test_validate_tuple_is_pair():
    pipeline = AbstractPipeline()
    targets = [{}]
    references = [{}]
    pipeline.validate_tuple_is_pair(targets, references)

    with pytest.raises(click.UsageError) as error:
        targets.append({})
        pipeline.validate_tuple_is_pair(targets, references)

    assert 'Target, reference pairs required' in str(error.value)


def test_validate_one_target_no_references():
    pipeline = AbstractPipeline()
    targets = [{}]
    references = []
    pipeline.validate_one_target_no_references(targets, references)

    with pytest.raises(click.UsageError) as error:
        references.append({})
        pipeline.validate_one_target_no_references(targets, references)

    assert 'References not allowed' in str(error.value)


def test_validate_atleast_onetarget_onereference():
    pipeline = AbstractPipeline()
    targets = [{}]
    references = [{}]
    pipeline.validate_at_least_one_target_one_reference(targets, references)

    with pytest.raises(click.UsageError) as error:
        targets = []
        pipeline.validate_at_least_one_target_one_reference(
            targets, references)

    assert 'References and targets required' in str(error.value)


def test_validate_targets_not_in_references():
    pipeline = AbstractPipeline()
    targets = [{'pk': 1, 'system_id': 1}]
    references = [{'pk': 2, 'system_id': 2}]
    pipeline.validate_targets_not_in_references(targets, references)

    with pytest.raises(click.UsageError) as error:
        references = targets
        pipeline.validate_targets_not_in_references(targets, references)

    assert '1 was also used as reference' in str(error.value)


def test_validate_dna_tuples():
    pipeline = AbstractPipeline()
    targets = [{'system_id': 1, 'technique': {'analyte': 'DNA'}}]
    references = [{'system_id': 2, 'technique': {'analyte': 'DNA'}}]
    pipeline.validate_dna_only(targets, references)

    with pytest.raises(click.UsageError) as error:
        targets[0]['technique']['analyte'] = 'RNA'
        pipeline.validate_dna_only(targets, references)

    assert 'analyte is not DNA' in str(error.value)


def test_validate_dna_pairs():
    pipeline = AbstractPipeline()
    targets = [{'system_id': 1, 'technique': {'analyte': 'DNA'}}]
    references = [{'system_id': 2, 'technique': {'analyte': 'DNA'}}]
    pipeline.validate_dna_pairs(targets, references)


def test_validate_same_technique():
    pipeline = AbstractPipeline()
    targets = [{'system_id': 1, 'technique': {'slug': '1'}}]
    references = [{'system_id': 2, 'technique': {'slug': '1'}}]
    pipeline.validate_same_technique(targets, references)

    with pytest.raises(click.UsageError) as error:
        targets = [{'system_id': 1, 'technique': {'slug': '2'}}]
        pipeline.validate_same_technique(targets, references)

    assert 'References technique differ from' in str(error.value)

    with pytest.raises(click.UsageError) as error:
        references.append({'system_id': 3, 'technique': {'slug': '2'}})
        pipeline.validate_same_technique(targets, references)

    assert 'Multiple references techniques:' in str(error.value)
