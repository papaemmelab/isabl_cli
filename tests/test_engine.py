from click.testing import CliRunner
import click
import pytest

from cli import _DEFAULTS
from cli import api
from cli.engine import AbstractPipeline

from . import factories

class TestPipeline(AbstractPipeline):

    NAME = 'A Test Pipeline'
    VERSION = 'A Test Version'

    @staticmethod
    def get_cli_help():
        return "This is a test pipeline"

    @staticmethod
    def get_cli_options():
        return [click.option('--filters', type=(str, str), multiple=True)]

    @staticmethod
    def get_tuples(filters, **cli_options):
        tuples = []
        filters = dict(filters)

        for i in api.get_instances('workflows', **filters):
            tuples.append(([i], [], []))

        return tuples

    def validate_tuple(self, targets, references, analyses):
        self.validate_onetarget_noreferences(targets, references, analyses)

        if targets[0]['center_id'] == '0':
            raise click.UsageError('Invalid Center ID')

        return True

    @staticmethod
    def get_command(analysis, settings):
        if analysis['targets'][0]['center_id'] == '1':
            return 'exit 1'
        return f"echo {analysis['targets'][0]['system_id']}"


def test_engine(tmpdir):
    data_storage_directory = tmpdir.mkdir('data_storage_directory')
    _DEFAULTS['BASE_STORAGE_DIRECTORY'] = data_storage_directory.strpath

    workflows = [factories.WorkflowFactory(center_id=str(i)) for i in range(3)]
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
    args = ['--filters', 'pk__in', pks, '--verbose']
    result = runner.invoke(command, args, catch_exceptions=False)

    assert 'FAILED' in result.output
    assert 'SUCCEEDED' in result.output
    assert 'SKIPPED 2' in result.output
    assert 'INVALID 1' in result.output

    args = ['--filters', 'pk__in', pks, '--commit', '--force']
    result = runner.invoke(command, args)
    assert '--commit is redundant with --force' in result.output

    args = ['--filters', 'pk__in', pks, '--force']
    result = runner.invoke(command, args)
    assert 'trashing:' in result.output


def test_validate_tuple_ispair():
    pipeline = AbstractPipeline()
    targets = [{}]
    references = [{}]
    assert pipeline.validate_tuple_ispair(targets, references, [])

    with pytest.raises(click.UsageError) as error:
        targets.append({})
        pipeline.validate_tuple_ispair(targets, references, [])

    assert 'Target, reference pairs required' in str(error.value)


def test_validate_onetarget_noreferences():
    pipeline = AbstractPipeline()
    targets = [{}]
    references = []
    assert pipeline.validate_onetarget_noreferences(targets, references, [])

    with pytest.raises(click.UsageError) as error:
        references.append({})
        pipeline.validate_onetarget_noreferences(targets, references, [])

    assert 'References not allowed' in str(error.value)


def test_validate_atleast_onetarget_onereference():
    pipeline = AbstractPipeline()
    targets = [{}]
    references = [{}]
    assert pipeline.validate_atleast_onetarget_onereference(targets, references, [])

    with pytest.raises(click.UsageError) as error:
        targets = []
        pipeline.validate_atleast_onetarget_onereference(targets, references, [])

    assert 'References and targets required' in str(error.value)


def test_validate_targets_notin_references():
    pipeline = AbstractPipeline()
    targets = [{'pk': 1, 'system_id': 1}]
    references = [{'pk': 2, 'system_id': 2}]
    assert pipeline.validate_targets_notin_references(targets, references, [])

    with pytest.raises(click.UsageError) as error:
        references = targets
        pipeline.validate_targets_notin_references(targets, references, [])

    assert '1 was also used as reference' in str(error.value)


def test_validate_dna_tuples():
    pipeline = AbstractPipeline()
    targets = [{'system_id': 1, 'technique': {'analyte': 'DNA'}}]
    references = [{'system_id': 2, 'technique': {'analyte': 'DNA'}}]
    assert pipeline.validate_dna_tuples(targets, references, [])

    with pytest.raises(click.UsageError) as error:
        targets[0]['technique']['analyte'] = 'RNA'
        pipeline.validate_dna_tuples(targets, references, [])

    assert 'analyte is not DNA' in str(error.value)


def test_validate_dna_pairs():
    pipeline = AbstractPipeline()
    targets = [{'system_id': 1, 'technique': {'analyte': 'DNA'}}]
    references = [{'system_id': 2, 'technique': {'analyte': 'DNA'}}]
    assert pipeline.validate_dna_pairs(targets, references, [])


def test_validate_same_technique():
    pipeline = AbstractPipeline()
    targets = [{'system_id': 1, 'technique': {'slug': '1'}}]
    references = [{'system_id': 2, 'technique': {'slug': '1'}}]
    assert pipeline.validate_same_technique(targets, references, [])

    with pytest.raises(click.UsageError) as error:
        targets = [{'system_id': 1, 'technique': {'slug': '2'}}]
        pipeline.validate_same_technique(targets, references, [])

    assert 'References technique differ from' in str(error.value)

    with pytest.raises(click.UsageError) as error:
        references.append({'system_id': 3, 'technique': {'slug': '2'}})
        pipeline.validate_same_technique(targets, references, [])

    assert 'Multiple references techniques:' in str(error.value)
