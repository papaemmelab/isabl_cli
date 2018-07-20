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
from cli.settings import system_settings
from cli.settings import _DEFAULTS


class TestPipeline(AbstractPipeline):

    NAME = 'HELLO_WORLD'
    VERSION = 'STILL_TESTING'
    ASSEMBLY = 'GRCh4000'
    SPECIES = 'HUMAN'

    cli_help = "This is a test pipeline"
    cli_options = [options.TARGETS]
    pipeline_settings = {'foo': 'bar'}
    pipeline_inputs = {'bar': 'foo'}

    def get_project_analysis_results(self, analysis):
        return {'project_result_key': None}

    def get_analysis_results(self, analysis):
        return {'analysis_result_key': None}

    def process_cli_options(self, targets):
        return [([i], []) for i in targets]

    def merge_project_analyses(self, storage_url, analyses):
        assert len(analyses) == 2, f'Expected 2, got: {len(analyses)}'

        with open(join(storage_url, 'test.merge'), 'w') as f:
            f.write(str(len(analyses)))

    def validate_workflows(self, targets, references):
        self.validate_one_target_no_references(targets, references)

        if targets[0]['center_id'] == '0':
            raise exceptions.ValidationError('Invalid Center ID')

        return True

    def get_dependencies(self, targets, references):
        return [], {'bar': 'foo'}

    def get_command(self, analysis, inputs, settings):
        if analysis['targets'][0]['center_id'] == '1':
            return 'exit 1'

        assert inputs['bar'] == 'foo'

        return f"echo {analysis['targets'][0]['system_id']}"


def test_pipeline_settings():
    pipeline = TestPipeline()
    pipeline.pipeline_settings = {
        'test_reference': 'reference_data_id:test_id',
        'needs_to_be_implemented': NotImplemented,
        'from_system_settings': None
        }

    _DEFAULTS['PIPELINES_SETTINGS'] = {
        f'{pipeline.NAME} {pipeline.VERSION} {pipeline.ASSEMBLY}':{
            'from_system_settings': 'BAR'
        }
    }

    pipeline.assembly['reference_data']['test_id'] = dict(url='FOO')
    assert pipeline.settings.test_reference == 'FOO'

    with pytest.raises(exceptions.MissingRequirementError) as error:
        pipeline.settings.needs_to_be_implemented

    assert 'is required' in str(error.value)
    assert pipeline.settings.system_settings == system_settings
    assert pipeline.settings.from_system_settings == 'BAR'


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
    tuples = [([i], []) for i in workflows]
    command = TestPipeline.as_cli_command()
    pipeline = TestPipeline()
    ran_analyses, _, __ = pipeline.run(tuples, commit=True)

    assert 'analysis_result_key' in ran_analyses[1][0]['results']
    assert 'analysis_result_key' in ran_analyses[2][0]['results']

    runner = CliRunner()
    result = runner.invoke(command, ["--help"])

    assert "This is a test pipeline" in result.output
    assert "--commit" in result.output
    assert "--force" in result.output
    assert "--verbose" in result.output

    pks = ','.join(str(i['pk']) for i in workflows)
    args = ['-fi', 'pk__in', pks, '--verbose']
    result = runner.invoke(command, args, catch_exceptions=False)
    analysis = pipeline.get_project_analysis(project)
    merged = join(analysis['storage_url'], 'test.merge')

    assert 'FAILED' in result.output
    assert 'SUCCEEDED' in result.output
    assert 'SKIPPED 3' in result.output
    assert 'INVALID 1' in result.output
    assert isfile(merged)
    assert 'project_result_key' in analysis['results']

    with open(merged) as f:
        assert f.read().strip() == '2'

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


def test_validate_reference_genome(tmpdir):
    reference = tmpdir.join('reference.fasta')
    required = ".fai", ".amb", ".ann", ".bwt", ".pac", ".sa"
    pipeline = AbstractPipeline()

    with pytest.raises(click.UsageError) as error:
        pipeline.validate_reference_genome(reference.strpath)

    assert 'Missing indexes please run' in str(error.value)

    for i in required:
        tmp = tmpdir.join('reference.fasta' + i)
        tmp.write('foo')

    with pytest.raises(click.UsageError) as error:
        pipeline.validate_reference_genome(reference.strpath)

    assert 'samtools dict -a' in str(error.value)


def test_validate_fastq_only():
    pipeline = AbstractPipeline()
    targets = [{'sequencing_data': [], 'system_id': 'FOO'}]
    references = []

    with pytest.raises(click.UsageError) as error:
        pipeline.validate_has_raw_sequencing_data(targets, references)

    assert 'FOO' in str(error.value)

    targets = [
        {'sequencing_data': [{'file_type': 'BAM'}], 'system_id': 'FOO'},
        {'sequencing_data': [{'file_type': 'FASTQ_R1'}], 'system_id': 'BAR'},
        ]

    with pytest.raises(click.UsageError) as error:
        pipeline.validate_single_data_type(targets, references)

    assert 'FOO' in str(error.value)

    targets = [{'sequencing_data': [{'file_type': 'BAM'}], 'system_id': 'FOO'}]

    with pytest.raises(click.UsageError) as error:
        pipeline.validate_fastq_only(targets, references)

    assert 'Only FASTQ supported' in str(error.value)


def test_validate_methods():
    pipeline = AbstractPipeline()
    targets = [{'technique': {'method': 'FOO'}, 'system_id': 'FOO BAR'}]
    references = []

    with pytest.raises(click.UsageError) as error:
        pipeline.validate_methods(targets, references, 'BAR')

    assert "Only 'BAR' sequencing method allowed" in str(error.value)


def test_validate_pdx_only():
    pipeline = AbstractPipeline()
    targets = [{'specimen': {'is_pdx': False}, 'system_id': 'FOO'}]
    references = []

    with pytest.raises(click.UsageError) as error:
        pipeline.validate_pdx_only(targets, references)

    assert 'is not PDX' in str(error.value)


def test_validate_dna_rna_only():
    pipeline = AbstractPipeline()
    targets = [{'technique': {'analyte': 'DNA'}, 'system_id': 'FOO'}]
    references = []

    with pytest.raises(click.UsageError) as error:
        pipeline.validate_rna_only(targets, references)

    assert 'is not RNA' in str(error.value)

    targets = [{'technique': {'analyte': 'RNA'}, 'system_id': 'FOO'}]

    with pytest.raises(click.UsageError) as error:
        pipeline.validate_dna_only(targets, references)

    assert 'is not DNA' in str(error.value)


def test_validate_species():
    pipeline = AbstractPipeline()
    targets = [{'specimen': {'individual': {'species': 'MOUSE'}}, 'system_id': 'FOO'}]
    references = []

    with pytest.raises(click.UsageError) as error:
        pipeline.validate_species(targets, references)

    assert 'species not supported' in str(error.value)


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
