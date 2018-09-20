from os.path import join
from os.path import isfile
from click.testing import CliRunner
import click
import pytest

from cli import api
from cli import factories
from cli import exceptions
from cli import options
from cli import AbstractApplication
from cli.settings import system_settings
from cli.settings import _DEFAULTS


class TestApplication(AbstractApplication):

    NAME = 'HELLO_WORLD'
    VERSION = 'STILL_TESTING'
    ASSEMBLY = 'GRCh4000'
    SPECIES = 'HUMAN'
    HOMEPAGE = 'http://www.fake-test-app.org'

    cli_help = "This is a test application"
    cli_options = [options.TARGETS]
    application_settings = {'foo': 'bar'}
    application_inputs = {'bar': None}

    def process_cli_options(self, targets):
        return [([i], []) for i in targets]

    def validate_experiments(self, targets, references):
        self.validate_one_target_no_references(targets, references)

        if targets[0]['center_id'] == '0':
            raise AssertionError('Invalid Center ID')

        return True

    def get_dependencies(self, targets, references, settings):
        return [], {'bar': 'foo'}

    def get_command(self, analysis, inputs, settings):
        if analysis['targets'][0]['center_id'] == '1':
            return 'exit 1'

        assert inputs['bar'] == 'foo'

        return f"echo {analysis['targets'][0]['system_id']}"

    def merge_project_analyses(self, storage_url, analyses):
        assert len(analyses) == 2, f'Expected 2, got: {len(analyses)}'

        with open(join(storage_url, 'test.merge'), 'w') as f:
            f.write(str(len(analyses)))

    def get_analysis_results(self, analysis):
        return {'analysis_result_key': None}

    def get_project_analysis_results(self, analysis):
        return {'project_result_key': None}


def test_application_settings():
    application = TestApplication()
    application.application_settings = {
        'test_reference': 'reference_data_id:test_id',
        'needs_to_be_implemented': NotImplemented,
        'from_system_settings': None
        }

    application.assembly['reference_data']['test_id'] = dict(url='FOO')
    assert application.settings.test_reference == 'FOO'

    with pytest.raises(exceptions.MissingRequirementError) as error:
        application.settings.needs_to_be_implemented

    assert 'is required' in str(error.value)
    assert application.settings.system_settings == system_settings


def test_engine(tmpdir):
    data_storage_directory = tmpdir.mkdir('data_storage_directory')
    _DEFAULTS['BASE_STORAGE_DIRECTORY'] = data_storage_directory.strpath

    individual = factories.IndividualFactory(species='HUMAN')
    sample = factories.SampleFactory(individual=individual)
    project = api.create_instance('projects', **factories.ProjectFactory())

    experiments = [
        factories.ExperimentFactory(
            center_id=str(i), sample=sample, projects=[project])
        for i in range(4)
        ]

    experiments = [api.create_instance('experiments', **i) for i in experiments]
    tuples = [([i], []) for i in experiments]

    command = TestApplication.as_cli_command()
    application = TestApplication()
    ran_analyses, _, __ = application.run(tuples, commit=True)

    assert 'analysis_result_key' in ran_analyses[1][0]['results']
    assert 'analysis_result_key' in ran_analyses[2][0]['results']

    runner = CliRunner()
    result = runner.invoke(command, ["--help"])

    assert "This is a test application" in result.output
    assert "--commit" in result.output
    assert "--force" in result.output
    assert "--verbose" in result.output
    assert "--homepage" in result.output

    runner = CliRunner()
    result = runner.invoke(command, ["--homepage"])

    assert 'http://www.fake-test-app.org' in result.output

    pks = ','.join(str(i['pk']) for i in experiments)
    args = ['-fi', 'pk__in', pks, '--verbose']
    result = runner.invoke(command, args, catch_exceptions=False)
    analysis = application.get_project_analysis(project)
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


def test_validate_is_pair():
    application = AbstractApplication()
    targets = [{}]
    references = [{}]
    application.validate_is_pair(targets, references)

    with pytest.raises(AssertionError) as error:
        targets.append({})
        application.validate_is_pair(targets, references)

    assert 'Pairs only.' in str(error.value)


def test_validate_reference_genome(tmpdir):
    reference = tmpdir.join('reference.fasta')
    required = ".fai", ".amb", ".ann", ".bwt", ".pac", ".sa"
    application = AbstractApplication()

    with pytest.raises(AssertionError) as error:
        application.validate_reference_genome(reference.strpath)

    assert 'Missing indexes please run' in str(error.value)

    for i in required:
        tmp = tmpdir.join('reference.fasta' + i)
        tmp.write('foo')

    with pytest.raises(AssertionError) as error:
        application.validate_reference_genome(reference.strpath)

    assert 'samtools dict -a' in str(error.value)


def test_validate_fastq_only():
    application = AbstractApplication()
    targets = [{'sequencing_data': [], 'system_id': 'FOO'}]

    with pytest.raises(AssertionError) as error:
        application.validate_has_raw_sequencing_data(targets)

    assert 'FOO' in str(error.value)

    targets = [
        {'sequencing_data': [{'file_type': 'BAM'}], 'system_id': 'FOO'},
        {'sequencing_data': [{'file_type': 'FASTQ_R1'}], 'system_id': 'BAR'},
        ]

    with pytest.raises(AssertionError) as error:
        application.validate_single_data_type(targets)

    assert 'FOO' in str(error.value)

    targets = [{'sequencing_data': [{'file_type': 'BAM'}], 'system_id': 'FOO'}]

    with pytest.raises(AssertionError) as error:
        application.validate_fastq_only(targets)

    assert 'Only FASTQ supported' in str(error.value)


def test_validate_methods():
    application = AbstractApplication()
    targets = [{'technique': {'method': 'FOO'}, 'system_id': 'FOO BAR'}]

    with pytest.raises(AssertionError) as error:
        application.validate_methods(targets, 'BAR')

    assert "Only 'BAR' sequencing method allowed" in str(error.value)


def test_validate_pdx_only():
    application = AbstractApplication()
    targets = [{'sample': {'is_pdx': False}, 'system_id': 'FOO'}]

    with pytest.raises(AssertionError) as error:
        application.validate_pdx_only(targets)

    assert 'is not PDX' in str(error.value)


def test_validate_dna_rna_only():
    application = AbstractApplication()
    targets = [{'technique': {'analyte': 'DNA'}, 'system_id': 'FOO'}]

    with pytest.raises(AssertionError) as error:
        application.validate_rna_only(targets)

    assert 'is not RNA' in str(error.value)

    targets = [{'technique': {'analyte': 'RNA'}, 'system_id': 'FOO'}]

    with pytest.raises(AssertionError) as error:
        application.validate_dna_only(targets)

    assert 'is not DNA' in str(error.value)


def test_validate_species():
    application = AbstractApplication()
    targets = [{'sample': {'individual': {'species': 'MOUSE'}}, 'system_id': 'FOO'}]

    with pytest.raises(AssertionError) as error:
        application.validate_species(targets)

    assert 'species not supported' in str(error.value)


def test_validate_one_target_no_references():
    application = AbstractApplication()
    targets = [{}]
    references = []
    application.validate_one_target_no_references(targets, references)

    with pytest.raises(AssertionError) as error:
        references.append({})
        application.validate_one_target_no_references(targets, references)

    assert 'No reference experiments' in str(error.value)


def test_validate_atleast_onetarget_onereference():
    application = AbstractApplication()
    targets = [{}]
    references = [{}]
    application.validate_at_least_one_target_one_reference(targets, references)

    with pytest.raises(AssertionError) as error:
        targets = []
        application.validate_at_least_one_target_one_reference(
            targets, references)

    assert 'References and targets required' in str(error.value)


def test_validate_targets_not_in_references():
    application = AbstractApplication()
    targets = [{'pk': 1, 'system_id': 1}]
    references = [{'pk': 2, 'system_id': 2}]
    application.validate_targets_not_in_references(targets, references)

    with pytest.raises(AssertionError) as error:
        references = targets
        application.validate_targets_not_in_references(targets, references)

    assert '1 was also used as reference' in str(error.value)


def test_validate_dna_tuples():
    application = AbstractApplication()
    targets = [{'system_id': 1, 'technique': {'analyte': 'DNA'}}]
    references = [{'system_id': 2, 'technique': {'analyte': 'DNA'}}]
    application.validate_dna_only(targets + references)

    with pytest.raises(AssertionError) as error:
        targets[0]['technique']['analyte'] = 'RNA'
        application.validate_dna_only(targets + references)

    assert 'analyte is not DNA' in str(error.value)


def test_validate_dna_pairs():
    application = AbstractApplication()
    targets = [{'system_id': 1, 'technique': {'analyte': 'DNA'}}]
    references = [{'system_id': 2, 'technique': {'analyte': 'DNA'}}]
    application.validate_dna_pairs(targets, references)


def test_validate_same_technique():
    application = AbstractApplication()
    targets = [{'system_id': 1, 'technique': {'slug': '1'}}]
    references = [{'system_id': 2, 'technique': {'slug': '1'}}]
    application.validate_same_technique(targets, references)

    with pytest.raises(AssertionError) as error:
        targets = [{'system_id': 1, 'technique': {'slug': '2'}}]
        application.validate_same_technique(targets, references)

    assert 'Same techniques required' in str(error.value)

    with pytest.raises(AssertionError) as error:
        references.append({'system_id': 3, 'technique': {'slug': '2'}})
        application.validate_same_technique(targets, references)

    assert 'Expected one technique, got:' in str(error.value)
