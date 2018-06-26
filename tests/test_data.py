import os
import re
import uuid

from click.testing import CliRunner
import click
import pytest
import yaml

from cli import api
from cli import data
from cli import factories
from cli.settings import _DEFAULTS


def test_trash_analysis_storage():
    # the rest is tested in test_engine
    with pytest.raises(click.UsageError) as error:
        data.trash_analysis_storage({'status': 'SUCCEEDED'})
    assert "You can't wipe a succeeded analysis" in str(error.value)


def test_make_storage_directory(tmpdir):
    i = data.make_storage_directory(tmpdir.strpath, 'test', 12345, use_hash=True)
    j = data.make_storage_directory(tmpdir.strpath, 'test', 12345, use_hash=False)
    assert '/test/23/45/12345' in i
    assert '/test/12345' in j


def test_import_reference_data(tmpdir):
    data_storage_directory = tmpdir.mkdir('data_storage_directory')
    _DEFAULTS['BASE_STORAGE_DIRECTORY'] = str(data_storage_directory)
    assembly_name = uuid.uuid4().hex
    species = 'HUMAN'
    reference_test = tmpdir.join('test.fasta')
    reference_test.write('foo')

    assembly = data.ReferenceDataImporter.import_data(
        assembly=assembly_name,
        species=species,
        data_src=reference_test.strpath,
        description='test description',
        data_id='reference_link',
        symlink=True)

    assert os.path.islink(assembly['reference_data']['reference_link']['url'])
    assert assembly['reference_data']['reference_link']['description'] == 'test description'

    assembly = data.ReferenceDataImporter.import_data(
        assembly=assembly_name,
        species=species,
        data_src=reference_test.strpath,
        description='test description',
        data_id='reference_move',
        symlink=False)

    assert os.path.isfile(assembly['reference_data']['reference_move']['url'])
    assert not os.path.islink(assembly['reference_data']['reference_move']['url'])

    reference_test.write('foo')
    command = data.ReferenceDataImporter.as_cli_command()
    runner = CliRunner()
    args = [
        '--assembly', assembly_name,
        '--assembly', species,
        '--data-src', reference_test.strpath,
        '--data-id', 'reference_move',
        '--description', 'Test',
        ]

    result = runner.invoke(command, args, catch_exceptions=False)
    assert 'has already reference data registered' in result.output


def test_import_bedfiles(tmpdir):
    data_storage_directory = tmpdir.mkdir('data_storage_directory')
    _DEFAULTS['BASE_STORAGE_DIRECTORY'] = str(data_storage_directory)
    technique = api.create_instance('techniques', **factories.TechniqueFactory())
    targets = tmpdir.join('targets.bed')
    baits = tmpdir.join('baits.bed')
    targets.write('2\t1\t2\n1\t1\t2\n')
    baits.write('2\t1\t2\n1\t1\t2\n')
    species = 'HUMAN'

    technique = data.BedImporter.import_bedfiles(
        species=species,
        technique_key=technique['pk'],
        targets_path=targets.strpath,
        baits_path=baits.strpath,
        assembly='AnAssembly',
        description='This are test bedfiles')

    for i in 'targets', 'baits':
        for j in '', '.gz', '.gz.tbi':
            assert os.path.isfile(technique['bedfiles']['AnAssembly'][i] + j)

        with open(technique['bedfiles']['AnAssembly'][i], 'r') as f:  # test bed is sorted
            assert next(f).startswith('1')

    command = data.BedImporter.as_cli_command()
    runner = CliRunner()
    args = [
        '--targets-path', targets.strpath,
        '--baits-path', baits.strpath,
        '--key', technique['pk'],
        '--assembly', 'AnAssembly',
        '--species', species,
        '--description', 'Test',
        ]
    result = runner.invoke(command, args, catch_exceptions=False)
    assert 'has registered bedfiles for' in result.output


def test_local_data_import(tmpdir):
    data_storage_directory = tmpdir.mkdir('data_storage_directory')
    _DEFAULTS['BASE_STORAGE_DIRECTORY'] = data_storage_directory.strpath

    projects = [api.create_instance('projects', **factories.ProjectFactory())]
    workflows = [factories.WorkflowFactory(projects=projects) for i in range(3)]
    workflows = [api.create_instance('workflows', **i) for i in workflows]
    keys = [i['pk'] for i in workflows]

    importer = data.DataImporter()
    _, summary = importer.import_data(directories=[tmpdir.strpath], pk__in=keys)
    obtained = len(summary.rsplit('no files matched'))
    assert obtained == 3 + 1

    with pytest.raises(click.UsageError) as error:
        path_1 = tmpdir.join(f'{workflows[0]["system_id"]}.fastq')
        path_1.write('foo')
        importer.import_data(directories=[tmpdir.strpath], pk__in=keys)

    path_1.remove()
    assert 'cant determine if read 1 or read 2' in str(error.value)

    path_1 = tmpdir.join(f'{workflows[0]["system_id"]}_R1_foo.fastq')
    path_2 = tmpdir.join(f'{workflows[0]["system_id"]}_R2_foo.fastq')
    path_1.write('foo')
    path_2.write('foo')

    _, summary = importer.import_data(directories=[tmpdir.strpath], pk__in=keys, commit=True)
    assert 'samples matched: 1' in summary

    with pytest.raises(click.UsageError) as error:
        path_1 = tmpdir.join(f'{workflows[1]["system_id"]}_1.fastq')
        path_2 = tmpdir.join(f'{workflows[1]["system_id"]}.bam')
        path_1.write('foo')
        path_2.write('foo')
        importer.import_data(directories=[tmpdir.strpath], pk__in=keys, commit=True)

    path_1.remove()
    path_2.remove()
    assert 'multiple formats' in str(error.value)

    with pytest.raises(click.UsageError) as error:
        api.patch_instance('workflows', workflows[1]['pk'], center_id='dup_id')
        api.patch_instance('workflows', workflows[2]['pk'], center_id='dup_id')
        importer.import_data(
            key=lambda x: x['center_id'],
            directories=[tmpdir.strpath],
            pk__in=keys)

    assert 'same identifier for' in str(error.value)

    path_1 = tmpdir.join(f'_{workflows[1]["system_id"]}_cram1_.cram')
    path_2 = tmpdir.join(f'_{workflows[1]["system_id"]}_cram2_.cram')
    path_3 = tmpdir.join(f'_{workflows[2]["system_id"]}_bam1_.bam')
    path_4 = tmpdir.join(f'_{workflows[2]["system_id"]}_bam2_.bam')

    path_1.write('foo')
    path_2.write('foo')
    path_3.write('foo')
    path_4.write('foo')

    imported, summary = importer.import_data(
        directories=[tmpdir.strpath],
        commit=True,
        symlink=True,
        pk__in=keys)

    project = api.get_instance('projects', projects[0]['pk'])
    assert project['storage_url']
    assert imported[0]['storage_usage'] > 0
    assert imported[0]['sequencing_data']
    assert imported[1]['sequencing_data']
    assert 'workflows' in imported[1]['storage_url']
    assert len(os.listdir(os.path.join(imported[1]['storage_url'], 'data'))) == 2
    assert 'samples matched: 2' in summary
    assert 'samples skipped: 1' in summary

    api.patch_instance('workflows', workflows[1]['pk'], sequencing_data=None)
    file_data = tmpdir.join('file_data.yaml')

    with open(file_data.strpath, 'w') as f:
        yaml.dump({
            os.path.basename(path_1.strpath): {'PU': 'TEST_PU'},
            os.path.basename(path_2.strpath): {'PU': 'TEST_PU'},
            }, f, default_flow_style=False)

    command = data.DataImporter.as_cli_command()
    runner = CliRunner()
    args = [
        '-di', tmpdir.strpath, '-id', 'system_id', '-fi', 'pk__in', keys,
        '--files-data', file_data.strpath, '--commit']

    result = runner.invoke(command, args, catch_exceptions=False)
    assert 'samples matched: 1' in result.output
    workflows[1] = api.get_instance('workflows', workflows[1]['pk'])
    assert workflows[1]['sequencing_data'][0]['file_data']['PU'] == 'TEST_PU'
    assert workflows[1]['sequencing_data'][1]['file_data']['PU'] == 'TEST_PU'

    args = ['-di', tmpdir.strpath, '-id', 'specimen', '-fi', 'pk__in', keys]
    result = runner.invoke(command, args)
    assert 'invalid type for identifier' in result.output


def test_get_dst():
    importer = data.DataImporter()
    bam_test = ['sample.bam']
    cram_test = ['sample.cram']
    fastq_test = [
        ('sample_R{}_moretext', 'sample_moretext_{}'),
        ('sample_R{}_', 'sample_{}'),
        ('sample_R{}', 'sample_{}'),
        ('sample.R{}.more_text', 'sample_more_text_{}'),
        ('sample.R{}.', 'sample_{}'),
        ('sample.R{}', 'sample_{}'),
        ('sample_{}', 'sample_{}'),
        ]

    for i in bam_test:
        assert re.search(importer.BAM_REGEX, i)
        assert not re.search(importer.BAM_REGEX, i + 'not a bam')

    for i in cram_test:
        assert re.search(importer.CRAM_REGEX, i)
        assert not re.search(importer.CRAM_REGEX, i + 'not a cram')

    for test, expected in fastq_test:
        for index in [1, 2]:
            for fastq in ['.fastq', '.fq']:
                for gzipped in ['', '.gz']:
                    path = test.format(index) + fastq + gzipped
                    obtained = importer.format_fastq_name(path)
                    assert obtained == expected.format(index) + '.fastq' + gzipped
