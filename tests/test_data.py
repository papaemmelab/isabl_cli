import pytest
import click

from cli import data
from cli import api
from cli import _DEFAULTS

from . import factories

def test_local_data_import(tmpdir):
    data_dir = tmpdir.mkdir('data_dir')
    _DEFAULTS['DATA_STORAGE_DIRECTORY'] = str(data_dir)

    keys = [1, 2, 3]
    workflows = api.get_instances('workflows', keys)

    if not workflows:
        workflows = [factories.WorkflowFactory() for i in range(3)]
        workflows = [api.create_instance('workflows', **i) for i in workflows]

    for i in workflows:
        api.patch_instance('workflows', i['pk'], data_type=None)

    importer = data.LocalDataImporter()
    importer.import_data(directories=[tmpdir.strpath], pk__in=keys)
    obtained = len(importer._get_summary().rsplit('no files matched'))
    assert obtained == 3 + 1

    # keys = [lambda x: x['system_id'], lambda x: x['specimen']['system_id']]
    # for key in keys:
    with pytest.raises(click.UsageError) as error:
        path = tmpdir.join(f'{workflows[0]["system_id"]}.fastq')
        path.write('foo')
        importer.import_data(directories=[tmpdir.strpath], pk__in=keys)

    path.remove()
    assert 'cant determine if read 1 or read 2' in str(error.value)

    read_1 = tmpdir.join(f'{workflows[0]["system_id"]}_R1_foo.fastq')
    read_2 = tmpdir.join(f'{workflows[0]["system_id"]}_R2_foo.fastq')
    read_1.write('foo')
    read_2.write('foo')

    importer.import_data(directories=[tmpdir.strpath], pk__in=keys, commit=True)


def test_get_dst():
    importer = data.LocalDataImporter()
    assert importer._get_fastq_dst('not a fastq') is None
    assert importer._get_bam_dst('not a bam') is None

    bam_test = [
        'sample.bam'
        'sample.bam.bai'
        'sample.bam.md5'
        ]

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
        assert importer._get_bam_dst(i) == i

    for test, expected in fastq_test:
        for index in [1, 2]:
            for fastq in ['.fastq', '.fq']:
                for gzipped in ['', '.gz']:
                    path = test.format(index) + fastq + gzipped
                    obtained = importer._get_fastq_dst(path)
                    assert obtained == expected.format(index) + '.fastq' + gzipped
