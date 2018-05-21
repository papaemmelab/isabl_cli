from click.testing import CliRunner
import click
import pytest
import requests

from cli import api
from cli import user_settings
from cli import system_settings

from . import factories


def test_get_token_headers():
    """Sample test for main command."""
    user_settings.token = None

    @click.command()
    def test():
        assert 'Authorization' in api.get_token_headers()

    runner = CliRunner()
    result = runner.invoke(test, input='admin\nadmin')
    assert not result.exception


def test_api_request():
    assert api.api_request('get', url=f'{system_settings.API_BASE_URL}/auth').ok
    with pytest.raises(requests.HTTPError):
        api.api_request('delete', url='http://google.com')


def test_api_methods():
    endpoint = 'diseases'
    diseases = [factories.DiseaseFactory() for _ in range(3)]
    created = [api.create_instance(endpoint, **i) for i in diseases]
    pk = created[0]['pk']
    pks = [i['pk'] for i in created[:2]]
    patched = api.patch_instance(endpoint, pk, data={'one': 1})

    assert patched['data']['one'] == 1
    assert api.get_instance(endpoint, pk)['pk'] == pk
    assert api.get_instances(endpoint, pk=pk)[0]['pk'] == pk
    assert api.get_instances_count(endpoint, pk=pk) == 1
    assert len(api.get_instances(endpoint)) == api.get_instances_count(endpoint)
    assert len(api.get_instances(endpoint, pks)) == 2
    assert len(api.get_instances(endpoint, pk__in=pks)) == 2

    for i in created:
        assert api.delete_instance(endpoint, i['pk']) is None


def test_patch_analyses_status():
    pipeline = factories.PipelineFactory()
    analyses = [factories.AnalysisFactory(pipeline=pipeline) for _ in range(2)]
    created = [api.create_instance('analyses', **i) for i in analyses]
    assert all([i['status'] == 'CREATED' for i in created])

    pks = [i['pk'] for i in created]
    api.patch_analyses_status(pks, 'SUCCEEDED')
    retrieved = api.get_instances('analyses', pks)
    assert all([i['status'] == 'SUCCEEDED' for i in retrieved])

    for i in created:
        api.delete_instance('analyses', i['pk'])

    api.delete_instance('pipelines', created[0]['pipeline']['pk'])


def test_system_id():
    data_a = factories.WorkflowFactory()
    data_b = factories.WorkflowFactory(specimen=data_a['specimen'])
    instance_a = api.create_instance('workflows', **data_a)
    instance_b = api.create_instance('workflows', **data_b)
    system_ids = [instance_a['system_id'], instance_b['system_id']]
    assert instance_a['specimen']['pk'] == instance_b['specimen']['pk']
    assert api.get_instance('workflows', system_ids[0])['pk'] == instance_a['pk']
    assert len(api.get_instances('workflows', system_ids)) == 2

    instance_a['specimen']['data']['key'] = 'value'
    instance_a['specimen']['notes'] = 'a note'
    patched = api.patch_instance('workflows', instance_a['pk'], specimen=instance_a['specimen'])
    assert patched['specimen']['data']['key'] == 'value'
    assert patched['specimen']['notes'] == 'a note'
