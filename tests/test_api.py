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
    """Sample test for main command."""
    assert api.api_request('get', url=f'{system_settings.API_V1_URL}/auth').ok

    with pytest.raises(requests.HTTPError):
        api.api_request('post', url='http://google.com')


def test_get_instances():
    endpoint = 'diseases'
    created = []
    diseases = [factories.DiseaseFactory() for _ in range(3)]
    url = f'{system_settings.API_V1_URL}/{endpoint}'

    for i in diseases:
        response = api.api_request('post', url=url, data=i)
        assert response.ok
        created.append(response.json())

    pk = created[0]['pk']
    pks = [i['pk'] for i in created[:2]]

    assert api.iterate(url, pk=pk)[0]['pk'] == pk
    assert [i['pk'] for i in api.iterate(url, pk__in=pks)] == pks
    assert api.get_instance(endpoint, pk)['pk'] == pk
    assert api.get_instances_count(endpoint, pk=pk) == 1
    assert len(api.get_instances(endpoint)) == api.get_instances_count(endpoint)
    assert len(api.get_instances(endpoint, pks)) == 2

    for i in created:
        assert api.api_request('delete', url=f'{url}/{i["pk"]}').ok
