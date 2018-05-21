"""Logic to interact with API."""

from itertools import islice
import time
import collections

from requests.packages.urllib3.exceptions import InsecureRequestWarning
import click
import requests

from cli import system_settings
from cli import user_settings

requests.packages.urllib3.disable_warnings(InsecureRequestWarning)  # pylint: disable=E1101


def chunks(array, size):
    """Yield successive n-sized chunks from ``array``."""
    array = iter(array)
    return iter(lambda: tuple(islice(array, size)), ())


def get_token_headers():
    """Get an API token and store it in user's home directory."""
    auth_url = f'{system_settings.API_BASE_URL}/auth'
    headers = {'Authorization': f'Token {user_settings.API_TOKEN}'}
    response = requests.get(url=auth_url, headers=headers)

    if not response.ok:
        data = {
            "username": click.prompt('username', type=str, hide_input=False),
            "password": click.prompt('password', type=str, hide_input=True)
            }

        response = requests.post(url=f'{auth_url}/token', data=data)
        response.raise_for_status()
        user_settings.API_TOKEN = response.json()['token']
        headers = {'Authorization': f'Token {user_settings.API_TOKEN}'}

    return headers


def api_request(method, **kwargs):
    """Perform any request operation using a naive retry implementation."""
    kwargs['headers'] = get_token_headers()
    kwargs['params'] = kwargs.get('params', {})
    kwargs['params']['format'] = 'json'

    for i in [0.2, 0.4, 0.6, 0.8, 5, 10]:  # attempt some retries
        response = getattr(requests, method)(verify=False, **kwargs)

        if not str(response.status_code).startswith('50'):
            break
        else:  # pragma: no cover
            time.sleep(i)

    if not response.ok:
        try:
            msg = click.style(str(response.json()), fg='red', blink=True)
        except Exception:  # pylint: disable=broad-except
            msg = ''

        click.echo(f'Request Error: {response.url}\n{msg}')
        response.raise_for_status()

    return response


def process_api_filters(**filters):
    """
    Process filters for API.

    Arguments:
        filters (dict): name, value pairs for API filtering.

    Returns:
        dict: a `requests` formatted `params` dict.
    """
    filters = filters or {}
    filters_dict = {}

    for key, value in filters.items():
        if isinstance(value, (str, int, float, type(None))):
            filters_dict[key] = value
        elif isinstance(value, collections.Iterable):
            is_in = key.endswith('__in') or key.endswith('__in!')
            value = list(map(str, value))
            filters_dict[key] = [','.join(value)] if is_in else value
        else:  # pragma: no cover
            raise click.UsageError(f'Invalid filter: {key}, {value}')

    return filters_dict


def iterate(url, limit=1000, **filters):
    """
    Iterate through a paginated API endpoint and yield instances.

    Arguments:
        url (str): API URL address.
        filters (dict): name, value pairs for API filtering.
        limit (int): max number of instances requested to API at the same time.

    Returns:
        list: of objects in the 'results' key of the API response.
    """
    filters = process_api_filters(limit=limit, **filters)
    nexturl = url
    objects = []

    while nexturl:
        filters = filters if nexturl == url else None  # nexturl includes params
        response = api_request('get', url=nexturl, params=filters)
        results = response.json()
        nexturl = results['next']
        objects += results['results']

    return objects


def get_instance(endpoint, identifier):
    """
    Get database instance.

    Arguments:
        identifier (str): a primary key, system_id, email or username.
        endpoint (str): endpoint without API base URL (e.g. `analyses`).

    Returns:
        types.SimpleNamespace: loaded with data returned from the API.
    """
    url = f'{system_settings.API_BASE_URL}/{endpoint}/{identifier}'
    data = api_request('get', url=url).json()
    return data


def create_instance(endpoint, **data):
    """
    Create database instance.

    Arguments:
        endpoint (str): endpoint without API base URL (e.g. `analyses`).
        data (dict): fields to be created.

    Returns:
        types.SimpleNamespace: loaded with data returned from the API.
    """
    url = f'{system_settings.API_BASE_URL}/{endpoint}'
    data = api_request('post', url=url, json=data).json()
    return data


def patch_instance(endpoint, identifier, **data):
    """
    Patch database instance.

    Arguments:
        identifier (str): a primary key, system_id, email or username.
        endpoint (str): endpoint without API base URL (e.g. `analyses`).
        data (dict): fields to be patched.

    Returns:
        types.SimpleNamespace: loaded with data returned from the API.
    """
    url = f'{system_settings.API_BASE_URL}/{endpoint}/{identifier}'
    data = api_request('patch', url=url, json=data).json()
    return data


def delete_instance(endpoint, identifier):
    """
    Delete database instance.

    Arguments:
        identifier (str): a primary key, system_id, email or username.
        endpoint (str): endpoint without API base URL (e.g. `analyses`).
    """
    url = f'{system_settings.API_BASE_URL}/{endpoint}/{identifier}'
    api_request('delete', url=url)


def get_instances(endpoint, identifiers=None, verbose=False, **filters):
    """
    Return instances from a list API endpoint.

    if not `identifiers` and not `filters` retrieves all objects in database.

    Arguments:
        endpoint (str): endpoint without API base URL (e.g. `analyses`).
        identifiers (list): List of identifiers.
        filters (dict): name, value pairs for API filtering.

    Raises:
        click.UsageError: if string identifier and endpoint not in individuals,
            specimens or workdflows.

    Returns:
        list: of types.SimpleNamespace objects loaded with dicts from API.
    """
    check_system_id = endpoint in {'individuals', 'specimens', 'workflows'}
    url = f'{system_settings.API_BASE_URL}/{endpoint}'
    instances = []
    keys = set()

    if verbose:
        count = get_instances_count(**filters) + len(identifiers or [])
        click.echo(f'Retrieving at least {count} {endpoint}...')

    if filters or identifiers is None:
        instances += iterate(url, limit=2000, **filters)
        keys = {i['pk'] for i in instances}

    for chunk in chunks(identifiers or [], 1000):
        primary_keys = set()
        system_ids = set()

        for i in map(str, chunk):
            if i.isdigit():
                primary_keys.add(i)
            elif check_system_id:
                system_ids.add(i)
            else:  # pragma: no cover
                msg = f'msg invalid identifier for {endpoint}: {i}'
                raise click.UsageError(msg)

        if primary_keys:
            kwargs = {'pk__in': ','.join(primary_keys), 'url': url}
            instances += [i for i in iterate(**kwargs) if i['pk'] not in keys]

        if system_ids:
            kwargs = {'system_id__in': ','.join(system_ids), 'url': url}
            instances += [i for i in iterate(**kwargs) if i['pk'] not in keys]

    return instances


def get_instances_count(endpoint, **filters):
    """
    Return the count of a list url.

    Arguments:
        endpoint (str): endpoint without API base URL (e.g. `analyses`).
        filters (dict): name, value pairs for API filtering.

    Returns:
        int: count of objects matching the provided filters.
    """
    url = f'{system_settings.API_BASE_URL}/{endpoint}'
    filters = process_api_filters(**filters)
    filters['limit'] = 1
    return int(api_request('get', url=url, params=filters).json()['count'])


def patch_analyses_status(primary_keys, status):
    """Patch the `status` of multiple analyses given their `primary_keys`."""
    url = f'{system_settings.API_BASE_URL}/analyses/status/'
    api_request('patch', url=url, json={'ids': primary_keys, 'status': status})
