"""Logic to interact with API."""

from itertools import islice
import time
import collections

from requests.packages.urllib3.exceptions import InsecureRequestWarning
import click
import requests

from cli.settings import system_settings
from cli.settings import user_settings
from cli import utils

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
            if key == 'fields' and 'pk' not in value:  # pk is required
                value = ','.join(value.split(',') + ['pk'])
            filters_dict[key] = value
        elif isinstance(value, collections.Iterable):
            is_in = key.endswith('__in') or key.endswith('__in!')
            value = list(map(str, value))
            filters_dict[key] = [','.join(value)] if is_in else value
        else:  # pragma: no cover
            raise click.UsageError(f'Invalid filter: {key}, {value}')

    return filters_dict


def iterate(url, **filters):
    """
    Iterate through a paginated API endpoint and yield instances.

    Arguments:
        url (str): API URL address.
        filters (dict): name, value pairs for API filtering.

    Returns:
        list: of objects in the 'results' key of the API response.
    """
    limit = int(filters.get('limit', 5000))  # default limit to 5000 per hit
    filters['limit'] = limit
    filters = process_api_filters(**filters)
    nexturl = url
    objects = []

    while nexturl:
        filters = filters if nexturl == url else {}  # nexturl includes params
        response = api_request('get', url=nexturl, params=filters)
        results = response.json()
        nexturl = results['next']
        objects += results['results']

        if len(objects) >= limit:  # pragma: no cover
            message = f'retrieved {len(objects)} ouf of {results["count"]}...'
            click.echo(message, err=True)

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

    if endpoint == 'analyses' and data.get('status') == 'SUCCEEDED':
        msg = 'storage_usage must be updated on success.'
        assert data.get('storage_usage') is not None, msg
        utils.check_admin()

    data = api_request('patch', url=url, json=data).json()

    if endpoint == 'analyses' and data.get('status'):
        _run_signals('analyses', data, system_settings.ON_STATUS_CHANGE)

    if endpoint == 'workflows' and data.get('sequencing_data'):
        _run_signals('workflows', data, system_settings.ON_DATA_IMPORT)

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
        verbose (bool): print to stderr how many instances will be retrieved.

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
        count = get_instances_count(endpoint, **filters)
        count += len(identifiers or [])
        ids_msg = ' at least ' if identifiers else ' '  # ids may be in filters
        click.echo(
            f'Retrieving{ids_msg}{count} from {endpoint} API endpoint...',
            err=True)

    if filters or identifiers is None:
        instances += iterate(url, **filters)
        keys = {i['pk'] for i in instances}

    for chunk in chunks(identifiers or [], 10000):
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


def patch_analyses_status(analyses, status):
    """
    Patch the `status` of multiple `analyses`.

    Arguments:
        analyses (list): dicts of analyses instances.
        status (str): status to be updated to.

    Raises:
        AssertionError: if status not in {'SUBMITTED', 'STAGED'}.

    Returns:
        list: of updated analyses.
    """
    data = {'ids': [], 'status': status, 'ran_by': system_settings.api_username}
    url = f'{system_settings.API_BASE_URL}/analyses/bulk_update/'
    assert status in {'SUBMITTED', 'STAGED'}, f'status not supported: {status}'

    for i in analyses:  # change locally, bulk_update doesn't return instances
        i['status'] = status
        i['ran_by'] = data.get('ran_by', i['ran_by'])
        data['ids'].append(i['pk'])

    api_request('patch', url=url, json=data)

    for i in analyses:
        _run_signals('analyses', i, system_settings.ON_STATUS_CHANGE)

    return analyses


def _run_signals(endpoint, instance, signals):
    errors = []
    on_failure = system_settings.ON_SIGNAL_FAILURE or (lambda **_: None)

    for signal in signals:
        try:
            signal(instance)
        except Exception as error:  # pragma: no cover pylint: disable=W0703
            errors.append(error)

            try:
                on_failure(endpoint, instance, signal, error)
            except Exception as on_failure_error:  # pylint: disable=W0703
                errors.append(on_failure_error)
