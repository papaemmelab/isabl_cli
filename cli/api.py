"""Logic to interact with API."""

from itertools import islice
import collections
import time
import shutil
import subprocess

from requests.packages.urllib3.exceptions import InsecureRequestWarning
import click
import requests

from cli import utils
from cli.settings import system_settings
from cli.settings import user_settings
from cli.settings import import_from_string

requests.packages.urllib3.disable_warnings(  # pylint: disable=E1101
    InsecureRequestWarning
)


def chunks(array, size):
    """Yield successive n-sized chunks from ``array``."""
    array = iter(array)
    return iter(lambda: tuple(islice(array, size)), ())


def get_token_headers():
    """Get an API token and store it in user's home directory."""
    url = f"{system_settings.API_BASE_URL}/auth"
    headers = {"Authorization": f"Token {user_settings.API_TOKEN}"}
    response = requests.get(url=url, headers=headers)

    if not response.ok:
        data = {
            "username": click.prompt("username", type=str, hide_input=False),
            "password": click.prompt("password", type=str, hide_input=True),
        }

        auth_url = f"{system_settings.API_BASE_URL}/rest-auth/login/"
        response = requests.post(url=f"{auth_url}", data=data)
        response.raise_for_status()
        user_settings.API_TOKEN = response.json()["key"]
        headers = {"Authorization": f"Token {user_settings.API_TOKEN}"}

    return headers


def api_request(method, **kwargs):
    """Perform any request operation using a naive retry implementation."""
    kwargs["headers"] = get_token_headers()
    kwargs["params"] = kwargs.get("params", {})
    kwargs["params"]["format"] = "json"

    for i in [0.2, 0.4, 0.6, 0.8, 5, 10]:  # attempt some retries
        response = getattr(requests, method)(verify=False, **kwargs)

        if not str(response.status_code).startswith("50"):
            break
        else:  # pragma: no cover
            time.sleep(i)

    if not response.ok:
        try:
            msg = click.style(str(response.json()), fg="red", blink=True)
        except Exception:  # pylint: disable=broad-except
            msg = ""

        click.echo(f"Request Error: {response.url}\n{msg}")
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
            if key == "fields" and "pk" not in value:  # pk is required
                value = ",".join(value.split(",") + ["pk"])
            filters_dict[key] = value
        elif isinstance(value, collections.Iterable):
            is_in = key.endswith("__in") or key.endswith("__in!")
            value = list(map(str, value))
            filters_dict[key] = [",".join(value)] if is_in else value
        else:  # pragma: no cover
            raise click.UsageError(f"Invalid filter: {key}, {value}")

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
    limit = int(filters.get("limit", 5000))  # default limit to 5000 per hit
    filters["limit"] = limit
    filters = process_api_filters(**filters)
    nexturl = url
    objects = []

    while nexturl:
        filters = filters if nexturl == url else {}  # nexturl includes params
        response = api_request("get", url=nexturl, params=filters)
        results = response.json()
        nexturl = results["next"]
        objects += results["results"]

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
    url = f"{system_settings.API_BASE_URL}/{endpoint}/{identifier}"
    return api_request("get", url=url).json()


def create_instance(endpoint, **data):
    """
    Create database instance.

    Arguments:
        endpoint (str): endpoint without API base URL (e.g. `analyses`).
        data (dict): fields to be created.

    Returns:
        types.SimpleNamespace: loaded with data returned from the API.
    """
    url = f"{system_settings.API_BASE_URL}/{endpoint}"
    return api_request("post", url=url, json=data).json()


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
    url = f"{system_settings.API_BASE_URL}/{endpoint}/{identifier}"
    instance = api_request("patch", url=url, json=data).json()

    if endpoint == "analyses" and instance.get("status"):
        _run_signals("analyses", instance, system_settings.ON_STATUS_CHANGE)

    if endpoint == "experiments" and instance.get("sequencing_data"):
        _run_signals("experiments", instance, system_settings.ON_DATA_IMPORT)

    return instance


def delete_instance(endpoint, identifier):
    """
    Delete database instance.

    Arguments:
        identifier (str): a primary key, system_id, email or username.
        endpoint (str): endpoint without API base URL (e.g. `analyses`).
    """
    url = f"{system_settings.API_BASE_URL}/{endpoint}/{identifier}"
    api_request("delete", url=url)


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
            samples or workdflows.

    Returns:
        list: of types.SimpleNamespace objects loaded with dicts from API.
    """
    check_system_id = endpoint in {"individuals", "samples", "experiments"}
    url = f"{system_settings.API_BASE_URL}/{endpoint}"
    instances = []
    keys = set()

    if verbose:
        count = get_instances_count(endpoint, **filters)
        count += len(identifiers or [])
        ids_msg = " at least " if identifiers else " "  # ids may be in filters
        count = f"Retrieving{ids_msg}{count} from {endpoint} API endpoint..."
        click.echo(count, err=True)

    if filters or identifiers is None:
        instances += iterate(url, **filters)
        keys = {i["pk"] for i in instances if i.get("pk")}

    for chunk in chunks(identifiers or [], 10000):
        primary_keys = set()
        system_ids = set()

        for i in map(str, chunk):
            if i.isdigit():
                primary_keys.add(i)
            elif check_system_id:
                system_ids.add(i)
            else:  # pragma: no cover
                msg = f"msg invalid identifier for {endpoint}: {i}"
                raise click.UsageError(msg)

        if primary_keys:
            kwargs = {"pk__in": ",".join(primary_keys), "url": url}
            instances += [i for i in iterate(**kwargs) if i["pk"] not in keys]

        if system_ids:
            kwargs = {"system_id__in": ",".join(system_ids), "url": url}
            instances += [i for i in iterate(**kwargs) if i["pk"] not in keys]

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
    url = f"{system_settings.API_BASE_URL}/{endpoint}"
    filters = process_api_filters(**filters)
    filters["limit"] = 1
    return int(api_request("get", url=url, params=filters).json()["count"])


def get_tree(identifier):
    """Get everything for an individual."""
    return get_instance("individuals/tree", identifier)


def get_trees(identifiers=None, **filters):
    """Get everything for multiple individuals."""
    return get_instances("individuals/tree", identifiers, **filters)


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
    data = {"ids": [], "status": status, "ran_by": system_settings.api_username}
    url = f"{system_settings.API_BASE_URL}/analyses/bulk_update/"
    assert status in {"SUBMITTED", "STAGED"}, f"status not supported: {status}"

    for i in analyses:  # change locally, bulk_update doesn't return instances
        i["status"] = status
        i["ran_by"] = data.get("ran_by", i["ran_by"])
        data["ids"].append(i["pk"])

    api_request("patch", url=url, json=data)

    for i in analyses:
        _run_signals("analyses", i, system_settings.ON_STATUS_CHANGE)

    return analyses


def patch_analysis_status(analysis, status):
    """
    Patch a successful analysis.

    Make sure analysis is owned by admin user and that results field is updated.

    Arguments:
        analysis (dict): analysis instance.
        status (dict): analysis status.

    Returns:
        dict: patched analysis instance.
    """
    data = {"status": status}
    application = analysis["application"]
    storage_url = analysis["storage_url"]

    if status in {"FAILED", "SUCCEEDED", "IN_PROGRESS"}:
        data["storage_usage"] = utils.get_tree_size(storage_url)

    # admin must own directory of regular analyses
    if status == "SUCCEEDED" and not analysis["project_level_analysis"]:
        utils.check_admin()

        if analysis["ran_by"] != system_settings.api_username:
            src = storage_url + "__tmp"
            shutil.move(storage_url, src)
            cmd = utils.get_rsync_command(src, storage_url, chmod="a-w")
            subprocess.check_call(cmd, shell=True)

    if status in {"SUCCEEDED", "IN_PROGRESS"}:
        try:
            application = import_from_string(application["application_class"])()
            get_results = application.get_analysis_results

            if analysis["project_level_analysis"]:
                get_results = application.get_project_analysis_results

            data["results"] = get_results(analysis)
        except ImportError:
            pass
        except Exception as error:  # pragma: no cover
            data["status"] = "FAILED"
            patch_instance("analyses", analysis["pk"], **data)
            raise error

    return patch_instance("analyses", analysis["pk"], **data)


def _run_signals(endpoint, instance, signals):
    errors = []
    on_failure = system_settings.ON_SIGNAL_FAILURE or (lambda *_, **__: None)

    for signal in signals:
        try:
            signal(instance)
        except Exception as error:  # pragma: no cover pylint: disable=W0703
            errors.append(error)

            try:
                on_failure(endpoint, instance, signal, error)
            except Exception as on_failure_error:  # pylint: disable=W0703
                errors.append(on_failure_error)

    # if errors:
    # errors = '\n'.join(click.style(str(i), fg='red') for i in errors)
    # raise RuntimeError('Errors occurred during signals run:\n' + errors)
