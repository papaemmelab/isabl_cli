"""Logic to interact with API."""

import time

from requests.packages.urllib3.exceptions import InsecureRequestWarning
import click
import requests

from cli import system_settings
from cli import user_settings

requests.packages.urllib3.disable_warnings(InsecureRequestWarning)  # pylint: disable=E1101


def get_token_headers():
    """Get an API token and store it in user's home directory."""
    token_url = system_settings.API_TOKEN_URL
    auth_url = f'{system_settings.API_V1_URL}/auth'
    headers = {'Authorization': f'Token {user_settings.token}'}
    response = requests.get(url=auth_url, headers=headers)

    if not response.ok:
        data = {
            "username": click.prompt('username', type=str, hide_input=False),
            "password": click.prompt('password', type=str, hide_input=True)
            }

        response = requests.post(url=token_url, data=data)
        response.raise_for_status()
        user_settings.token = response.json()['token']
        headers = {'Authorization': f'Token {user_settings.token}'}

    return headers


def api_request(method, **kwargs):
    """Perform any request operation using a naive retry implementation."""
    kwargs['headers'] = get_token_headers()

    for i in [0.2, 0.4, 0.6, 0.8, 5, 10]:  # attempt some retries
        response = getattr(requests, method)(verify=False, **kwargs)

        if not str(response.status_code).startswith('50'):
            break
        else:
            time.sleep(i)

    if not response.ok:
        try:
            msg = click.style(str(response.json()), fg='red', blink=True)
        except Exception:  # pylint: disable=broad-except
            msg = ''

        click.echo(f'API error: {response.url}\n{msg}')
        response.raise_for_status()

    return response
