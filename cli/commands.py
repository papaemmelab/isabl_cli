"""commands logic."""

import subprocess
import shutil
import click

from cli import system_settings
from cli import utils
from cli import api


@click.command()
@click.option(
    "--key",
    show_default=True,
    type=click.INT,
    required=True)
@click.option(
    "--status",
    show_default=True,
    type=click.Choice([
        'CREATED', 'FAILED', 'FINISHED', 'IN_PROGRESS',
        'STAGED', 'STARTED', 'SUBMITTED', 'SUCCEEDED']),
    required=True)
def patch_status(key, status):
    """Patch status of a given analysis."""
    analysis = api.get_instance('analyses', key)
    storage_url = analysis['storage_url']

    if analysis['status'] == 'SUCCEEDED':
        raise click.UsageError('Analysis already SUCCEEDED, cannot be patched')

    if status == 'SUCCEEDED' and True:
        utils.check_admin()

    api.patch_instance(
        endpoint='analyses',
        identifier=analysis['pk'],
        status=status,
        storage_url=storage_url,
        storage_usage=utils.get_tree_size(storage_url))


def processed_finished():
    utils.check_admin()

    for i in api.get_instances('analyses', status='FINISHED'):
        get_storage_url = system_settings.GET_STORAGE_DIRECTORY_FUNCTION
        storage_url = get_storage_url('analyses', i['pk'])
        src = i['storage_url']

        if storage_url != src:  # true when BASE_RUN_DIRECTORY is set
            command = utils.get_rsync_command(src, storage_url, chmod='a-w')
        elif i['ran_by'] == system_settings.api_username:
            command = f'chmod -R a-w {storage_url}'
        else:  # need to copy if ran by other user than admin
            src = src + '__tmp'
            shutil.move(storage_url, src)
            command = utils.get_rsync_command(src, storage_url, chmod='a-w')

        subprocess.check_call(command, shell=True)

        api.patch_instance(
            endpoint='analyses',
            identifier=i['pk'],
            storage_url=storage_url,
            status='SUCCEEDED',
            storage_usage=utils.get_tree_size(storage_url))
