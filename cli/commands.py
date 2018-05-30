"""commands logic."""

from getpass import getuser
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

    if status == 'SUCCEEDED':
        utils.check_admin()
        get_storage_url = system_settings.GET_STORAGE_DIRECTORY_FUNCTION
        storage_url = get_storage_url('analyses', analysis['pk'])
        src = analysis['storage_url']

        if storage_url != src:  # true when BASE_RUN_DIRECTORY is set
            command = utils.get_rsync_command(src, storage_url, chmod='a-w')
        elif analysis['ran_by'] == getuser():
            command = f'chmod -R a-w {storage_url}'
        else:  # need to copy if ran by other user than admin
            src = src + '__tmp'
            shutil.move(storage_url, src)
            command = utils.get_rsync_command(src, storage_url, chmod='a-w')

        subprocess.check_call(command, shell=True)

    api.patch_instance(
        endpoint='analyses',
        identifier=analysis['pk'],
        storage_url=storage_url,
        storage_usage=utils.get_tree_size(storage_url))
