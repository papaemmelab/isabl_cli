"""commands logic."""

import subprocess
import shutil
import click

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

    api.patch_instance(
        endpoint='analyses',
        identifier=analysis['pk'],
        status=status,
        storage_url=storage_url,
        storage_usage=utils.get_tree_size(storage_url))


@click.command()
def processed_finished():
    """Process and update all finished analyses that were not ran by admin."""
    utils.check_admin()
    for i in api.get_instances('analyses', status='FINISHED'):
        dst = i['storage_url']
        src = dst + '__tmp'
        shutil.move(dst, src)
        command = utils.get_rsync_command(src, dst, chmod='a-w')
        subprocess.check_call(command, shell=True)
        api.patch_instance('analyses', i['pk'], status='SUCCEEDED')
