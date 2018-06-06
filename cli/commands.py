"""commands logic."""

from glob import glob
from os.path import join
import json
import shutil
import subprocess

import click

from cli import api
from cli import data
from cli import options
from cli import utils


@click.command()
def processed_finished():
    """Process and update finished analyses."""
    utils.check_admin()
    for i in api.get_instances('analyses', status='FINISHED'):
        dst = i['storage_url']
        src = dst + '__tmp'
        shutil.move(dst, src)
        command = utils.get_rsync_command(src, dst, chmod='a-w')
        subprocess.check_call(command, shell=True)
        api.patch_instance('analyses', i['pk'], status='SUCCEEDED')


@click.command()
@options.PRIMARY_KEY
@options.ANALYSIS_STATUS
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
@options.ENDPOINT
@options.FIELDS
@options.FILTERS
@options.NO_HEADERS
def get_fields(endpoint, field, filters, no_headers):
    """Get instances database fields."""
    result = [] if no_headers else ['\t'.join('.'.join(i) for i in field)]

    for i in api.get_instances(endpoint, verbose=True, **filters):
        values = []
        for j in field:
            value = i
            for k in j:
                try:
                    value = value.get(k, f'INVALID FIELD ({k})') or None
                except AttributeError:
                    value = f'INVALID FIELD ({k})'

            if isinstance(value, dict):
                value = json.dumps(value)

            values.append(str(value))
        result.append('\t'.join(values))
    click.echo('\n'.join(result).expandtabs(30))


@click.command()
@options.ENDPOINT
@options.FILE_PATTERN
@options.FILTERS
def get_paths(endpoint, pattern, filters):
    """Get paths from storage directories."""
    for i in api.get_instances(endpoint, verbose=True, **filters):
        if i['storage_url']:
            click.echo('\n'.join(glob(join(i['storage_url'], pattern))))


@click.command()
@options.ENDPOINT
@options.FILTERS
def get_dirs(endpoint, filters):
    """Get instances storage directory."""
    for i in api.get_instances(endpoint, verbose=True, **filters):
        if i['storage_url']:
            click.echo(i['storage_url'])

@click.command()
@options.PRIMARY_KEY
@options.BEDFILE
def import_bed(primary_key, bedfile):
    """Import a technique bedfile."""
    data.import_bed(primary_key, bedfile)
