"""commands logic."""

from glob import glob
from os.path import join
import shutil
import subprocess

import click

from cli import api
from cli import data
from cli import options
from cli import utils
from cli.settings import system_settings
from cli.settings import import_from_string


@click.command()
@click.option('--project', help='primary key of project to merge by', type=int)
@click.option('--pipeline', help='analyses pipeline primary key', type=int)
def merge_analyses(project, pipeline):
    """Merge analyses by project primary key."""
    project = api.get_instance('projects', project)
    pipeline = api.get_instance('pipelines', pipeline)
    analyses = api.get_instances(
        endpoint='analyses',
        pipeline=pipeline['pk'],
        project=project['pk'],
        status='SUCCEEDED')

    if analyses:
        pipeline_class = import_from_string(pipeline['pipeline_class'])()
        merge_function = pipeline_class.merge_analyses_by_project

        try:
            merge_function.__isabstractmethod__  # pylint: disable=W0104
            raise NotImplementedError('No merging logic available!')
        except AttributeError:
            pass

        analysis = api.create_instance(
            endpoint='analyses',
            project_level_analysis=project,
            pipeline=pipeline)

        if not analysis['storage_url']:
            analysis = data.update_storage_url(
                endpoint='analyses',
                identifier=analysis['pk'],
                use_hash=True)

        try:
            merge_function(analysis['storage_url'], analyses)
        except Exception as error:
            api.patch_instance('analyses', analysis['pk'], status='FAILED')
            raise error


@click.command()
@options.FILTERS
def processed_finished(filters):
    """Process and update finished analyses."""
    utils.check_admin()
    filters.update(status='FINISHED')

    for i in api.get_instances('analyses', **filters):
        if i['ran_by'] != system_settings.api_username:  # admin must own dir
            src = i['storage_url'] + '__tmp'
            shutil.move(i['storage_url'], src)
            cmd = utils.get_rsync_command(src, i['storage_url'], chmod='a-w')
            subprocess.check_call(cmd, shell=True)

        api.patch_instance(
            endpoint='analyses',
            identifier=i['pk'],
            status='SUCCEEDED',
            storage_usage=utils.get_tree_size(i['storage_url']))


@click.command()
@options.ANALYSIS_PRIMARY_KEY
@options.ANALYSIS_STATUS
def patch_status(key, status):
    """Patch status of a given analysis."""
    analysis = api.get_instance('analyses', key)
    api.patch_instance(
        endpoint='analyses',
        identifier=analysis['pk'],
        status=status,
        storage_usage=utils.get_tree_size(analysis['storage_url']))


@click.command()
@options.ENDPOINT
@options.FIELDS
@options.FILTERS
@options.NO_HEADERS
def get_attributes(endpoint, field, filters, no_headers):
    """Get database attributes from API."""
    result = [] if no_headers else ['\t'.join('.'.join(i) for i in field)]
    filters['fields'] = ','.join(i[0] for i in field)

    for i in api.get_instances(endpoint, verbose=True, **filters):
        values = [utils.traverse_dict(i, j, serialize=True) for j in field]
        result.append('\t'.join(values))

    click.echo('\n'.join(result).expandtabs(30))


@click.command()
@options.ENDPOINT
@options.FILE_PATTERN
@options.FILTERS
def get_paths(endpoint, pattern, filters):
    """Get storage directories, use `pattern` to match files inside dirs."""
    filters.update(fields='storage_url', limit=100000)
    for i in api.get_instances(endpoint, verbose=True, **filters):
        if i['storage_url']:
            if pattern:
                click.echo('\n'.join(glob(join(i['storage_url'], pattern))))
            else:
                click.echo(i['storage_url'])
