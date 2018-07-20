"""commands logic."""

from glob import glob
from os.path import join

import click

from cli import api
from cli import options
from cli import utils
from cli.settings import import_from_string


@click.command()
@click.option('--project', help='primary key of project to merge by', type=int)
@click.option('--pipeline', help='analyses pipeline primary key', type=int)
def merge_project_analyses(project, pipeline):  # pragma: no cover
    """Merge analyses by project primary key."""
    project = api.get_instance('projects', project)
    pipeline = api.get_instance('pipelines', pipeline)
    pipeline = import_from_string(pipeline['pipeline_class'])()
    pipeline.run_project_merge(project)


@click.command()
@options.FILTERS
def processed_finished(filters):
    """Process and update finished analyses."""
    utils.check_admin()
    filters.update(status='FINISHED')

    for i in api.get_instances('analyses', **filters):
        api.patch_analysis_status(i, 'SUCCEEDED')


@click.command()
@options.FILTERS
def patch_results(filters):
    """Update the results field of many analyses."""
    utils.check_admin()

    for i in api.get_instances('analyses', **filters):
        pipeline = import_from_string(i['pipeline']['pipeline_class'])()
        results = pipeline.get_analysis_results(i)
        api.patch_instance('analyses', i['pk'], results=results)


@click.command()
@options.ANALYSIS_PRIMARY_KEY
@options.ANALYSIS_STATUS
def patch_status(key, status):
    """Patch status of a given analysis."""
    analysis = api.get_instance('analyses', key)
    api.patch_analysis_status(analysis, status)


@click.command()
@options.ENDPOINT
@options.FIELDS
@options.NULLABLE_FILTERS
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
@options.NULLABLE_FILTERS
def get_count(endpoint, filters):
    """Get count of database instances."""
    click.echo(api.get_instances_count(endpoint, **filters))


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
