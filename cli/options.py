"""Click options decorators."""

import click


IDENTIFIER = click.option(
    '--identifier', '-id',
    show_default=True,
    help='identifier to be used (traverse with dot, e.g. `specimen.system_id`)',
    callback=lambda _, __, i: i.split('.'),
    required=True)

FIELDS = click.option(
    '--field', '-f',
    show_default=True,
    multiple=True,
    help='fields to be retrieved (traverse with dot, e.g. `specimen.disease`)',
    callback=lambda _, __, i: [j.split('.') for j in i],
    required=True)

ANALYSIS_STATUS = click.option(
    '--status',
    show_default=True,
    help='analysis status',
    type=click.Choice([
        'CREATED', 'FAILED', 'FINISHED', 'IN_PROGRESS',
        'STAGED', 'STARTED', 'SUBMITTED', 'SUCCEEDED']),
    required=True)

ENDPOINT = click.option(
    '--endpoint', '-e',
    help='API endpoint to be queried',
    default='analyses',
    show_default=True,
    type=click.Choice([
        'workflows',
        'analyses',
        'projects',
        'techniques',
        'specimens',
        'individuals',
        'assemblies']),
    required=True)

ANALYSIS_PRIMARY_KEY = click.option(
    '--key',
    help='analysis primary key',
    show_default=True,
    type=click.INT,
    required=True)

TECHNIQUE_PRIMARY_KEY = click.option(
    '--key',
    help='technique primary key',
    show_default=True,
    type=click.INT,
    required=True)

COMMIT = click.option(
    '--commit',
    help='commit results',
    show_default=True,
    is_flag=True)

NO_HEADERS = click.option(
    '--no-headers',
    help='do not output headers',
    show_default=True,
    is_flag=True)

FILE_PATTERN = click.option(
    '--pattern',
    help='paths pattern to be matched',
    default=None,
    show_default=True,
    type=click.STRING,
    required=False)

FILTERS = click.option(
    '--filters', '-fi',
    multiple=True,
    help='API filters',
    show_default=True,
    type=(str, str),
    callback=lambda _, __, i: dict(i),
    required=True)

DIRECTORIES = click.option(
    '--directories',
    '-di',
    show_default=True,
    type=click.Path(
        exists=True,
        file_okay=False,
        dir_okay=True,
        writable=False,
        readable=True,
        resolve_path=True),
    multiple=True,
    help='directories to be processed')

SYMLINK = click.option(
    '--symlink',
    help='use symlink',
    show_default=True,
    is_flag=True)

REFERENCE_DATA_SOURCE = click.option(
    '--data-src',
    show_default=True,
    type=click.Path(
        resolve_path=True,
        exists=True,
        readable=True),
    help='path to reference data')

TARGETS_PATH = click.option(
    '--targets-path',
    show_default=True,
    type=click.Path(
        resolve_path=True,
        exists=True,
        file_okay=True,
        dir_okay=False,
        writable=False,
        readable=True),
    help='path to targets bedfile')

BAITS_PATH = click.option(
    '--baits-path',
    show_default=True,
    type=click.Path(
        resolve_path=True,
        exists=True,
        file_okay=True,
        dir_okay=False,
        writable=False,
        readable=True),
    help='path to baits bedfile')

FILES_DATA = click.option(
    '--files-data',
    help='a yaml file with extra annotation for imported files.',
    show_default=True,
    default=None,
    type=click.Path(
        resolve_path=True,
        exists=True,
        file_okay=True,
        dir_okay=False,
        writable=False,
        readable=True))
