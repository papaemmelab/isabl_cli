"""Data import logic."""

from datetime import datetime
from getpass import getuser
from os.path import basename
from os.path import isdir
from os.path import join
import os
import re
import shutil
import subprocess

import click

from cli import api
from cli import options
from cli import system_settings
from cli import utils


def symlink_workflow_to_projects(workflow):
    """Create symlink from workflow directory and projects directories."""
    for i in workflow['projects']:
        storage_url = i['storage_url']

        if not storage_url:
            get_dir = system_settings.GET_STORAGE_DIRECTORY
            storage_url = get_dir('projects', i['pk'], use_hash=False)
            os.makedirs(storage_url, exist_ok=True)
            api.patch_instance('projects', i['pk'], storage_url=storage_url)

        utils.force_symlink(
            workflow['storage_url'],
            join(storage_url, workflow['system_id']))


def symlink_analysis_to_targets(analysis):
    """Create symlink from workflow directory and projects directories."""
    for i in analysis['targets']:
        storage_url = i['storage_url']

        if not storage_url:
            get_dir = system_settings.GET_STORAGE_DIRECTORY
            storage_url = get_dir('workflows', i['pk'])
            os.makedirs(storage_url, exist_ok=True)
            api.patch_instance('workflows', i['pk'], storage_url=storage_url)

        dst = '__'.join([
            analysis['pipeline']['name'].lower().replace(' ', '_'),
            analysis['pipeline']['version'].lower().replace(' ', '_'),
            str(analysis['pk'])])

        utils.force_symlink(analysis['storage_url'], join(storage_url, dst))


def trash_analysis_storage(analysis):
    """Move analysis `storage_url` to a trash directory."""
    if analysis['status'] == 'SUCCEEDED':
        raise click.UsageError("You can't wipe a succeeded analysis")

    if isdir(analysis['storage_url']):
        slug = f'primary_key_{analysis["pk"]}__user_{getuser()}__date_'
        slug += datetime.now(system_settings.TIME_ZONE).isoformat()
        trash_dir = get_storage_directory('.analyses_trash', analysis['pk'])
        os.makedirs(trash_dir, exist_ok=True)
        dst = join(trash_dir, slug)
        click.echo(f"\ntrashing: {analysis['storage_url']} -> {dst}\n")
        shutil.move(analysis['storage_url'], dst)


def get_storage_directory(endpoint, primary_key, root=None, use_hash=True):
    """
    Get path to instance's data directory.

    If `use_hash` a naming system using the primay key is used. For example
    given `endpoint` analyses with primary keys 12345 and 2345:

        {root}/analyses/23/45/2345
        {root}/analyses/23/45/12345

    If `use_hash` is False:

        {root}/projects/100
        {root}/techniques/2

    Arguments:
        endpoint (str): instance's API endpoint.
        primary_key (str): instance's primary key.
        root (str): default is system_settings.BASE_STORAGE_DIRECTORY.
        use_hash (bool): hash primary key for directories.

    Returns:
        str: path to instance's data directory.
    """
    root = root or system_settings.BASE_STORAGE_DIRECTORY
    hash_1 = f'{primary_key:04d}' [-4:-2]
    hash_2 = f'{primary_key:04d}' [-2:]

    if not root:  # pragma: no cover
        raise click.UsageError('Setting `BASE_STORAGE_DIRECTORY` not defined.')

    if use_hash:
        path = os.path.join(endpoint, hash_1, hash_2)
    else:
        path = os.path.join(endpoint)

    return os.path.join(root, path, str(primary_key))


def import_bed(technique_primary_key, input_bed_path):
    """
    Register input_bed_path in technique data directory.

    This method will sort the bed and compress + tabix it. Both gzipped and
    uncompressed versions are kept.

    Instance's `storage_url`, `storage_usage`, `bed_url` are updated, setting
    the latter to the path to the uncompressed bedfile.

    Arguments:
        technique_primary_key (int): technique primary key.
        input_bed_path (str): path to incoming bedfile.

    Returns:
        dict: updated technique instance as retrieved from API.
    """
    if not input_bed_path.endswith('.bed'):  # pragma: no cover
        raise click.UsageError(f'No .bed suffix: {input_bed_path}')

    instance = api.get_instance('techniques', technique_primary_key)
    data_dir_fn = system_settings.GET_STORAGE_DIRECTORY
    data_dir = data_dir_fn('techniques', instance['pk'])

    if instance['bed_url']:
        raise click.UsageError(
            f'{instance["slug"]} has a bed registered: {instance["bed_url"]}')

    os.makedirs(data_dir, exist_ok=True)
    bed_path = join(data_dir, f"{instance['slug']}.bed")
    sorted_bed = subprocess.check_output(
        ['sort', '-k1,1V', '-k2,2n', input_bed_path])

    with open(bed_path, '+w') as f:
        f.write(sorted_bed.decode('utf-8'))

    subprocess.check_call(['bgzip', bed_path])
    subprocess.check_call(['tabix', '-p', 'bed', bed_path + '.gz'])

    with open(bed_path, '+w') as f:  # write uncompressed file again
        f.write(sorted_bed.decode('utf-8'))

    return api.patch_instance(
        endpoint='techniques',
        identifier=instance['pk'],
        storage_url=data_dir,
        storage_usage=utils.get_tree_size(data_dir),
        bed_url=bed_path)


class LocalDataImporter():

    """
    A Data import engine for workflows.

    Attributes:
        FASTQ_REGEX (str): a regex pattern used to match fastq files.
        BAM_REGEX (str): a regex pattern to match bams.
    """

    BAM_REGEX = r'\.bam(\.[a-z0-9]+)?$'
    FASTQ_REGEX = r'(([_.]R{0}[_.].+)|([_.]R{0}\.)|(_{0}\.))f(ast)?q(\.gz)?$'

    def __init__(self):
        """Initialize cache to None."""
        self.cache = None

    def import_data(
            self, directories, symlink=False, commit=False,
            key=lambda x: x['system_id'], **filters):
        """
        Import raw data for multiple workflows.

        Workflows's `storage_url`, `storage_usage`, `data_type` are updated,
        setting the latter to the data type found (e.g. FASTQ, BAM).

        Arguments:
            directories (list): list of directories to be recursively explored.
            symlink (bool): if True symlink instead of moving.
            commit (bool): if True perform import operation.
            key (function): given a workflow dict returns id to match.
            filters (dict): key value pairs to use as API query params.

        Raises:
            click.UsageError: if `key` returns the same identifier for multiple
                workflows. If a workflow matches both fastq and bam files.
                if cant determine read 1 or read 2 from matched fastq files.

        Returns:
            list: of workflows dicts for which data has been matched.
        """
        utils.check_admin()
        workflows_matched, files_matched = [], 0
        data_storage_dir = system_settings.BASE_STORAGE_DIRECTORY
        pattern = self._build_cache(key, filters)

        if pattern:
            label = f'Exploring directories...'
            for directory in directories:
                with click.progressbar(os.walk(directory), label=label) as bar:
                    for root, _, files in bar:
                        if root.startswith(data_storage_dir):
                            continue

                        for i in files:
                            path = join(root, i)
                            matched = self._update_src_dst(path, pattern)
                            files_matched += 1 if matched else 0

        if commit and files_matched:
            label = f'Processing {files_matched} matched files...'
            with click.progressbar(self.cache.values(), label=label) as bar:
                for i in sorted(bar, key=lambda x: x['workflow']['pk']):
                    if i['src_dst_tuples']:
                        for src, dst in i['src_dst_tuples']:
                            if symlink:
                                self.symlink(src, dst)
                            else:
                                self.move(src, dst)

                        workflows_matched.append(api.patch_instance(
                            endpoint='workflows',
                            identifier=i['workflow']['pk'],
                            storage_url=i['storage_url'],
                            storage_usage=utils.get_tree_size(i['storage_url']),
                            data_type=i['data_type']))

        elif files_matched:  # pragma: no cover
            for i in self.cache.values():
                if i['src_dst_tuples']:
                    workflows_matched.append(i['workflow'])

        click.echo(self.get_summary())
        return workflows_matched

    def get_summary(self):
        """Get a summary of the matched, skipped, and missing files."""
        skipped, missing, matched, total_matched, nl = [], [], [], 0, '\n'

        for i in self.cache.values():
            if i['workflow']['data_type']:
                msg = click.style(f"skipped {i['uid']}\t", fg='cyan')
                msg += f"{i['workflow']['data_type']} exists"
                skipped.append(msg)
            elif i['src_dst_tuples']:
                msg = click.style(f"found {i['uid']}\n\t\t", fg='green')
                msg += '\n\t\t'.join(src for src, _ in i['src_dst_tuples'])
                total_matched += len(i['src_dst_tuples'])
                matched.append(msg)
            else:
                msg = click.style(f"missing {i['uid']}\t", fg='red')
                msg += 'no files matched'
                missing.append(msg)

        return (
            f"{nl.join([nl] + skipped) if skipped else ''}"
            f"{nl.join([nl] + missing) if missing else ''}"
            f"{nl.join([nl] + matched) if matched else ''}"
            f'\n\ntotal samples: {len(self.cache)}'
            f'\nsamples skipped: {len(skipped)}'
            f'\nsamples missing: {len(missing)}'
            f'\nsamples matched: {len(matched)}'
            f'\ntotal files matched: {total_matched}')

    @classmethod
    def as_cli_command(cls):
        """Get data importerls as click command line interface."""
        @click.command(help='local data import', name='import_data')
        @options.DIRECTORIES
        @options.IDENTIFIER
        @options.FILTERS
        @options.COMMIT
        @options.SYMLINK
        def command(identifier, commit, filters, directories, symlink):
            """Click command to be used in the CLI."""
            def key(workflow):
                value, types = workflow, (int, str, type(None))
                for i in identifier:
                    value = value.get(i)
                if not isinstance(value, types):
                    raise click.UsageError(
                        f'invalid type for identifier '
                        f'`{".".join(identifier)}`: {type(value)}')
                return value

            matched = cls().import_data(
                directories, symlink, commit, key, **filters)

            if not commit and matched:  # pragma: no cover
                utils.echo_add_commit_message()

        return command

    @staticmethod
    def symlink(src, dst):
        """Create symlink from `src` to `dst`."""
        return utils.force_symlink(os.path.realpath(src), dst)

    @staticmethod
    def move(src, dst):
        """Rename `src` to `dst`."""
        return os.rename(os.path.realpath(src), dst)

    @staticmethod
    def get_regex_pattern(group_name, identifier):
        """
        Get regex pattern for `identifier` group as `group_name`.

        This pattern treats dashes, underscores and dots equally.

        Arguments:
            group_name (str): regex pattern group name.
            identifier (str): identifier to be matched by regex.

        Returns:
            str: a regex pattern.
        """
        pattern = re.sub(r'[-_.]', r'[-_.]', identifier)
        return r'(?P<{}>[-_.]?{}[-_.])'.format(group_name, pattern)

    def _build_cache(self, key, filters):
        """
        Build cache dictionary and return a regex pattern to match identifiers.

        Workflows for which `data_type` is set will not be included in the
        regex pattern.

        Set self.cache to a dict that looks like:

            {
                "primary_key_1": {
                    "workflow": workflow dict as returned from API,
                    "src_dst_tuples": list of src, dst tuples to store matches,
                    "dtype": key used to store data type matched,
                    "data_dir": data dir assigned to workflow,
                    },

                "primary_key_2": ...
            }

        Arguments:
            key (function): given a workflow dict returns id to match.
            filters (dict): key value pairs to use as API query params.

        Returns:
            str: compiled regex pattern to match identifiers to workflows.
        """
        self.cache = {}
        patterns = []
        identifiers = {}
        get_dir = system_settings.GET_STORAGE_DIRECTORY

        for i in api.get_instances('workflows', verbose=True, **filters):
            index = f"primary_key_{i['pk']}"
            identifier = key(i) or 'IDENTIFIER IS NULL'
            storage_url = i['storage_url'] or get_dir('workflows', i['pk'])
            data_dir = join(storage_url, 'data')

            if not i['storage_url']:
                os.makedirs(data_dir, exist_ok=True)

            if identifier in identifiers:  # duplicated identifiers not valid
                raise click.UsageError(
                    f"Can't use same identifier for {i['system_id']} "
                    f'and {identifiers[identifier]}: {identifier}')
            elif identifier != 'IDENTIFIER IS NULL':
                identifiers[identifier] = i['system_id']

            self.cache[index] = {}
            self.cache[index]['workflow'] = i
            self.cache[index]['identifier'] = identifier
            self.cache[index]['src_dst_tuples'] = []
            self.cache[index]['data_type'] = None
            self.cache[index]['storage_url'] = storage_url
            self.cache[index]['data_dir'] = data_dir
            self.cache[index]['uid'] = f"{i['system_id']} (using {identifier})"

            if not i['data_type']:
                patterns.append(self.get_regex_pattern(index, identifier))

        # see http://stackoverflow.com/questions/8888567
        return re.compile('|'.join(patterns))

    def _update_src_dst(self, path, pattern):
        """Match `path` with `pattern` and update cache if fastq or bam."""
        try:
            matches = pattern.finditer(path)
            index = next(matches).lastgroup
            assert index is not None  # happens when pattern is empty
            uid = self.cache[index]['uid']
            data_dir = self.cache[index]['data_dir']
            system_id = self.cache[index]['workflow']['system_id']
            src_dst_tuples = self.cache[index]['src_dst_tuples']
        except (StopIteration, AssertionError):  # pragma: no cover
            return None

        dst = None
        bam_dst = path if re.search(self.BAM_REGEX, path) else None
        fastq_dst = self._format_fastq_path(path)

        if bam_dst or fastq_dst:
            dst = basename(bam_dst or fastq_dst)
            dst = dst if dst.startswith(system_id) else f'{system_id}_{dst}'
            data_type = 'BAM' if bam_dst else 'FASTQ'
            src_dst_tuples.append((path, join(data_dir, dst)))

            if self.cache[index]['data_type'] in {None, data_type}:
                self.cache[index]['data_type'] = data_type
            else:
                raise click.UsageError(
                    f'{uid} matched different data types: '
                    f"{', '.join(i for i, _ in src_dst_tuples)}")

        return dst is not None

    def _format_fastq_path(self, path):
        """Format fastq file name if `path` is valid fastq path."""
        for i in [1, 2]:
            letter_index_fastq = r'[_.]R{}([_.])?\.f(ast)?q'.format(i)
            number_index_fastq = r'[_.]{}([_.])?\.f(ast)?q'.format(i)
            letter_index_anyix = r'[_.]R{}[_.]'.format(i)

            if re.search(self.FASTQ_REGEX.format(i), path):
                suffix = f'_{system_settings.FASTQ_READ_PREFIX}{i}.fastq'
                dst = re.sub(letter_index_fastq, '.fastq', path)
                dst = re.sub(number_index_fastq, '.fastq', dst)
                dst = re.sub(letter_index_anyix, '_', dst)
                return re.sub(r'[_.]f(ast)?q', suffix, dst)

        if re.search(r'\.f(ast)?q(\.gz)?$', path):
            msg = f'cant determine if read 1 or read 2 from: {path}'
            raise click.UsageError(msg)

        return None
