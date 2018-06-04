"""Data import logic."""

from datetime import datetime
from os.path import basename
from os.path import join
from os.path import isdir
import shutil
import os
import re
import subprocess
from getpass import getuser

import click

from cli import api
from cli import system_settings
from cli import utils


def trash_analysis_storage(analysis):
    """Move analysis `storage_url` to a trash directory."""
    if analysis['status'] == 'SUCCEEDED':
        raise click.UsageError("You can't wipe a succeeded analysis")

    slug = f'primary_key_{analysis["pk"]}__user_{getuser()}__date_'
    slug += datetime.now(system_settings.TIME_ZONE).isoformat()

    if isdir(analysis['storage_url']):
        assert 'analyses' in analysis['storage_url']
        trash_directory = get_storage_directory(
            endpoint='.trashed_analyses',
            primary_key=analysis['pk'],
            base_directory=analysis['storage_url'].split('analyses')[0])

        os.makedirs(trash_directory, exist_ok=True)
        dst = join(trash_directory, slug)
        click.echo(f"\n\ntrashing: {analysis['storage_url']} -> {dst}")
        shutil.move(analysis['storage_url'], dst)


def get_storage_directory(endpoint, primary_key, base_directory=None):
    """
    Get path to instance's data directory.

    A hash naming system using the primay key is used for individuals,
    spcimens, workflows and analyses. For example given analyses with
    primary key 12345 and 2345:

        {base_directory}/analyses/23/45/2345
        {base_directory}/analyses/23/45/12345

    For other models such as projects or techniques the primary key is used
    naively:

        {base_directory}/projects/100
        {base_directory}/techniques/2

    Arguments:
        endpoint (str): instance's API endpoint.
        primary_key (str): instance's primary key.
        base_directory (str): default is system_settings.BASE_STORAGE_DIRECTORY.

    Returns:
        str: path to instance's data directory.
    """
    base_directory = base_directory or system_settings.BASE_STORAGE_DIRECTORY
    hash_1 = f'{primary_key:04d}' [-4:-2]
    hash_2 = f'{primary_key:04d}' [-2:]
    dont_use_hash = {'projects', 'pipelines', 'techniques'}

    if not base_directory:  # pragma: no cover
        raise click.UsageError('Setting `BASE_STORAGE_DIRECTORY` not defined.')

    if endpoint in dont_use_hash:
        path = os.path.join(endpoint)
    else:
        path = os.path.join(endpoint, hash_1, hash_2)

    return os.path.join(base_directory, path, str(primary_key))


def import_bedfile(technique_primary_key, input_bed_path):
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
    data_dir_fn = system_settings.GET_STORAGE_DIRECTORY_FUNCTION
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
    """

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
        imported, files_matched = [], 0
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
                        os.makedirs(i['data_dir'], exist_ok=True)

                        for src, dst in i['src_dst_tuples']:
                            if symlink:
                                self.symlink(src, dst)
                            else:
                                self.move(src, dst)

                        imported.append(api.patch_instance(
                            endpoint='workflows',
                            identifier=i['workflow']['pk'],
                            storage_url=i['data_dir'],
                            storage_usage=utils.get_tree_size(i['data_dir']),
                            data_type=i['data_type']))

        summary = self.get_summary()
        click.echo(summary)

        return imported

    def format_bam_dst(self, path):
        """Return `path` if its a valid bam file name."""
        return path if re.search(r'\.bam(\.[a-z0-9]+)?$', path) else None

    def format_fastq_dst(self, path):
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
        data_dir_fn = system_settings.GET_STORAGE_DIRECTORY_FUNCTION

        for i in api.get_instances('workflows', verbose=True, **filters):
            index = f"primary_key_{i['pk']}"
            identifier = key(i)

            if identifier in identifiers:  # duplicated identifiers not valid
                raise click.UsageError(
                    f"{key} returned the same identifier for {i['system_id']} "
                    f'and {identifiers[identifier]}: {identifier}')
            elif not identifier:
                identifier = 'IDENTIFIER IS NULL'
            else:
                identifiers[identifier] = i['system_id']

            self.cache[index] = {}
            self.cache[index]['workflow'] = i
            self.cache[index]['identifier'] = identifier
            self.cache[index]['src_dst_tuples'] = []
            self.cache[index]['data_type'] = None
            self.cache[index]['data_dir'] = data_dir_fn('workflows', i['pk'])
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
            uid = self.cache[index]['uid']
            data_dir = self.cache[index]['data_dir']
            system_id = self.cache[index]['workflow']['system_id']
            src_dst_tuples = self.cache[index]['src_dst_tuples']
        except StopIteration:  # pragma: no cover
            return

        dst = None
        bam_dst = self.format_bam_dst(path)
        fastq_dst = self.format_fastq_dst(path)

        if bam_dst or fastq_dst:
            dst = bam_dst or fastq_dst
            dst = join(data_dir, system_id + '_' + basename(dst))
            data_type = 'BAM' if bam_dst else 'FASTQ'
            src_dst_tuples.append((path, dst))

            if self.cache[index]['data_type'] in {None, data_type}:
                self.cache[index]['data_type'] = data_type
            else:
                raise click.UsageError(
                    f'{uid} matched different data types: '
                    f"{', '.join(i for i, _ in src_dst_tuples)}")

        return dst is not None
