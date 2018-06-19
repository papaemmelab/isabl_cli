"""Data import logic."""

from collections import defaultdict
from datetime import datetime
from getpass import getuser
from os.path import basename
from os.path import getsize
from os.path import isdir
from os.path import join
import os
import re
import shutil
import subprocess

import click
import yaml

from cli import api
from cli import options
from cli import system_settings
from cli import utils


def symlink_workflow_to_projects(workflow):
    """Create symlink from workflow directory and projects directories."""
    for i in workflow['projects']:
        if not i['storage_url']:
            i = make_storage_diectory('projects', i['pk'])

        utils.force_symlink(
            workflow['storage_url'],
            join(i['storage_url'], workflow['system_id']))


def symlink_analysis_to_targets(analysis):
    """Create symlink from workflow directory and projects directories."""
    src = analysis['storage_url']
    dst = '__'.join([
        analysis['pipeline']['name'].lower().replace(' ', '_'),
        analysis['pipeline']['version'].lower().replace(' ', '_'),
        str(analysis['pk'])])

    for i in analysis['targets']:
        if not i['storage_url']:
            i = make_storage_diectory('workflows', i['pk'])
        utils.force_symlink(src, join(i['storage_url'], dst))

    if analysis['project_level_analysis']:
        i = analysis['project_level_analysis']
        if not i['storage_url']:
            i = make_storage_diectory('projects', i['pk'])
        utils.force_symlink(src, join(i['storage_url'], dst))


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


def make_storage_diectory(endpoint, primary_key):
    """Make storage directory and return patched instance."""
    get_dir = system_settings.GET_STORAGE_DIRECTORY
    storage_url = get_dir(endpoint, primary_key, use_hash=False)
    os.makedirs(storage_url, exist_ok=True)
    return api.patch_instance(endpoint, primary_key, storage_url=storage_url)


class LocalBedImporter():

    """An import engine for techniques' bedfiles."""

    @classmethod
    def as_cli_command(cls):
        """Get bed importer as click command line interface."""
        @click.command(name='import_bed')
        @options.PRIMARY_KEY
        @options.BEDFILE
        @click.option('--assembly', help='name of reference genome')
        def command(key, bedfile, assembly):
            """
            Register a `bedfile` in technique`s data directory.

            The incoming bedfile will be compressed and tabixed.
            Both gzipped and uncompressed versions are kept.

            Instance's `storage_url`, `storage_usage`, `data` are updated,
            setting the latter's `bedfiles` key to:

                'bedfiles': {
                    <assembly>: {
                        'uncompressed': path/to/bedfile.bed,
                        'gzipped': path/to/bedfile.bed.gz
                    }
                }
            """
            cls().import_bed(key, bedfile, assembly)
        return command

    @staticmethod
    def import_bed(technique_primary_key, input_bed_path, assembly):
        """
        Register input_bed_path in technique's storage dir and update `data`.

        Arguments:
            technique_primary_key (int): technique primary key.
            input_bed_path (str): path to incoming bedfile.
            assembly (str): name of reference genome for bedfile.

        Returns:
            dict: updated technique instance as retrieved from API.
        """
        utils.check_admin()
        instance = api.get_instance('techniques', technique_primary_key)

        if not input_bed_path.endswith('.bed'):  # pragma: no cover
            raise click.UsageError(f'No .bed suffix: {input_bed_path}')

        if not instance['storage_url']:
            instance = make_storage_diectory('techniques', instance['pk'])

        if instance['data'].get('bedfiles', {}).get(assembly):
            raise click.UsageError(
                f'{instance["slug"]} has a bed registered for {assembly}: '
                f'{instance["data"]["bedfiles"][assembly]}')

        if not instance['data'].get('bedfiles'):
            instance['data']['bedfiles'] = {}

        bed_dir = join(instance['storage_url'], assembly)
        bed_path = join(bed_dir, f"{instance['slug']}.bed")
        os.makedirs(bed_dir, exist_ok=True)
        sorted_bed = subprocess.check_output(
            ['sort', '-k1,1V', '-k2,2n', input_bed_path])

        with open(bed_path, '+w') as f:
            f.write(sorted_bed.decode('utf-8'))

        subprocess.check_call(['bgzip', bed_path])
        subprocess.check_call(['tabix', '-p', 'bed', bed_path + '.gz'])

        with open(bed_path, '+w') as f:  # write uncompressed file again
            f.write(sorted_bed.decode('utf-8'))

        instance['data']['bedfiles'][assembly] = {}
        instance['data']['bedfiles'][assembly]['uncompressed'] = bed_path
        instance['data']['bedfiles'][assembly]['gzipped'] = bed_path + '.gz'

        return api.patch_instance(
            endpoint='techniques',
            identifier=instance['pk'],
            storage_usage=utils.get_tree_size(instance['storage_url']),
            data=instance['data'])


class LocalDataImporter():

    """
    A Data import engine for workflows.

    Attributes:
        FASTQ_REGEX (str): a regex pattern used to match fastq files.
        BAM_REGEX (str): a regex pattern to match bams.
    """

    BAM_REGEX = r'\.bam$'
    CRAM_REGEX = r'\.cram$'
    FASTQ_REGEX = r'(([_.]R{0}[_.].+)|([_.]R{0}\.)|(_{0}\.))f(ast)?q(\.gz)?$'

    def import_data(
            self, directories, symlink=False, commit=False,
            key=lambda x: x['system_id'], files_data={}, **filters):
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
            files_data (dict): keys are files basenames and values are
                dicts with extra annotations such as PL, LB, or any other.

        Raises:
            click.UsageError: if `key` returns the same identifier for multiple
                workflows. If a workflow matches both fastq and bam files.
                if cant determine read 1 or read 2 from matched fastq files.

        Returns:
            list, str: list of workflows for which data has been matched and a
                summary of the operation.
        """
        utils.check_admin()
        workflows_matched = []
        cache = defaultdict(dict)
        patterns = []
        identifiers = {}

        for i in api.get_instances('workflows', verbose=True, **filters):
            index = f"primary_key_{i['pk']}"
            using_id = f"{i['system_id']} (Skipped, identifier is NULL)"
            identifier = key(i)

            if identifier in identifiers:  # duplicated identifiers not valid
                raise click.UsageError(
                    f"Can't use same identifier for {i['system_id']} "
                    f'and {identifiers[identifier]}: {identifier}')

            if identifier and not i['sequencing_data']:
                identifiers[identifier] = i['system_id']
                patterns.append(self.get_regex_pattern(index, identifier))
                using_id = f"{i['system_id']} (using {identifier})"

            cache[index]['using_id'] = using_id
            cache[index]['instance'] = i
            cache[index]['files'] = []

        if patterns:
            # see http://stackoverflow.com/questions/8888567
            data_storage_dir = system_settings.BASE_STORAGE_DIRECTORY
            pattern = re.compile('|'.join(patterns))
            label = f'Exploring directories...'
            for directory in directories:
                with click.progressbar(os.walk(directory), label=label) as bar:
                    for root, _, files in bar:
                        if not root.startswith(data_storage_dir):
                            for i in files:
                                path = join(root, i)
                                index = self.match_path(path, pattern)
                                if index:
                                    cache[index]['files'].append(path)

            label = 'Processing...'
            bar = sorted(cache.values(), key=lambda x: x['instance']['pk'])
            with click.progressbar(bar, label=label) as bar:
                for i in bar:
                    if commit and i['files']:
                        workflows_matched.append(self.import_files(
                            instance=i['instance'],
                            files=i['files'],
                            symlink=symlink,
                            files_data=files_data))
                    elif i['files']:  # pragma: no cover
                        workflows_matched.append(i['instance'])

        return workflows_matched, self.get_summary(cache)

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

    def match_path(self, path, pattern):
        """Match `path` with `pattern` and update cache if fastq or bam."""
        try:
            matches = pattern.finditer(path)
            index = next(matches).lastgroup
            assert index is not None  # happens when pattern is empty

            # check if valid data type
            valid = True if re.search(self.BAM_REGEX, path) else False
            valid |= True if re.search(self.CRAM_REGEX, path) else False

            for i in [1, 2]:
                if re.search(self.FASTQ_REGEX.format(i), path):
                    valid = True

            if re.search(r'\.f(ast)?q(\.gz)?$', path) and not valid:
                msg = f'cant determine if read 1 or read 2 from: {path}'
                raise click.UsageError(msg)

            assert valid
        except (StopIteration, AssertionError):  # pragma: no cover
            return None

        return index

    def import_files(self, instance, files, files_data, symlink):
        """
        Move/link files into instance's `storage_url` and update database.

        Arguments:
            instance (dict): workflow instance.
            files (dict): list of files to be imported.
            symlink (dict): whether to symlink or move the data.
            files_data (dict): keys are files basenames and values are
                dicts with extra annotations such as PL, LB, or any other.

        Raises:
            click.UsageError: if multiple data formats are found.

        Returns:
            dict: patched workflow instance.
        """
        sequencing_data = []
        src_dst = []
        get_dir = system_settings.GET_STORAGE_DIRECTORY

        if not instance['storage_url']:
            instance['storage_url'] = get_dir('workflows', instance['pk'])

        data_dir = join(instance['storage_url'], 'data')
        os.makedirs(data_dir, exist_ok=True)

        for src in files:
            file_name = basename(src)

            if re.search(self.BAM_REGEX, src):
                file_type = 'BAM'
            elif re.search(self.CRAM_REGEX, src):
                file_type = 'CRAM'
            else:
                file_type = 'FASTQ'
                file_name = self.format_fastq_name(file_name)

            if not file_name.startswith(instance['system_id']):
                file_name = f'{instance["system_id"]}__{file_name}'

            dst = join(data_dir, file_name)
            src_dst.append((src, dst))
            sequencing_data.append(dict(
                file_url=dst,
                file_type=file_type,
                file_data=files_data.get(basename(src), {}),
                hash_value=getsize(src),
                hash_method="os.path.getsize"))

        if len({i['file_type'] for i in sequencing_data}) > 1:
            raise click.UsageError(
                'We should have catched this earlier, but multiple formats are '
                f'not supported, these were found for {instance["system_id"]}: '
                f'{",".join(i["file_url"] for i in sequencing_data)}')

        for src, dst in src_dst:
            if symlink:
                self.symlink(src, dst)
            else:
                self.move(src, dst)

        return api.patch_instance(
            endpoint='workflows',
            identifier=instance['pk'],
            storage_url=instance['storage_url'],
            storage_usage=utils.get_tree_size(instance['storage_url']),
            sequencing_data=sequencing_data)

    def format_fastq_name(self, file_name):
        """Return destination file name."""
        suffix = None

        for i in [1, 2]:
            if re.search(self.FASTQ_REGEX.format(i), file_name):
                suffix = f'_{system_settings.FASTQ_READ_PREFIX}{i}.fastq'
                break

        assert suffix, f"Couldn't determine read 1 or read 2 from {file_name}"
        letter_index_fastq = r'[_.]R{}([_.])?\.f(ast)?q'.format(i)
        number_index_fastq = r'[_.]{}([_.])?\.f(ast)?q'.format(i)
        letter_index_any_location = r'[_.]R{}[_.]'.format(i)
        file_name = re.sub(letter_index_fastq, '.fastq', file_name)
        file_name = re.sub(number_index_fastq, '.fastq', file_name)
        file_name = re.sub(letter_index_any_location, '_', file_name)
        return re.sub(r'[_.]f(ast)?q', suffix, file_name)

    @staticmethod
    def symlink(src, dst):
        """Create symlink from `src` to `dst`."""
        return utils.force_symlink(os.path.realpath(src), dst)

    @staticmethod
    def move(src, dst):
        """Rename `src` to `dst`."""
        return os.rename(os.path.realpath(src), dst)

    @staticmethod
    def get_summary(cache):
        """Get a summary of the matched, skipped, and missing files."""
        skipped, missing, matched, total_matched, nl = [], [], [], 0, '\n'

        for i in cache.values():
            if i['instance']['sequencing_data']:
                msg = click.style(f"skipped {i['using_id']}\t", fg='cyan')
                skipped.append(msg)
            elif i['files']:
                msg = click.style(f"found {i['using_id']}\n\t\t", fg='green')
                total_matched += len(i['files'])
                matched.append(msg + '\n\t\t'.join(i['files']))
            else:
                msg = click.style(f"missing {i['using_id']}\t", fg='red')
                missing.append(msg + 'no files matched')

        return (
            f"{nl.join([nl] + skipped) if skipped else ''}"
            f"{nl.join([nl] + missing) if missing else ''}"
            f"{nl.join([nl] + matched) if matched else ''}"
            f'\n\ntotal samples: {len(cache)}'
            f'\nsamples skipped: {len(skipped)}'
            f'\nsamples missing: {len(missing)}'
            f'\nsamples matched: {len(matched)}'
            f'\ntotal files matched: {total_matched}')

    @classmethod
    def as_cli_command(cls):
        """Get data importer as a click command line interface."""
        @click.command(name='import_data')
        @options.DIRECTORIES
        @options.IDENTIFIER
        @options.FILTERS
        @options.COMMIT
        @options.SYMLINK
        @options.FILES_DATA
        def command(
                identifier, commit, filters, directories, symlink, files_data):
            """
            Find and import data for multiple workflows from many directories.

            Search is recursive and any cram, bam or fastq file that matches
            the workflow identifier will be imported.

            Its possible to provide custom annotation per file (e.g. PL, PU, or
            LB in the case of fastq data). In order to do so, provide a yaml
            file using the `--files-data` argument. Such file must look
            like this, please note that keys are file names, not full paths:

            \b
                1.fq:
                    ID: 9
                    LB: Library_id
                    PL: ILLUMINA
                    PM: HiSeq-XTen
                    PU: MICHELLE
                2.fq:
                    LB: Library2
                    ...
            """
            def key(workflow):
                value, types = workflow, (int, str, type(None))
                for i in identifier:
                    value = value.get(i)
                if not isinstance(value, types):
                    raise click.UsageError(
                        f'invalid type for identifier '
                        f'`{".".join(identifier)}`: {type(value)}')
                return value

            if files_data:
                with open(files_data) as f:
                    files_data = yaml.load(f.read())
            else:
                files_data = {}

            matched, summary = cls().import_data(
                directories=directories, symlink=symlink, commit=commit,
                key=key, files_data=files_data, **filters)

            click.echo(summary)

            if not commit and matched:  # pragma: no cover
                utils.echo_add_commit_message()

        return command
