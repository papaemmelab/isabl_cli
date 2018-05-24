"""Data import logic."""

from os.path import basename
from os.path import join
from os.path import dirname
import collections
import os
import re

import click

from cli import api
from cli import system_settings
from cli import utils


def import_bedfile():  # pragma: no cover
    """Register bedfile in technique data directory."""
    # we need to make sure bedfiles don't have the chr prefix
    bedfiles = []

    for i in bedfiles:
        with open(i, 'r') as f:
            if f.readline().startswith('chr'):
                msg = "Bedfiles can't have 'chr' prefixes: " + i
                raise click.UsageError(msg)


def get_data_dir(endpoint, primary_key):
    base = system_settings.DATA_STORAGE_DIRECTORY
    hash_1 = f'{primary_key:04d}' [-4:-2]
    hash_2 = f'{primary_key:04d}' [-2:]
    use_hash = {'individuals', 'spcimens', 'workflows', 'analyses'}

    if not base:  # pragma: no cover
        return 'Setting `DATA_STORAGE_DIRECTORY` not defined.'

    if endpoint in use_hash and base:
        path = os.path.join(endpoint, hash_1, hash_2)
    else:
        path = os.path.join(endpoint)

    return os.path.join(base, 'workflow', path, str(primary_key))


class LocalDataImporter():

    FASTQ_REGEX = r'(([_.]R{0}[_.].+)|([_.]R{0}\.)|(_{0}\.))f(ast)?q(\.gz)?$'

    def __init__(self):
        self.cache = None

    def import_data(
            self, directories, symlink=False, commit=False,
            key=lambda x: x['system_id'], **filters):
        utils.check_admin()
        patterns, imported, files_matched = [], [], 0
        data_dir_fn = system_settings.GET_DATA_DIR_FUNCTION
        self.cache = {}

        for i in api.get_instances('workflows', verbose=True, **filters):
            index = f"primary_key_{i['pk']}"
            identifier = str(key(i))
            self.cache[index] = {}
            self.cache[index]['workflow'] = i
            self.cache[index]['identifier'] = identifier
            self.cache[index]['src_dst_tuples'] = []
            self.cache[index]['dtype'] = None
            self.cache[index]['data_dir'] = data_dir_fn('workflows', i['pk'])
            self.cache[index]['uid'] = f"{i['system_id']} (using {identifier})"

            if not i['data_type']:
                patterns.append(self._get_regex_pattern(index, identifier))

        if patterns:
            label = f'Exploring directories...'
            pattern = re.compile('|'.join(patterns))

            for directory in directories:
                with click.progressbar(os.walk(directory), label=label) as bar:
                    # see http://stackoverflow.com/questions/8888567
                    for root, _, files in bar:
                        for i in files:
                            path = join(root, i)
                            matched = self._update_src_dst(path, pattern)
                            files_matched += 1 if matched else 0

        if commit and files_matched:
            label = f'Processing {files_matched} matched files...'
            with click.progressbar(self.cache.values(), label=label) as bar:
                for i in bar:
                    if i['src_dst_tuples']:
                        os.makedirs(i['data_dir'], exist_ok=True)

                        for src, dst in i['src_dst_tuples']:
                            if symlink:
                                self._symlink(src, dst)
                            else:
                                self._move(src, dst)

                        imported.append(api.patch_instance(
                            endpoint='workflows',
                            identifier=i['workflow']['pk'],
                            storage_url=i['data_dir'],
                            storage_usage=utils.get_tree_size(i['data_dir']),
                            data_type=i['dtype']))

        summary = self._get_summary()
        click.echo(summary)

        return imported

    @staticmethod
    def _symlink(src, dst):
        return utils.force_symlink(os.path.realpath(src), dst)

    @staticmethod
    def _move(src, dst):
        return os.rename(os.path.realpath(src), dst)

    @staticmethod
    def _get_regex_pattern(group_name, identifier):
        """
        Get regex pattern for `identifier` that sets `group_name` as group.

        This pattern treats dashes, underscores and dots equally.

        Arguments:
            group_name (str): regex pattern group name.
            identifier (str): identifier to be matched by regex.

        Returns:
            str: a regex pattern.
        """
        pattern = re.sub(r'[-_.]', r'[-_.]', identifier)
        return r'(?P<{}>[-_.]?{}[-_.])'.format(group_name, pattern)

    def _update_src_dst(self, path, pattern):
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
        bam_dst = self._get_bam_dst(path)
        fastq_dst = self._get_fastq_dst(path)

        if bam_dst or fastq_dst:
            try:  # raise error if multiple samples match the same file.
                uid2 = self.cache[next(matches).lastgroup]['uid']
                raise click.UsageError(f'{uid} and {uid2} matched: {path}')
            except StopIteration:  # pragma: no cover
                pass

            dst = bam_dst or fastq_dst
            dst = join(data_dir, system_id + '_' + basename(dst))
            dtype = 'BAM' if bam_dst else 'FASTQ'
            src_dst_tuples.append((path, dst))

            if self.cache[index]['dtype'] in {None, dtype}:
                self.cache[index]['dtype'] = dtype
            else:
                raise click.UsageError(
                    f'{uid} matched different data types: '
                    f"{', '.join(i for i, _ in src_dst_tuples)}")


        return dst is not None

    def _get_bam_dst(self, path):
        return path if re.search(r'\.bam(\.[a-z0-9]+)?$', path) else None

    def _get_fastq_dst(self, path):
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

    def _get_summary(self):
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
