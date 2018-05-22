from os.path import basename
from os.path import join
import collections
import os
import re

import click

from cli import api
from cli import system_settings
from cli import utils


def import_bedfile():
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

    if not base:
        return 'Setting `DATA_STORAGE_DIRECTORY` not defined.'

    if endpoint in use_hash and base:
        path = os.path.join(endpoint, hash_1, hash_2)
    else:
        path = os.path.join(endpoint)

    return os.path.join(base, 'workflow', path, str(primary_key))


class LocalDataImporter():

    def __init__(self):
        self.cache = None

    def import_raw_data(
            self, filters, directories,
            symlink=False, commit=False, key=lambda x: x['system_id']):
        utils._check_admin()
        patterns, imported, files_matched = [], [], 0
        self.cache = collections.defaultdict(dict)

        for i in api.get_instances(verbose=True, **filters):
            pk = i['pk']
            identifier = str(key(i))
            data_dir = system_settings.GET_DATA_DIR_FUNCTION('workflows', pk)
            self.cache[pk]['workflow']: i
            self.cache[pk]['identifier']: identifier
            self.cache[pk]['src_dst_tuples']: []
            self.cache[pk]['dtype']: None
            self.cache[pk]['data_dir']: data_dir
            self.cache[pk]['uid']: f"{i['system_id']} (using {identifier})"

            if not i['raw_data']:
                patterns.append(self.get_regex_pattern(i['pk'], identifier))

        label = f'Exploring directories...'
        for directory in directories:
            with click.progressbar(os.walk(directory), label=label) as bar:
                # see http://stackoverflow.com/questions/8888567
                pattern = re.compile('|'.join(patterns))
                for root, _, files in bar:
                    for i in files:
                        path = join(root, i)
                        matched = self._update_src_dst(path, pattern)
                        files_matched += 1 if matched else 0

        if commit:
            label = f'Processing {files_matched} matched files...'
            with click.progressbar(self.cache.items(), label=label) as bar:
                for pk, i in bar:
                    for src, dst in i['src_dst_tuples']:
                        if symlink:
                            self.symlink(src, dst)
                        else:
                            self.move(src, dst)

                    imported.append(api.patch_instance(
                        endpoint='workflows',
                        identifier=pk,
                        data_url=i['data_dir'],
                        data_usage=utils.get_tree_size(i['data_dir']),
                        data_type=i['workflow']['dtype']))

        self._echo_summary()
        return imported

    @staticmethod
    def symlink(src, dst):
        return utils.force_symlink(os.path.realpath(src), dst)

    @staticmethod
    def move(src, dst):
        return os.rename(os.path.realpath(src), dst)

    def is_valid_fastq(self, path):
        return path.endswith('.fastq.gz')

    def is_valid_bam(self, path):
        return path.endswith('.bam') or '.bam.' in path

    @staticmethod
    def get_regex_pattern(primary_key, identifier):
        """
        Get regex pattern for `identifier` that uses `primary_key` as group.

        This pattern treats dashes, underscores and dots equally.

        Arguments:
            primary_key (int): primary key used as group in regex pattern.
            identifier (str): identifier to be matched by regex.

        Returns:
            str: a regex pattern.
        """
        pattern = re.sub(r'[-_.]', r'[-_.]', identifier)
        return r'(?P<{}>[-_.]?{}[-_.])'.format(primary_key, pattern)

    def _update_src_dst(self, path, pattern):
        dst = None
        matches = pattern.finditer(path)

        try:
            pk = next(matches).lastgroup
            uid = self.cache[pk]['uid']
            src_dst_tuples = self.cache[pk]['src_dst_tuples']
        except StopIteration:
            return

        try:
            # raise error if multiple samples match the same file.
            uid2 = self.cache[next(matches).lastgroup]['uid']
            raise click.UsageError(f'{uid} and {uid2} matched: {path}')
        except StopIteration:
            pass

        if self.is_valid_bam(path):
            src_dst_tuples.append((path, self._get_bam_dst(pk, path)))
            dtype = 'BAM'
        elif self.is_valid_fastq(path):
            src_dst_tuples.append((path, self._get_fastq_dst(pk, path)))
            dtype = 'FASTQ'

        if self.cache['pk']['dtype'] in {None, dtype}:
            self.cache['pk']['dtype'] = dtype
        else:
            raise click.UsageError(
                f'{uid} matched different data types: '
                f"{', '.join(i for i, _ in src_dst_tuples)}")

        return dst is not None

    def _get_bam_dst(self, pk, path):
        return join(
            self.cache[pk]['data_dir'],
            self.cache[pk]['workflow']['system_id'] + '_' + basename(path))

    def _get_fastq_dst(self, pk, path):
        if re.search('_[1|2].fastq', path):
            dst = path
        elif re.search('[_.]R1[_.]', path):
            dst = re.sub('[_.]R1[_.]', '_', path)
            dst = re.sub('[-_.]fastq', '_1.fastq', dst)
        elif re.search('[_.]R2[_.]', path):
            dst = re.sub('[_.]R2[_.]', '_', path)
            dst = re.sub('[-_.]fastq', '_2.fastq', dst)
        else:
            raise click.UsageError('Fastq names missing `_R1_` or `_R2_`')

        return join(
            self.cache[pk]['data_dir'],
            self.cache[pk]['workflow']['system_id'] + '_' + basename(dst))

    def _echo_summary(self):
        """Get a summary of the matched, skipped, and missing files."""
        skipped, missing, matched, total_matched, nl = [], [], [], 0, '\n'

        for i in self.cache.values():
            if i['workflow']['raw_data']:
                msg = click.style(f"skipped {i['uid']}\t", fg='cyan')
                msg += f"{i['workflow']['raw_data']} exists"
                skipped.append(msg)
            elif i['src_dst_tuples']:
                msg = click.style(f"found {i['uid']}\t\t", fg='green')
                msg += '\n\t\t'.join(src for src, _ in i['src_dst_tuples'])
                total_matched += len(i['src_dst_tuples'])
                matched.append(msg)
            else:
                msg = click.style(f"missing {i['uid']}\t", fg='red')
                msg += 'no files matched'
                missing.append(msg)

        click.echo(
            f"{nl.join([nl] + skipped) if skipped else ''}"
            f"{nl.join([nl] + missing) if missing else ''}"
            f"{nl.join([nl] + matched) if matched else ''}"
            f'\ntotal samples: {len(self.cache)}'
            f'\nsamples skipped: {len(skipped)}'
            f'\nsamples missing: {len(missing)}'
            f'\nsamples matched: {len(matched)}'
            f'\ntotal files matched: {total_matched}')
