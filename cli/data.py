from os.path import join
from os.path import basename
import os
import re
import click
import collections

from cached_property import cached_property

from cli import api
from cli import user_settings
from cli.utils import force_symlink


def import_bedfile():
    # we need to make sure bedfiles don't have the chr prefix
    bedfiles = []

    for i in bedfiles:
        with open(i, 'r') as f:
            if f.readline().startswith('chr'):
                msg = "Bedfiles can't have 'chr' prefixes: " + i
                raise click.UsageError(msg)


class LocalDataImporter():

    def __init__(self, filters, directories, dtype, key):
        key = key if key else lambda x: x['system_id']
        self.dtype = dtype
        self.filters = filters
        self.directories = directories
        self.workflows = collections.defaultdict(dict)
        self.patterns = None

        for i in api.get_instances(verbose=True, **filters):
            self.workflows[i['pk']]['data']: i
            self.workflows[i['pk']]['identifier']: key(i)
            self.workflows[i['pk']]['result']: []

    @cached_property
    def regex(self):
        """Compile regex, see http://stackoverflow.com/questions/8888567."""
        patterns = []

        for pk, i in self.workflows.items():
            if not i['raw_data']:  # skip workflows that already have data
                pattern = re.sub(r'[-_.]', r'[-_.]', i['identifier'])
                pattern = r'(?P<{}>[-_.]?{}[-_.])'.format(pk, pattern)
                patterns.append(pattern)

        return re.compile('|'.join(patterns))

    def import_raw_data(self, link=False, commit=False):
        _check_admin()
        label = f'Exploring {len(self.directories)} directories...'

        for directory in self.directories:
            with click.progressbar(os.walk(directory), label=label) as walk:
                for root, _, files in walk:
                    for i in files:
                        self.process_path(join(root, i))

        return

    def process_path(self, path):
        try:
            pk = next(finditer).lastgroup
        except StopIteration:
            return

        try:
            # raise error if multiple samples match the same file.
            id_1 = self.workflows[pk]['data']['system_id']
            id_1_matched = self.workflows[pk]['identifier']
            pk_2 = next(finditer).lastgroup
            id_2 = self.workflows[pk_2]['data']['system_id']
            id_2_matched = self.workflows[pk_2]['identifier']

            if id_1 != id_1_matched:
                id_1 = f'{id_1} ({id_1_matched})'
                id_2 = f'{id_2} ({id_2_matched})'

            raise click.UsageError(f'{id_1} and {id_2} matched: {path}')
        except StopIteration:
            pass

        if self.dtype == 'fastq':
            dst = self.get_bam_dst(pk, path)
        elif self.dtype == 'bam':
            dst = self.get_fastq_dst(pk, path)

        self.workflows[pk]['result'].append((path, dst))

    def get_bam_dst(self, pk, path):
        if self.is_valid_bam(path):
            suffix = path.split(".bam")[-1]
            dst = join(w.bamdir, str(w.leukid) + ".bam" + suffix)

    def get_fastq_dst(self, pk, path):
        if self.is_valid_fastq(path):
            if re.search("_[1|2].fastq", path):
                dst = path
            elif re.search("[_.]R1[_.]", path):
                dst = re.sub("[_.]R1[_.]", "_", path)
                dst = re.sub("[-_.]fastq", "_1.fastq", dst)
            elif re.search("[_.]R2[_.]", path):
                dst = re.sub("[_.]R2[_.]", "_", path)
                dst = re.sub("[-_.]fastq", "_2.fastq", dst)
            else:
                msg = "Fastq file names must contain `_R1_` or `_R2_`."
                raise click.ClickException(msg)

            dst = "_".join([w.slug, dst])
            dst = join(w.fastqdir, dst)
            matched[w.slug]["result"].append((src, dst))
            return dst

    def is_valid_fastq(self, path):
        return path.endswith('.fastq.gz')

    def is_valid_bam(self, path):
        return path.endswith('.bam') or '.bam.' in path

    def _matched_summary(self, matched, id_type):
        """Get a summary of the matched, skipped, and missing files."""
        skipped, missing, found, imported, files = [], [], [], [], 0
        summary = "\ntotal samples: " + str(len(matched))
        summary += "\nsamples skipped: %s\nsamples missing: %s\nsamples found: %s"
        summary += "\nfiles found: %s"

        for slug in matched:
            w = matched[slug]["workflow"]

            if matched[slug]["result"] == "skipped":
                msg = "skipped %s (using %s) " % (w.leukid, getattr(w, id_type))
                msg = click.style(msg, fg="cyan")
                msg += "Workflow has data."
                skipped.append(msg)

            elif matched[slug]["result"] == "missing":
                msg = "missing %s (using %s) " % (w.leukid, getattr(w, id_type))
                msg = click.style(msg, fg="red")
                msg += "not a single file matched."
                missing.append(msg)

            else:
                imported.append(w)
                msg = "found %s (using %s)\n\t" % (w.leukid, getattr(w, id_type))
                msg = click.style(msg, fg="green")
                msg += "\n\t".join(i[0] for i in matched[w.slug]["result"]) + "\n"
                files += len(matched[w.slug]["result"])
                found.append(msg)

        if skipped: click.echo("\n\n" + "\n".join(skipped))
        if missing: click.echo("\n\n" + "\n".join(missing))
        if found: click.echo("\n\n" + "\n".join(found))
        click.echo(summary % (len(skipped), len(missing), len(found), files))

        return imported
