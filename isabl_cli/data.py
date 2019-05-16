"""Data import logic."""

from collections import defaultdict
from datetime import datetime
from getpass import getuser
from glob import glob
from os.path import dirname
from os.path import basename
from os.path import getsize
from os.path import isdir
from os.path import join
import os
import re
import shutil
import subprocess

from slugify import slugify
import click
import yaml

from isabl_cli import api
from isabl_cli import options
from isabl_cli import utils
from isabl_cli.settings import system_settings
from isabl_cli.settings import import_from_string


def update_experiment_bam_file(experiment, assembly_name, analysis_pk, bam_url):
    """
    Update default bam for a experiment given the assembly.

    Arguments:
        experiment (dict): experiment dict.
        assembly_name (str): assembly name.
        analysis_pk (int): analysis primary key.
        bam_url (str): bam url.

    Returns:
        dict: patched experiment instance
    """
    utils.check_admin()
    pk = experiment["pk"]
    bam_files = experiment["bam_files"]

    if bam_files.get(assembly_name, None):
        raise click.UsageError(f"Experiment {pk} already has {assembly_name} bam")

    bam_files[assembly_name] = {"url": bam_url, "analysis": analysis_pk}
    return api.patch_instance("experiments", pk, bam_files=bam_files)


def symlink_experiment_to_projects(experiment):
    """Create symlink from experiment directory and projects directories."""
    for i in experiment["projects"]:
        if not i["storage_url"]:  # pragma: no cover
            i = update_storage_url("projects", i["pk"])

        experiments_dir = join(i["storage_url"], "experiments")
        experiment_dir = join(experiments_dir, experiment["system_id"])
        os.makedirs(experiments_dir, exist_ok=True)
        utils.force_symlink(experiment["storage_url"], experiment_dir)


def symlink_analysis_to_targets(analysis):
    """Create symlink from experiment directory and projects directories."""
    if analysis["status"] != "SUCCEEDED":
        return

    src = analysis["storage_url"]
    dst = "__".join(
        [
            analysis["application"]["name"].lower().replace(" ", "_"),
            analysis["application"]["version"].lower().replace(" ", "_"),
            str(analysis["pk"]),
        ]
    )

    for i in analysis["targets"]:
        if not i["storage_url"]:  # pragma: no cover
            i = update_storage_url("experiments", i["pk"])

        analyses_dir = join(i["storage_url"], "analyses")
        os.makedirs(analyses_dir, exist_ok=True)
        utils.force_symlink(src, join(analyses_dir, dst))

    if analysis["project_level_analysis"]:
        i = analysis["project_level_analysis"]
        if not i["storage_url"]:
            i = update_storage_url("projects", i["pk"])

        analyses_dir = join(i["storage_url"], "analyses")
        os.makedirs(analyses_dir, exist_ok=True)
        utils.force_symlink(src, join(analyses_dir, dst))


def trigger_analyses_merge(analysis):
    """Submit project level analyses merge if neccessary."""
    if analysis["status"] not in {"SUCCEEDED", "FAILED"}:
        return

    try:
        application = import_from_string(analysis["application"]["application_class"])()
    except ImportError:
        return

    if not hasattr(application.merge_project_analyses, "__isabstractmethod__"):
        projects = {j["pk"]: j for i in analysis["targets"] for j in i["projects"]}

        for i in projects.values():
            pending = api.get_instances_count(
                endpoint="analyses",
                status__in="STARTED,SUBMITTED",
                application=analysis["application"]["pk"],
                projects=i["pk"],
            )

            if not pending:
                application.submit_merge_analysis(i)

    if not hasattr(application.merge_individual_analyses, "__isabstractmethod__"):
        individuals = {
            i.sample.individual.pk: i.sample.individual for i in analysis["targets"]
        }

        for i in individuals.values():
            pending = api.get_instances_count(
                endpoint="analyses",
                status__in="STARTED,SUBMITTED",
                application=analysis["application"]["pk"],
                targets__sample__individual__pk=i.pk,
            )

            if not pending:
                application.submit_merge_analysis(i)


def trash_analysis_storage(analysis):
    """Move analysis `storage_url` to a trash directory."""
    if analysis["status"] == "SUCCEEDED":
        raise click.UsageError("You can't wipe a succeeded analysis")

    if isdir(analysis["storage_url"]):
        trash_dir = system_settings.MAKE_STORAGE_DIRECTORY(
            root=system_settings.BASE_STORAGE_DIRECTORY,
            base=join(".analyses_trash", getuser()),
            identifier=analysis["pk"],
            use_hash=True,
        )

        slug = f'primary_key_{analysis["pk"]}__date_'
        slug += datetime.now(system_settings.TIME_ZONE).isoformat()
        dst = join(trash_dir, slug)
        click.echo(f"\ntrashing: {analysis['storage_url']} -> {dst}\n")
        subprocess.check_call(["mv", analysis["storage_url"], dst])


def _make_storage_directory(root, base, identifier, use_hash=False):
    """
    Get and create path to a data directory.

    The path is set to:

        <root>/<base>/<identifier>

    If `use_hash`, identifier must be integer and path is build using the
    four last digits of the identifier. Say identifiers is 12345, then path is:

        <root>/<base>/23/45/12345

    Arguments:
        base (str): instance's API base.
        identifier (str): instance's primary key.
        root (str): root directory.
        use_hash (bool): hash integer identifier for directories.

    Returns:
        str: path to instance's data directory.
    """
    if not root:  # pragma: no cover
        raise click.UsageError("Base storage root directory is required.")

    if use_hash:
        if not str(identifier).isdigit():  # pragma: no cover
            raise click.UsageError("`use_hash` only supported for integers.")

        hash_1 = f"{identifier:04d}"[-4:-2]
        hash_2 = f"{identifier:04d}"[-2:]
        path = join(base, hash_1, hash_2)
    else:
        path = join(base)

    storage_directory = join(root, path, str(identifier))
    os.makedirs(storage_directory, exist_ok=True, mode=0o770)
    return storage_directory


def get_storage_url(endpoint, identifier, use_hash=False):
    """Make storage directory and return patched instance."""
    return system_settings.MAKE_STORAGE_DIRECTORY(
        root=system_settings.BASE_STORAGE_DIRECTORY,
        base=endpoint,
        identifier=identifier,
        use_hash=use_hash,
    )


def update_storage_url(endpoint, identifier, use_hash=False, **data):
    """Make storage directory and return patched instance."""
    data["storage_url"] = get_storage_url(endpoint, identifier, use_hash)
    return api.patch_instance(endpoint, identifier, **data)


class BaseImporter:
    @staticmethod
    def symlink(src, dst):
        """Create symlink from `src` to `dst`."""
        return utils.force_symlink(os.path.realpath(src), dst)

    @staticmethod
    def move(src, dst):
        """Rename `src` to `dst`."""
        return os.rename(os.path.realpath(src), dst)


class LocalReferenceDataImporter(BaseImporter):

    """An import engine for assemblies' reference_data."""

    @classmethod
    def import_data(
        cls, assembly, species, data_src, data_id, symlink, description, sub_dir=None
    ):
        """
        Register reference resources for a given assembly.

        Arguments:
            assembly (str): name of assembly.
            species (str): name of species.
            data_src (str): path to reference data.
            data_id (str): identifier that will be used for reference data.
            symlink (str): symlink instead of move.
            description (str): reference data description.
            sub_dir (str): target sub dir for the resource, default is data_id.

        Returns:
            dict: updated assembly instance as retrieved from API.
        """
        utils.check_admin()
        data_id = slugify(data_id, separator="_")
        click.echo(f'`data_id` set to: {click.style(data_id, fg="green")}')
        assembly = api.create_instance(
            endpoint="assemblies", name=assembly, species=species
        )

        if data_id in assembly["reference_data"]:
            raise click.UsageError(
                f"Assembly '{assembly['name']}' "
                f"has already reference data registered with id '{data_id}':\n"
                f'\n\t{assembly["reference_data"][data_id]}'
            )

        if not assembly["storage_url"]:
            assembly = update_storage_url("assemblies", assembly["name"])

        data_dir = join(assembly["storage_url"], sub_dir or data_id)
        data_dst = join(data_dir, basename(data_src))
        os.makedirs(data_dir, exist_ok=True)

        if symlink:
            click.echo(f"\nLinking:\n\t{data_src}\n\tto {data_dst}")
            cls.symlink(data_src, data_dst)
        else:
            click.echo(f"\nMoving:\n\t{data_src}\n\tto {data_dst}")
            cls.move(data_src, data_dst)

        click.secho(f'\nSuccess! patching {assembly["name"]}...', fg="green")
        assembly["reference_data"][data_id] = {}
        assembly["reference_data"][data_id]["url"] = data_dst
        assembly["reference_data"][data_id]["description"] = description

        return api.patch_instance(
            endpoint="assemblies",
            identifier=assembly["pk"],
            storage_usage=utils.get_tree_size(assembly["storage_url"]),
            reference_data=assembly["reference_data"],
        )

    @classmethod
    def as_cli_command(cls):
        """Get reference data importer as click command line interface."""
        # build isabl_cli command and return it
        @click.command(name="import-reference-data")
        @click.option("--assembly", help="name of reference genome")
        @click.option("--species", help="species of reference genome")
        @click.option("--description", help="reference data description")
        @click.option("--data-id", help="data identifier (will be slugified)")
        @click.option("--sub-dir", help="target resource sub dir, default is data-id")
        @options.REFERENCE_DATA_SOURCE
        @options.SYMLINK
        def cmd(assembly, data_id, symlink, description, data_src, species, sub_dir):
            """
            Register reference data in the assembly's data directory.

            If you want to import a Reference Genome please refer to the
            import_reference_genome command. Incoming data (files or directories) will
            moved unless `--symlink` is provided.

            Assembly's `storage_url`, `storage_usage` and `reference_data`
            fields are updated, setting the latter to:

                'reference_data': {
                    <data identifier>: {
                        'url': path/to/targets_bedfile.bed,
                        'description': path/to/baits_bedfile.bed
                        },
                    ...
                    }
            """
            cls().import_data(
                data_id=data_id,
                symlink=symlink,
                data_src=data_src,
                assembly=assembly,
                species=species,
                description=description,
                sub_dir=sub_dir,
            )

        return cmd


class LocalReferenceGenomeImporter:

    """An import engine for an assembly reference genome."""

    @classmethod
    def as_cli_command(cls):
        """Get reference data importer as click command line interface."""
        # build isabl_cli command and return it
        @click.command(name="import-reference-genome")
        @click.option("--assembly", help="name of reference genome")
        @click.option("--species", help="species of reference genome")
        @click.option(
            "--genome-path",
            show_default=True,
            type=click.Path(resolve_path=True, exists=True, readable=True),
            help="path to reference genome (.fq, .fasta, .fasta.gz)",
        )
        @click.option(
            "--dont-index",
            is_flag=True,
            show_default=True,
            help="Do not attempt to generate indexes, I'll do it myself",
        )
        @options.SYMLINK
        def cmd(assembly, symlink, genome_path, species, dont_index):
            """
            Register an assembly reference genome.

            By default, an attempt to create indexes will be perfomed.
            """
            assembly = LocalReferenceDataImporter.import_data(
                data_id="genome_fasta",
                symlink=symlink,
                data_src=genome_path,
                assembly=assembly,
                species=species,
                description="Reference Genome Fasta File.",
            )

            genome_fasta = assembly["reference_data"]["genome_fasta"]["url"]
            genome_dir = dirname(genome_fasta)
            commands = [
                ["bwa", "index", genome_fasta],
                ["samtools", "faidx", genome_fasta],
                [
                    "samtools",
                    "dict",
                    genome_fasta,
                    "-a",
                    assembly["name"],
                    "-s",
                    assembly["species"],
                    "-o",
                    join(genome_fasta + ".dict"),
                ],
            ]

            for i in commands:
                if dont_index:
                    click.secho(f"Skipping indexing:\n\n\t{' '.join(i)}", fg="orange")
                    continue

                try:
                    subprocess.check_call(i)
                except subprocess.CalledProcessError:
                    click.secho(
                        f"INDEX FAILED, MUST BE FIXED:\n\n\t{' '.join(i)}", fg="red"
                    )

            indexes = {
                "bwa index": ["amb", "ann", "bwt", "pac", "sa"],
                "samtools faidx": ["fai"],
                "samtools dict": ["dict"],
            }

            for i, indexes in indexes.items():
                for j in indexes:
                    assembly["reference_data"][f"genome_fasta_{j}"] = {
                        "url": join(genome_fasta + f".{j}"),
                        "description": f"Index generated by: {i}",
                    }

            for i in glob(genome_fasta.split(".", 1)[0] + "*"):
                dst = join(genome_dir, assembly["name"] + "." + i.split(".", 1)[-1])
                if i != dst:
                    utils.force_symlink(i, dst)

            api.patch_instance(
                endpoint="assemblies",
                identifier=assembly["pk"],
                storage_usage=utils.get_tree_size(assembly["storage_url"]),
                reference_data=assembly["reference_data"],
            )

        return cmd


class LocalBedImporter:

    """An import engine for techniques' bed_files."""

    @staticmethod
    def process_bedfile(path):
        """Sort, tabix and gzip a bedfile."""
        command = ["sort", "-k1,1V", "-k2,2n", path]
        sorted_bed = subprocess.check_output(command)

        with open(path, "+w") as f:
            f.write(sorted_bed.decode("utf-8"))

        subprocess.check_call(["bgzip", path])
        subprocess.check_call(["tabix", "-p", "bed", path + ".gz"])

        with open(path, "+w") as f:  # write uncompressed file again
            f.write(sorted_bed.decode("utf-8"))

    @classmethod
    def import_bedfiles(
        cls, technique, targets_path, baits_path, assembly, species, description=None
    ):
        """
        Register input_bed_path in technique's storage dir and update `data`.

        Arguments:
            technique (str): technique slug.
            targets_path (str): path to targets bedfile.
            baits_path (str): path to baits bedfile.
            assembly (str): name of reference genome for bedfile.
            species (str): name of genome species.
            description (str): a description of the bed_files.

        Returns:
            dict: updated technique instance as retrieved from API.
        """
        utils.check_admin()
        technique = api.get_instance("techniques", technique)

        if assembly in technique["bed_files"]:
            raise click.UsageError(
                f"Technique '{technique['slug']}' "
                f"has registered bed_files for '{assembly}':\n"
                f'\n\t{technique["bed_files"][assembly]["targets"]}'
                f'\n\t{technique["bed_files"][assembly]["baits"]}'
            )

        if not technique["storage_url"]:
            technique = update_storage_url("techniques", technique["pk"])

        api.create_instance("assemblies", name=assembly, species=species)
        beds_dir = join(technique["storage_url"], "bed_files", assembly)
        base_name = slugify(f'{technique["slug"]}.{assembly}')
        targets_dst = join(beds_dir, f"{base_name}.targets.bed")
        baits_dst = join(beds_dir, f"{base_name}.baits.bed")
        os.makedirs(beds_dir, exist_ok=True)

        for src, dst in [(targets_path, targets_dst), (baits_path, baits_dst)]:
            click.echo(f"\nCopying:\n\t{src}\n\tto {dst}")
            shutil.copy(src, dst)
            click.secho(f"\nProcessing {basename(dst)}...", fg="blue")
            cls.process_bedfile(dst)

        click.secho(f'\nSuccess! patching {technique["slug"]}...', fg="green")
        technique["bed_files"][assembly] = {}
        technique["bed_files"][assembly]["targets"] = targets_dst + ".gz"
        technique["bed_files"][assembly]["baits"] = baits_dst + ".gz"
        technique["bed_files"][assembly]["description"] = description

        return api.patch_instance(
            endpoint="techniques",
            identifier=technique["pk"],
            storage_usage=utils.get_tree_size(technique["storage_url"]),
            bed_files=technique["bed_files"],
        )

    @classmethod
    def as_cli_command(cls):
        """Get bed importer as click command line interface."""
        # build isabl_cli command and return it
        @click.command(name="import-bedfiles")
        @options.TECHNIQUE_IDENTIFIER
        @options.TARGETS_PATH
        @options.BAITS_PATH
        @click.option("--assembly", help="name of reference genome", required=True)
        @click.option("--species", help="name of species", required=True)
        @click.option("--description", help="bed_files description")
        def cmd(technique, assembly, targets_path, baits_path, description, species):
            """
            Register targets and baits bed_files in technique's data directory.

            Incoming bed_files will be compressed and tabixed.
            Both gzipped and uncompressed versions are kept.

            Instance's storage_url, storage_usage and bed_files fields
            are updated, setting the latter to:

            \b
                'bed_files': {
                    <assembly name>: {
                        'targets': path/to/targets_bedfile.bed,
                        'baits': path/to/baits_bedfile.bed
                        },
                    ...
                    }
            """
            cls().import_bedfiles(
                technique=technique,
                targets_path=targets_path,
                baits_path=baits_path,
                assembly=assembly,
                species=species,
                description=description,
            )

        return cmd


class LocalDataImporter(BaseImporter):

    """
    A Data import engine for experiments.

    Attributes:
        FASTQ_REGEX (str): a regex pattern used to match fastq files.
        BAM_REGEX (str): a regex pattern to match bams.
        CRAM_REGEX (str): a regex pattern to match crams.
    """

    BAM_REGEX = r"\.bam$"
    CRAM_REGEX = r"\.cram$"
    FASTQ_REGEX = (
        (r"(([_.]R{0}[_.].+)|([_.]R{0}\.)|(_{0}\.))f(ast)?q(\.gz)?$", "R"),
        (r"(([_.]I{0}[_.].+)|([_.]I{0}\.))f(ast)?q(\.gz)?$", "I"),
    )

    def import_data(
        self,
        directories,
        symlink=False,
        commit=False,
        key=lambda x: x["system_id"],
        files_data=None,
        **filters,
    ):
        """
        Import raw data for multiple experiments.

        Experiments's `storage_url`, `storage_usage`, `sequencing_data` are
        updated.

        Arguments:
            directories (list): list of directories to be recursively explored.
            symlink (bool): if True symlink instead of moving.
            commit (bool): if True perform import operation.
            key (function): given a experiment dict returns id to match.
            filters (dict): key value pairs to use as API query params.
            files_data (dict): keys are files basenames and values are
                dicts with extra annotations such as PL, LB, or any other.

        Raises:
            click.UsageError: if `key` returns the same identifier for multiple
                experiments. If a experiment matches both fastq and bam files.
                if cant determine read 1 or read 2 from matched fastq files.

        Returns:
            tuple: list of experiments for which data has been matched and a
                summary of the operation.
        """
        utils.check_admin()
        files_data = files_data or {}
        experiments_matched = []
        cache = defaultdict(dict)
        patterns = []
        identifiers = {}

        # validate files_data
        for i, j in files_data.items():
            if not isinstance(j, dict):  # pragma: no cover
                raise click.UsageError(f"Invalid file data, expected dict {i}: {j}")

        # get experiments and load cache dictionary
        for i in api.get_instances("experiments", verbose=True, **filters):
            index = f"primary_key_{i['pk']}"
            using_id = f"{i['system_id']} (Skipped, identifier is NULL)"
            identifier = key(i)

            if identifier in identifiers:  # duplicated identifiers not valid
                raise click.UsageError(
                    f"Can't use same identifier for {i['system_id']} "
                    f"and {identifiers[identifier]}: {identifier}"
                )

            if i["sequencing_data"]:
                using_id = f"{i['system_id']} (Skipped, experiment has data)"
            elif identifier:
                identifiers[identifier] = i["system_id"]
                patterns.append(self.get_regex_pattern(index, identifier))
                using_id = f"{i['system_id']} (using {identifier})"

            cache[index]["using_id"] = using_id
            cache[index]["instance"] = i
            cache[index]["files"] = []

        if patterns:
            # see http://stackoverflow.com/questions/8888567 for pattern
            pattern = re.compile("|".join(patterns))
            data_storage_dir = system_settings.BASE_STORAGE_DIRECTORY
            label = f"Exploring directories..."

            # explore dirs
            for directory in set(directories):
                with click.progressbar(os.walk(directory), label=label) as bar:
                    for root, _, files in bar:
                        if not root.startswith(data_storage_dir):
                            for i in files:
                                if len(patterns) > 500:
                                    click.echo(
                                        f"Matching {i} against "
                                        f"{len(patterns)} experiments..."
                                    )

                                path = join(root, i)
                                index = self.match_path(path, pattern)

                                if index:
                                    cache[index]["files"].append(path)

            # process files if needed
            label = "Processing..."
            bar = sorted(cache.values(), key=lambda x: x["instance"]["pk"])
            with click.progressbar(bar, label=label) as bar:
                for i in bar:
                    if commit and i["files"]:
                        experiments_matched.append(
                            self.import_files(
                                instance=i["instance"],
                                files=i["files"],
                                symlink=symlink,
                                files_data=files_data,
                            )
                        )
                    elif i["files"]:  # pragma: no cover
                        experiments_matched.append(i["instance"])

        return experiments_matched, self.get_summary(cache)

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
                for j, _ in self.FASTQ_REGEX:
                    if re.search(j.format(i), path):
                        valid = True

            if re.search(r"\.f(ast)?q(\.gz)?$", path) and not valid:
                msg = f"cant determine fastq type from: {path}"
                raise click.UsageError(msg)

            assert valid
        except (StopIteration, AssertionError):  # pragma: no cover
            return None

        return index

    def import_files(self, instance, files, files_data, symlink):
        """
        Move/link files into instance's `storage_url` and update database.

        Arguments:
            instance (dict): experiment instance.
            files (dict): list of files to be imported.
            symlink (dict): whether to symlink or move the data.
            files_data (dict): keys are files basenames and values are
                dicts with extra annotations such as PL, LB, or any other.

        Raises:
            click.UsageError: if multiple data formats are found.

        Returns:
            dict: patched experiment instance.
        """
        dtypes = set()
        sequencing_data = []
        src_dst = []

        if not instance["storage_url"]:
            instance = update_storage_url(
                endpoint="experiments", identifier=instance["pk"], use_hash=True
            )

        data_dir = join(instance["storage_url"], "data")
        os.makedirs(data_dir, exist_ok=True)

        for src in files:
            file_name = basename(src)
            file_data = files_data.get(file_name, {})

            if re.search(self.BAM_REGEX, src):
                file_type = "BAM"
                dtypes.add(file_type)
            elif re.search(self.CRAM_REGEX, src):
                file_type = "CRAM"
                dtypes.add(file_type)
            else:
                file_type = self.get_fastq_type(file_name)
                dtypes.add("FASTQ")

            # make sure there are no duplicate file names
            if not file_name.startswith(instance["system_id"]):
                file_hash = hex(abs(hash(dirname(src))))[2:]
                file_name = f'{instance["system_id"]}_{file_hash}_{file_name}'

            # make sure we don't add the same file twice
            if all(i != src for i, _ in src_dst):
                dst = join(data_dir, file_name)
                src_dst.append((src, dst))
                sequencing_data.append(
                    dict(
                        file_url=dst,
                        file_type=file_type,
                        file_data=file_data,
                        hash_value=getsize(src),
                        hash_method="os.path.getsize",
                    )
                )

        if len(dtypes) > 1:
            raise click.UsageError(
                "We should have catched this earlier, but multiple formats are "
                f'not supported, these were found for {instance["system_id"]}: '
                f'{",".join(i["file_url"] for i in sequencing_data)}'
            )

        for src, dst in src_dst:
            if symlink:
                self.symlink(src, dst)
            else:
                self.move(src, dst)

        return api.patch_instance(
            endpoint="experiments",
            identifier=instance["pk"],
            storage_url=instance["storage_url"],
            storage_usage=utils.get_tree_size(instance["storage_url"]),
            sequencing_data=sorted(sequencing_data, key=lambda i: i["file_url"]),
        )

    def get_fastq_type(self, file_name):
        """Return destination file name."""
        for index in [1, 2]:
            for pattern, fq_type in self.FASTQ_REGEX:
                if re.search(pattern.format(index), file_name):
                    return f"FASTQ_{fq_type}{index}"

        raise AssertionError(f"Couldn't determine R1, R2 or I1 from {file_name}")

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
        pattern = re.sub(r"[-_. ]", r"[-_. ]", identifier)
        return r"(?P<{}>(^|[-_. ])?{}[-_. ])".format(group_name, pattern)

    @staticmethod
    def get_summary(cache):
        """Get a summary of the matched, skipped, and missing files."""
        skipped, missing, matched, total_matched, nl = [], [], [], 0, "\n"

        for i in cache.values():
            if i["instance"]["sequencing_data"]:
                msg = click.style(f"skipped {i['using_id']}\t", fg="cyan")
                skipped.append(msg)
            elif i["files"]:
                msg = click.style(f"found {i['using_id']}\n\t\t", fg="green")
                total_matched += len(i["files"])
                matched.append(msg + "\n\t\t".join(i["files"]))
            else:
                msg = click.style(f"missing {i['using_id']}\t", fg="red")
                missing.append(msg + "no files matched")

        return (
            f"{nl.join([nl] + skipped) if skipped else ''}"
            f"{nl.join([nl] + missing) if missing else ''}"
            f"{nl.join([nl] + matched) if matched else ''}"
            f"\n\ntotal samples: {len(cache)}"
            f"\nsamples skipped: {len(skipped)}"
            f"\nsamples missing: {len(missing)}"
            f"\nsamples matched: {len(matched)}"
            f"\ntotal files matched: {total_matched}"
        )

    @classmethod
    def as_cli_command(cls):
        """Get data importer as a click command line interface."""
        # build isabl_cli command and return it
        @click.command(name="import-data")
        @options.DIRECTORIES
        @options.IDENTIFIER
        @options.FILTERS
        @options.COMMIT
        @options.SYMLINK
        @options.FILES_DATA
        def cmd(identifier, commit, filters, directories, symlink, files_data):
            """
            Find and import experiments data from many directories.

            Search is recursive and any cram, bam or fastq file that matches
            the experiment identifier will be imported. Multiple data types for
            same experiment is not currently supported.

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
            # function to get nested attributes from dictionary
            def key(experiment):
                value, types = experiment, (int, str, type(None))
                for i in identifier:
                    value = value.get(i)
                if not isinstance(value, types):
                    raise click.UsageError(
                        f"invalid type for identifier "
                        f'`{".".join(identifier)}`: {type(value)}'
                    )
                return value

            if files_data:
                with open(files_data) as f:
                    files_data = yaml.load(f.read())
            else:
                files_data = {}

            matched, summary = cls().import_data(
                directories=directories,
                symlink=symlink,
                commit=commit,
                key=key,
                files_data=files_data,
                **filters,
            )

            click.echo(summary)

            if not commit and matched:  # pragma: no cover
                utils.echo_add_commit_message()

        return cmd
