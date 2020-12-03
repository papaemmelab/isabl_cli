"""Data import logic."""

from collections import defaultdict
from datetime import datetime
from getpass import getuser
from glob import glob
from os.path import basename
from os.path import dirname
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
from isabl_cli import exceptions
from isabl_cli import options
from isabl_cli import utils
from isabl_cli.settings import import_from_string
from isabl_cli.settings import system_settings


def raw_data_inspector(path):
    """Determine if a path is a supported raw data file."""
    for i, j in [
        # sequencing
        (r"\.bam$", "BAM"),
        (r"\.cram$", "CRAM"),
        # imaging
        (r"\.png$", "PNG"),
        (r"\.jp(e)?g$", "JPEG"),
        (r"\.tiff$", "TIFF"),
        (r"\.dicom$", "DICOM"),
        # text
        (r"\.tsv(\.gz)?$", "TSV"),
        (r"\.csv(\.gz)?$", "CSV"),
        (r"\.txt(\.gz)?$", "TXT"),
        # other
        (r"\.pdf$", "PDF"),
        (r"\.html$", "HTML"),
        (r"\.md5$", "MD5"),
        (r"\.y[a]?ml$", "YAML"),
    ]:
        if re.search(i, path):
            return j

    for i in [1, 2]:
        for pattern, fq_type in (
            (r"(([_.]R{0}[_.].+)|([_.]R{0}\.)|(_{0}\.))f(ast)?q(\.gz)?$", "R"),
            (r"(([_.]I{0}[_.].+)|([_.]I{0}\.))f(ast)?q(\.gz)?$", "I"),
        ):
            if re.search(pattern.format(i), path):
                return f"FASTQ_{fq_type}{i}"

    if re.search(r"\.f(ast)?q(\.gz)?$", path):
        raise click.UsageError(f"cant determine fastq type from: {path}")

    return None


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

    if bam_files.get(assembly_name, None):  # pragma: no cover
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
        utils.makedirs(experiments_dir)
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
        utils.makedirs(analyses_dir)
        utils.force_symlink(src, join(analyses_dir, dst))

    if analysis["project_level_analysis"]:
        i = analysis["project_level_analysis"]
        if not i["storage_url"]:
            i = update_storage_url("projects", i["pk"])

        analyses_dir = join(i["storage_url"], "analyses")
        utils.makedirs(analyses_dir)
        utils.force_symlink(src, join(analyses_dir, dst))


def trigger_analyses_merge(analysis):
    """Submit project level analyses merge if neccessary."""
    if analysis["status"] not in {"SUCCEEDED", "FAILED"}:
        return

    try:
        application = import_from_string(analysis["application"]["application_class"])()
    except ImportError:
        return

    def _echo_action(instance, application, pending):
        click.secho(
            ("Skipping " if pending else "Submitting ")
            + ("individual " if "species" in instance else "project ")
            + f"merge for {instance} and application {application}",
            fg="green" if pending else "yellow",
        )

    if application.has_project_auto_merge:
        projects = {j["pk"]: j for i in analysis["targets"] for j in i["projects"]}

        for i in projects.values():
            pending = api.get_instances_count(
                endpoint="analyses",
                status__in="STARTED,SUBMITTED",
                application=analysis["application"]["pk"],
                projects=i["pk"],
            )

            _echo_action(i, analysis.application, pending)

            if not pending:
                application.submit_merge_analysis(i)

    if application.has_individual_auto_merge:
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

            _echo_action(i, analysis.application, pending)

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
    original_umask = os.umask(0)
    utils.makedirs(storage_directory)
    os.umask(original_umask)
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
    def echo_src_dst(msg, src, dst):
        """Print a src to dst msg."""
        click.echo(f"\n{msg}:\n\t{src}\n\t{click.style('->', fg='green')} {dst}")

    @staticmethod
    def symlink(src, dst):
        """Create symlink from `src` to `dst`."""
        utils.force_symlink(os.path.realpath(src), dst)

    @staticmethod
    def move(src, dst):
        """Rename `src` to `dst`."""
        os.rename(os.path.realpath(src), dst)

    @staticmethod
    def copy(src, dst):
        """Copy files from `src` to `dst`."""
        shutil.copy(src, dst)


class LocalReferenceDataImporter(BaseImporter):

    """An import engine for assemblies' reference_data."""

    @classmethod
    def import_data(
        cls,
        identifier,
        data_src,
        data_id,
        symlink,
        description,
        sub_dir=None,
        model="assemblies",
    ):
        """
        Register reference resources for a given assembly.

        Arguments:
            identifier (str): name of assembly or technique.
            model (str): either `techniques` or `assemblies`.
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
        instance = api.get_instance(model, identifier)

        if data_id in instance["reference_data"]:
            raise click.UsageError(
                f"{instance['name']} has already reference data registered with id "
                f'"{data_id}":\n\n\t{instance["reference_data"][data_id]}'
            )

        if not instance["storage_url"]:
            instance = update_storage_url(model, instance["name"])

        data_dir = join(instance["storage_url"], sub_dir or data_id)
        data_dst = join(data_dir, basename(data_src))
        utils.makedirs(data_dir)

        if symlink:
            cls.echo_src_dst("Linking", data_src, data_dst)
            cls.symlink(data_src, data_dst)
        else:
            cls.echo_src_dst("Moving", data_src, data_dst)
            cls.move(data_src, data_dst)

        click.secho(f'\nSuccess! patching {instance["name"]}...', fg="green")
        instance["reference_data"][data_id] = {}
        instance["reference_data"][data_id]["url"] = data_dst
        instance["reference_data"][data_id]["description"] = description
        return api.patch_instance(
            endpoint=model,
            instance_id=instance["pk"],
            storage_usage=utils.get_tree_size(instance["storage_url"]),
            reference_data=instance["reference_data"],
        )

    @classmethod
    def as_cli_command(cls):
        """Get reference data importer as click command line interface."""
        # build isabl_cli command and return it
        @click.command(name="import-reference-data")
        @click.option(
            "--identifier",
            help="name of technique or assembly, see --model",
            required=True,
        )
        @click.option(
            "--model",
            help="default model is assemblies",
            type=click.Choice(["assemblies", "techniques"]),
            default="assemblies",
        )
        @click.option("--description", help="reference data description")
        @click.option("--data-id", help="data identifier (will be slugified)")
        @click.option("--sub-dir", help="target resource sub dir, default is data-id")
        @options.REFERENCE_DATA_SOURCE
        @options.SYMLINK
        def cmd(identifier, model, data_id, symlink, description, data_src, sub_dir):
            """
            Register reference data for assemblies (default) or techniques.

            If you want to import a Reference Genome please refer to the
            import-reference-genome command. Incoming data (files or directories) will
            moved unless `--symlink` is provided.
            """
            cls().import_data(
                data_id=data_id,
                data_src=data_src,
                description=description,
                identifier=identifier,
                model=model,
                sub_dir=sub_dir,
                symlink=symlink,
                copy=copy,
            )

        return cmd


class LocalReferenceGenomeImporter:

    """An import engine for an assembly reference genome."""

    @classmethod
    def as_cli_command(cls):
        """Get reference data importer as click command line interface."""
        # build isabl_cli command and return it
        @click.command(name="import-reference-genome")
        @click.option("--assembly", help="name of genome assembly")
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
        def cmd(assembly, symlink, genome_path, dont_index):
            """
            Register an assembly reference genome.

            By default, an attempt to create indexes will be perfomed.
            """
            assembly = LocalReferenceDataImporter.import_data(
                data_id="genome_fasta",
                symlink=symlink,
                data_src=genome_path,
                identifier=assembly,
                model="assemblies",
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
                    click.secho(f"Skipping indexing:\n\n\t{' '.join(i)}", fg="yellow")
                    continue

                try:  # pragma: no cover
                    subprocess.check_call(i)
                except (
                    FileNotFoundError,
                    subprocess.CalledProcessError,
                ):  # pragma: no cover
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
                instance_id=assembly["pk"],
                storage_usage=utils.get_tree_size(assembly["storage_url"]),
                reference_data=assembly["reference_data"],
            )

        return cmd


class LocalBedImporter(BaseImporter):

    """An import engine for Technique BED files."""

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
            description (str): a description of the BED files.

        Returns:
            dict: updated technique instance as retrieved from API.
        """
        utils.check_admin()
        technique = api.get_instance("techniques", technique)
        targets_key = f"{assembly}_targets_bedfile"
        baits_key = f"{assembly}_baits_bedfile"

        if targets_key in technique["reference_data"]:
            raise click.UsageError(
                f"Technique '{technique['slug']}' "
                f"has registered BED files for '{assembly}':\n"
                f'\n\t{technique["reference_data"][targets_key]}'
                f'\n\t{technique["reference_data"][baits_key]}'
            )

        if not technique["storage_url"]:
            technique = update_storage_url("techniques", technique["pk"])

        api.create_instance("assemblies", name=assembly, species=species)
        beds_dir = join(technique["storage_url"], "bed_files", assembly)
        base_name = slugify(f'{technique["slug"]}.{assembly}')
        targets_dst = join(beds_dir, f"{base_name}.targets.bed")
        baits_dst = join(beds_dir, f"{base_name}.baits.bed")
        utils.makedirs(beds_dir)

        for src, dst in [(targets_path, targets_dst), (baits_path, baits_dst)]:
            cls.echo_src_dst("Copying", src, dst)
            shutil.copy(src, dst)
            click.secho(f"\nProcessing {basename(dst)}...", fg="blue")
            cls.process_bedfile(dst)

        click.secho(f'\nSuccess! patching {technique["slug"]}...', fg="green")

        for i, j in [(targets_key, targets_dst), (baits_key, baits_dst)]:
            technique["reference_data"][i] = {
                "url": j + ".gz",
                "description": description,
            }

        return api.patch_instance(
            endpoint="techniques",
            instance_id=technique["pk"],
            storage_usage=utils.get_tree_size(technique["storage_url"]),
            reference_data=technique["reference_data"],
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
        @click.option("--description", help="BED files description")
        def cmd(technique, assembly, targets_path, baits_path, description, species):
            """
            Register targets and baits BED files in a technique's data directory.

            Incoming BED files will be compressed and tabixed. Both gzipped and
            uncompressed versions are kept. Instance's storage_url, storage_usage,
            and reference_data fields are updated.
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
        RAW_DATA_INSPECTORS (list): list of functions which argument is a path and
            are tasked to return the data type of supported formats.
    """

    RAW_DATA_INSPECTORS = [raw_data_inspector]

    def annotate_file_data(
        self, experiment, file_data, file_type, src, dst
    ):  # pylint: disable=no-self-use,unused-argument
        """Overwrite this method to provide custom file_data annotations."""
        return file_data

    def import_data(
        self,
        directories,
        symlink=False,
        commit=False,
        key=lambda x: x["system_id"],
        files_data=None,
        dtypes=None,
        iexact=False,
        ignore_ownership=False,
        copy=False,
        **filters,
    ):
        """
        Import raw data for multiple experiments.

        Experiments's `storage_url`, `storage_usage`, `raw_data` are
        updated.

        Arguments:
            directories (list): list of directories to be recursively explored.
            copy (bool): if True files are copied intead of moved.
            symlink (bool): if True symlink instead of moving.
            commit (bool): if True perform import operation.
            key (function): given a experiment dict returns id to match.
            filters (dict): key value pairs to use as API query params.
            dtypes (list): data types that should be matched (e.g. BAM, PNG. etc.).
            iexact (bool): case insensitive match of identifiers.
            ignore_ownership (bool): raise error if files not owned by admin user.
            files_data (dict): keys are files basenames and values are
                dicts with extra annotations such as PL, LB, or any other,
                see also annotate_file_data.

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
        dtypes = set(dtypes or [])

        # validate files_data
        for i, j in files_data.items():
            if not isinstance(j, dict):  # pragma: no cover
                raise click.UsageError(f"Invalid file data, expected dict {i}: {j}")

        # get experiments and load cache dictionary
        for i in api.get_instances("experiments", verbose=True, **filters):
            index = f"primary_key_{i['pk']}"
            using_id = f"{i['system_id']} (Skipped, identifier is NULL)"

            try:
                identifier = key(i)
            except Exception as error:  # pylint: disable=broad-except
                raise click.UsageError(f"Failed to retieve id using {key}: {error}")

            if identifier in identifiers:  # duplicated identifiers not valid
                raise click.UsageError(
                    f"Can't use same identifier for {i['system_id']} "
                    f"and {identifiers[identifier]}: {identifier}"
                )

            if i["raw_data"] or i["bam_files"]:
                using_id = f"{i['system_id']} (Skipped, experiment has raw data)"
            elif identifier:
                identifiers[identifier] = i["system_id"]
                patterns.append(self.get_regex_pattern(index, identifier, iexact))
                using_id = f"{i['system_id']} (using {identifier})"

            cache[index]["using_id"] = using_id
            cache[index]["instance"] = i
            cache[index]["files"] = []

        # forgive me gods of complexity
        overlaps = []
        for i in identifiers:
            for j in identifiers:
                if i != j and (i in j or j in i):
                    overlaps.append((i, j))

        if overlaps:
            raise click.UsageError(
                "Overlapping identifiers are not allowed, use a more specific field "
                "or further filter your experiments:\n\n"
                + f"\n\n".join(
                    f"  {identifiers[i]}\t{click.style(i, fg='red')}\n"
                    + f"  {identifiers[j]}\t{click.style(j, fg='red')}"
                    for i, j in overlaps
                ).expandtabs(20)
            )

        if patterns:
            # see http://stackoverflow.com/questions/8888567 for pattern
            pattern = re.compile("|".join(patterns))
            data_storage_dir = system_settings.BASE_STORAGE_DIRECTORY
            label = f"Exploring directories..."

            # explore dirs
            for directory in set(directories):
                with click.progressbar(
                    os.walk(directory, followlinks=True), label=label
                ) as bar:
                    for root, _, files in bar:
                        if not root.startswith(data_storage_dir):
                            for i in files:
                                if len(patterns) > 500:  # pragma: no cover
                                    click.echo(
                                        f"Matching {i} against "
                                        f"{len(patterns)} experiments..."
                                    )

                                path = join(root, i)
                                match = self.match_path(path, pattern)

                                if match and (not dtypes or match["dtype"] in dtypes):
                                    cache[match.pop("index")]["files"].append(match)

            # check ownership if needed
            if not ignore_ownership and not symlink and not copy:
                self.check_ownership(cache)

            # check files are readable
            self.check_are_readable(cache)

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
                                copy=copy,
                                files_data=files_data,
                            )
                        )
                    elif i["files"]:  # pragma: no cover
                        experiments_matched.append(i["instance"])

        return experiments_matched, self.get_summary(cache)

    @staticmethod
    def check_ownership(cache):
        """Make sure files matched are owned by current user."""
        label = "Checking ownership, ignore with --ignore-ownership..."
        bar = sorted(cache.values(), key=lambda x: x["instance"]["pk"])
        owner_mismatch = []

        with click.progressbar(bar, label=label) as bar:
            for i in bar:
                for j in i["files"]:
                    try:
                        utils.assert_same_owner(j["path"])
                    except AssertionError:
                        owner_mismatch.append(j["path"])

        if owner_mismatch:
            raise click.UsageError(
                click.style(
                    "The following files are not owned by current user "
                    "(consider using --ignore-ownership):\n\t",
                    fg="red",
                )
                + "\n\t".join(owner_mismatch)
            )

    @staticmethod
    def check_are_readable(cache):
        """Make sure files matched can be accessed and read."""
        label = "Checking files are readable..."
        bar = sorted(cache.values(), key=lambda x: x["instance"]["pk"])
        unreadable_files = []

        with click.progressbar(bar, label=label) as bar:
            for i in bar:
                for j in i["files"]:
                    try:
                        assert os.access(j["path"], os.R_OK)
                    except AssertionError:
                        unreadable_files.append(j["path"])

        if unreadable_files:
            raise click.UsageError(
                click.style(
                    "The following files are not readable by current user:\n\t",
                    fg="red",
                )
                + "\n\t".join(unreadable_files)
            )

    def match_path(self, path, pattern):
        """Match `path` with `pattern` and update cache if fastq or bam."""
        try:
            # first match ids in path
            matches = pattern.finditer(path)
            index = next(matches).lastgroup
            assert index is not None  # happens when pattern is empty

            # determine if valid data type
            dtypes = [i(path) for i in self.RAW_DATA_INSPECTORS]
            dtypes = set(i for i in dtypes if i)
            assert dtypes  # happens if no data type matched

            # raise error if multiple data types matched
            if len(dtypes) != 1:
                raise exceptions.ImplementationError(  # pragma: no cover
                    f"Conflicting data types ({dtypes}) were "
                    f"identified by inspectors: {self.RAW_DATA_INSPECTORS} "
                )

            return dict(index=index, path=path, dtype=dtypes.pop())
        except (StopIteration, AssertionError):  # pragma: no cover
            return None

    def import_files(self, instance, files, files_data, symlink, copy):
        """
        Move/link files into instance's `storage_url` and update database.

        Arguments:
            instance (dict): experiment instance.
            files (dict): list of files to be imported.
            copy (bool): if True files are copied intead of moved.
            symlink (dict): whether to symlink or move the data.
            files_data (dict): keys are files basenames and values are
                dicts with extra annotations such as PL, LB, or any other.

        Raises:
            click.UsageError: if multiple data formats are found.

        Returns:
            dict: patched experiment instance.
        """
        raw_data = []
        src_dst = []

        if not instance["storage_url"]:
            instance = update_storage_url(
                endpoint="experiments", identifier=instance["pk"], use_hash=True
            )

        data_dir = join(instance["storage_url"], "data")
        utils.makedirs(data_dir)

        for src, file_type in [(i["path"], i["dtype"]) for i in files]:
            file_name = basename(src)
            file_data = files_data.get(file_name, {})

            # make sure there are no duplicate file names
            if not file_name.startswith(instance["system_id"]):
                file_hash = hex(abs(hash(dirname(src))))[2:]
                file_name = f'{instance["system_id"]}_{file_hash}_{file_name}'

            # make sure we don't add the same file twice
            if all(i != src for i, _ in src_dst):
                dst = join(data_dir, file_name)
                src_dst.append((src, dst))
                raw_data.append(
                    dict(
                        hash_value=getsize(src),
                        hash_method="os.path.getsize",
                        file_url=dst,
                        file_type=file_type,
                        file_data=self.annotate_file_data(
                            experiment=instance,
                            file_type=file_type,
                            file_data=file_data,
                            src=src,
                            dst=dst,
                        ),
                    )
                )

        for src, dst in src_dst:
            if copy:
                self.copy(src, dst)
            elif symlink:
                self.symlink(src, dst)
            else:
                self.move(src, dst)

            try:
                subprocess.check_call(["chmod", "a-w", dst])
            except subprocess.CalledProcessError:
                pass

        return api.patch_instance(
            endpoint="experiments",
            instance_id=instance["pk"],
            storage_url=instance["storage_url"],
            storage_usage=utils.get_tree_size(instance["storage_url"]),
            raw_data=sorted(raw_data, key=lambda i: i["file_url"]),
        )

    @staticmethod
    def get_regex_pattern(group_name, identifier, iexact=False):
        """
        Get regex pattern for `identifier` group as `group_name`.

        This pattern treats dashes, underscores and dots equally.

        Arguments:
            group_name (str): regex pattern group name.
            identifier (str): identifier to be matched by regex.
            iexact (bool): use case insensitive match.

        Returns:
            str: a regex pattern.
        """
        return r"(?P<{}>{}(^|[-_. ])?{}[-_. ])".format(
            group_name,
            r"(?i)" if iexact else "",
            re.sub(r"[-_. ]", r"[-_. ]", identifier),
        )

    @staticmethod
    def get_summary(cache):
        """Get a summary of the matched, skipped, and missing files."""
        skipped, missing, matched, total_matched, nl = [], [], [], 0, "\n"

        for i in cache.values():
            if i["instance"]["raw_data"]:
                msg = click.style(f"skipped {i['using_id']}\t", fg="cyan")
                skipped.append(msg)
            elif i["files"]:
                total_matched += len(i["files"])
                matched.append(
                    click.style(f"found {i['using_id']}\n\t\t", fg="green")
                    + "\n\t\t".join([f"{j['dtype']} - {j['path']}" for j in i["files"]])
                )
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
        @click.option(
            "--dtypes", help="Limit data types to be imported.", multiple=True
        )
        @click.option(
            "--iexact", help="Case insensitive match of identifiers.", is_flag=True
        )
        @click.option("--ignore-ownership", help="Don't check ownership.", is_flag=True)
        @click.option(
            "--copy", help="Copy files, instead of moving them.", is_flag=True
        )
        def cmd(
            identifier,
            commit,
            filters,
            directories,
            symlink,
            files_data,
            dtypes,
            iexact,
            ignore_ownership,
            copy,
        ):
            """
            Find and import experiments data from many directories.

            A recursive search will match raw data with database identifiers. Use
            --dtypes to limit the data types that will be matched at a given import.

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

            if symlink and copy:
                raise click.UsageError(f"Use either --copy or --symlink.")

            if symlink and ignore_ownership:
                click.secho("--ignore-ownership isnt used when --symlink.", fg="yellow")

            matched, summary = cls().import_data(
                directories=directories,
                symlink=symlink,
                commit=commit,
                key=key,
                files_data=files_data,
                dtypes=dtypes,
                iexact=iexact,
                ignore_ownership=ignore_ownership,
                copy=copy,
                **filters,
            )

            click.echo(summary)

            if not commit and matched:  # pragma: no cover
                utils.echo_add_commit_message()

        return cmd


class LocalYamlDataImporter(LocalDataImporter):
    def import_data_from_yaml(
        self,
        symlink=False,
        commit=False,
        files_data=None,
        ignore_ownership=False,
        **filters,
    ):
        utils.check_admin()
        experiments_matched = []
        files_to_import = []

        # retrieve the experiment
        experiments = api.get_instances("experiments", verbose=True, **filters)
        assert (
            len(experiments) == 1
        ), f"{len(experiments)} experiments retrieved when expected only 1"
        experiment = experiments[0]

        assert (
            experiment["raw_data"] is None and len(experiment["bam_files"]) == 0
        ), f"Experiment already has data"

        # to reuse import_files method, creating a new dictionary w/o absolute paths
        modified_files_data = {}

        # get files to import
        with open(files_data) as file:
            files_data_yaml = yaml.load(file, Loader=yaml.FullLoader)

            for file_to_be_imported in files_data_yaml.keys():
                # verify files listed in yaml actually exist
                assert os.path.exists(
                    file_to_be_imported
                ), f"File {file_to_be_imported} does not exist"

                # check file ownership, if applicable
                if not ignore_ownership and not symlink:
                    utils.assert_same_owner(file_to_be_imported)

                # determine if valid data type
                dtypes = [i(file_to_be_imported) for i in self.RAW_DATA_INSPECTORS]
                dtypes = set(i for i in dtypes if i)
                assert len(dtypes) == 1  # happens if no data type matched
                files_to_import.append(
                    {"path": file_to_be_imported, "dtype": dtypes.pop()}
                )

                modified_files_data[basename(file_to_be_imported)] = files_data_yaml[
                    file_to_be_imported
                ]

        # import experiment if commit is set to true and there are files to import
        if commit and files_to_import:
            experiments_matched.append(
                self.import_files(
                    instance=experiment,
                    files=files_to_import,
                    symlink=symlink,
                    files_data=modified_files_data,
                )
            )
        else:
            experiments_matched = [experiment]

        return self.get_summary(
            files_to_import, experiment, commit, experiments_matched
        )

    @staticmethod
    def get_summary(files, experiment, commit, matched):
        summary = (
            f"\n\nFor experiment '{experiment.system_id}' "
            f"linked to sample '{experiment.sample.identifier}', "
            f"{click.style('IMPORTED ' if commit else 'WOULD HAVE IMPORTED ', fg='cyan', bold=True)}"
        )
        summary += f"the following {len(files)} files:"

        for file in files:
            summary += f"\n\t {click.style('->', fg='yellow', bold=True)} {file}"

        if not commit and matched:
            summary += click.style(
                "\n\nAdd --commit to proceed.\n", fg="green", blink=True
            )

        return summary

    @classmethod
    def as_cli_command(cls):
        """Get data importer as a click command line interface."""

        # build isabl_cli command and return it
        @click.command(name="import-data-from-yaml")
        @options.FILTERS
        @options.COMMIT
        @options.SYMLINK
        @options.FILES_DATA
        @click.option("--ignore-ownership", help="Don't check ownership.", is_flag=True)
        def cmd(commit, filters, symlink, files_data, ignore_ownership):
            """
            Import data into an experiment by specifying a path to a
            files_data yaml file. The files_data yaml file must contain
            absolute file paths.

            files_data.yaml structure:
            The top level key needs to be an absolute path to a file,
            while the values can be whatever data is relevant to the person importing.

            \b
            \b
            Example of files_data yaml structure:
                \b
                absolute/path/to/file/file_name_1.some_extension      <-- top level key
                    key_1: value_1                                    <-- top level value
                    key_2: value_2                                    <-- top level value
                absolute/path/to/file/file_name_2.some_extension      <-- top level key
                    key_1: value_1                                    <-- top level value
                    key_2: value_2                                    <-- top level value
                ...
            """

            # verify yaml file exits
            assert os.path.exists(
                files_data
            ), f"The following files_data yaml path '{files_data}' does not exist."

            if symlink and ignore_ownership:
                click.secho("--ignore-ownership isnt used when --symlink.", fg="yellow")

            summary = cls().import_data_from_yaml(
                symlink=symlink,
                commit=commit,
                files_data=files_data,
                ignore_ownership=ignore_ownership,
                **filters,
            )

            click.echo(summary)

        return cmd
