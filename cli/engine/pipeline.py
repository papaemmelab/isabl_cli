"""Abstract Pipeline."""
# pylint: disable=R0201,C0103

from collections import defaultdict
from os.path import isdir
from os.path import isfile
from os.path import join
from os import environ
import abc
import os
import subprocess
import sys

from cached_property import cached_property
from click import progressbar
from slugify import slugify
import click

from cli import api
from cli import data
from cli import exceptions
from cli import utils
from cli.exceptions import ValidationError
from cli.settings import PipelineSettings
from cli.settings import system_settings


COMMIT = click.option(
    "--commit",
    show_default=True,
    is_flag=True,
    help="Commit results.")


VERBOSE = click.option(
    "--verbose",
    show_default=True,
    is_flag=True,
    help="Verbose output.")


FORCE = click.option(
    "--force",
    help="Wipe all analyses and start from scratch.",
    is_flag=True)


class AbstractPipeline:

    """
    An Abstract pipeline.

    Attributes:
        NAME (object): TODO.
        VERSION (object): TODO.
        ASSEMBLY (object): TODO.
        SPECIES (object): TODO.
        import_strings (object): TODO.
        engine_settings (object): TODO.
        pipeline_settings (object): TODO.
        cli_help (object): TODO.
        cli_options (object): TODO.
        skip_status (object): TODO.
        skip_exceptions (object): TODO.
    """

    NAME = None
    VERSION = None
    ASSEMBLY = None
    SPECIES = None

    import_strings = {}
    engine_settings = {'raise_error': False}
    pipeline_settings = {}
    cli_help = ""
    cli_options = []
    skip_status = {'FAILED', 'FINISHED', 'STARTED', 'SUBMITTED', 'SUCCEEDED'}
    skip_exceptions = (
        click.UsageError,
        exceptions.MissingRequirementError,
        exceptions.ConfigurationError,
        exceptions.MissingOutputError)

    # -----------------------------
    # USER REQUIRED IMPLEMENTATIONS
    # -----------------------------

    def get_tuples(self, **cli_options):  # pylint: disable=W9008
        """
        Must return list of tuples given the parsed options.

        Arguments:
            cli_options (dict): parsed command line options.

        Returns:
            list: of (targets, references, analyses) tuples.
        """
        raise NotImplementedError()

    def validate_tuple(self, targets, references, analyses):  # pylint: disable=W9008
        """
        Must raise UsageError if tuple combination isnt valid else return True.

        Arguments:
            targets (list): list of targets dictionaries.
            references (list): list of references dictionaries.
            analyses (list): list of analyses dictionaries.

        Raises:
            click.UsageError: if tuple is invalid.

        Returns:
            bool: True if (targets, references, analyses) combination is ok.
        """
        raise NotImplementedError()

    def get_command(self, analysis, settings):  # pylint: disable=W9008
        """
        Must return command and final analysis status.

        Allowed final status are SUCCEEDED and IN_PROGRESS.

        Arguments:
            analysis (dict): an analysis object as retrieved from API.
            settings (dict): pipelines settings.

        Returns:
            tuple: command (str), final status (str)
        """
        raise NotImplementedError()

    @abc.abstractmethod  # add __isabstractmethod__ property to method
    def merge_analyses(self, storage_url, analyses):
        """
        Merge analyses on a project level basis.

        If implemented, a new project level analyses will be created. This
        function will only be called if no other analysis of the same pipeline
        is currently running or submitted.

        Arguments:
            storage_url (str): path of the project level analysis directory.
            analyses (list): list of succeeded analyses instances.
        """
        raise NotImplementedError('No merging logic available!')

    def validate_settings(self, settings):
        """Validate settings."""
        return

    def get_status(self, analysis):  # pylint: disable=W0613
        """Possible values are FINISHED and IN_PROGRESS."""
        return 'FINISHED'

    # ------------------------------
    # MERGE BY PROJECT CONFIGURATION
    # ------------------------------

    def merge_project_analyses(self, project):
        """
        Merge analyses on a project level basis.

        Arguments:
            project (dict): project instance.
        """
        analyses = api.get_instances(
            endpoint='analyses',
            pipeline=self.pipeline['pk'],
            projects=project['pk'],
            status='SUCCEEDED')

        if analyses:
            try:
                analysis = self.get_project_level_analysis(project)
                self.merge_analyses(analysis['storage_url'], analyses)
            except NotImplementedError as error:
                raise error
            except Exception as error:
                api.patch_instance('analyses', analysis['pk'], status='FAILED')
                raise error

    def get_project_level_analysis(self, project):
        """
        Get or create project level analysis.

        Arguments:
            project (dict): project instance.

        Returns:
            dict: analysis instance.
        """
        analysis = api.create_instance(
            endpoint='analyses',
            project_level_analysis=project,
            pipeline=self.pipeline)

        if not analysis['storage_url']:
            analysis = data.update_storage_url(
                endpoint='analyses',
                identifier=analysis['pk'],
                use_hash=True)

        return analysis

    def submit_merge_project_analyses(self, project):
        self.merge_project_analyses(project)

    # ----------
    # PROPERTIES
    # ----------

    @cached_property
    def settings(self):
        """Return the pipeline settings."""
        defaults = self.pipeline_settings.copy()

        for i, j in self.engine_settings.items():
            if i in self.pipeline_settings:
                msg = f"System setting '{i}' can't be used, pick other name"
                raise exceptions.ConfigurationError(msg)
            defaults[i] = j

        return PipelineSettings(
            pipeline=self,
            defaults=defaults,
            import_strings=self.import_strings)

    @cached_property
    def pipeline(self):
        """Get pipeline database object."""
        pipeline_class = f'{self.__module__}.{self.__class__.__name__}'

        if not all([self.NAME, self.VERSION, self.ASSEMBLY, self.SPECIES]):
            raise NotImplementedError(
                f"NAME must be set: {self.NAME}\n"
                f"VERSION must be set: {self.VERSION}\n"
                f"ASSEMBLY must be set: {self.ASSEMBLY}\n"
                f"SPECIES must be set: {self.SPECIES}\n")

        pipeline = api.create_instance(
            endpoint='pipelines',
            name=self.NAME,
            version=self.VERSION,
            assembly={'name': self.ASSEMBLY, 'species': self.SPECIES},
            pipeline_class=pipeline_class)

        if pipeline['pipeline_class'] != pipeline_class:
            api.patch_instance(
                endpoint='pipelines',
                identifier=pipeline['pk'],
                pipeline_class=pipeline_class)

        return pipeline

    @cached_property
    def assembly(self):
        """Get assembly database object."""
        return self.pipeline['assembly']

    # ----------------------------
    # COMMAND LINE INTERFACE LOGIC
    # ----------------------------

    @classmethod
    def as_cli_command(cls):
        """Get pipeline as click command line interface."""
        pipe = cls()

        @click.command(help=cls.cli_help, name=pipe.get_cli_command_name())
        @utils.apply_decorators(cls.cli_options + [COMMIT, FORCE, VERBOSE])
        def command(commit, force, verbose, **cli_options):
            """Click command to be used in the CLI."""
            if commit and force:
                raise click.UsageError(
                    '--commit is redundant with --force, simply use --force')

            pipe.run_tuples(
                tuples=pipe.get_tuples(**cli_options),
                commit=commit,
                force=force,
                verbose=verbose)

            if not (commit or force):
                utils.echo_add_commit_message()

        return command

    def get_cli_command_name(self):
        """Get name for cli command."""
        name = f"{self.NAME} {self.VERSION} {self.ASSEMBLY}"
        return slugify(name, separator='_')

    # ------------------------
    # ANALYSES EXECUTION LOGIC
    # ------------------------

    @staticmethod
    def get_patch_status_command(key, status):
        """Return a command to patch the `status` of a given analysis `key`."""
        return f'cli patch_status --key {key} --status {status}'

    def run_tuples(self, tuples, commit, force=False, verbose=True):
        """
        Run a list of tuples.

        Arguments:
            force (bool): if true, analyses are wiped before being submitted.
            commit (bool): if true, analyses are started (`force` overwrites).
            verbose (bool): whether or not verbose output should be printed.
            tuples (list): list of (targets, references, analyses) tuples.
                elements can be objects or identifiers.

        Returns:
            tuple: command_tuples, skipped_tuples, invalid_tuples
        """
        commit = True if force else commit
        tuples = list(tuples)
        utils.echo_title(
            f'Running {len(tuples)} tuples for {self.NAME} ({self.VERSION})')

        if not tuples:  # pragma: no cover
            return None

        # make sure required settings are set before building commands
        for i in self.pipeline_settings:
            getattr(self.settings, i)

        # run extra settings validation
        self.validate_settings(self.settings)

        # create and run analyses
        analyses, invalid_tuples = self.get_or_create_analyses(tuples)
        command_tuples, skipped_tuples = self.run_analyses(
            analyses=analyses,
            commit=commit,
            force=force)

        if verbose:
            self.echo_summary(command_tuples, skipped_tuples, invalid_tuples)

        click.echo(
            f'RAN {len(command_tuples)} | '
            f'SKIPPED {len(skipped_tuples)} | '
            f'INVALID {len(invalid_tuples)}\n')

        return command_tuples, skipped_tuples, invalid_tuples

    def run_analyses(self, analyses, commit, force):
        """Run a list of analyses."""
        skipped_tuples = []
        command_tuples = []

        with progressbar(analyses, label='Building commands...\t\t') as bar:
            for i in bar:
                if force and i['status'] not in {'SUCCEEDED', 'FINISHED'}:
                    system_settings.TRASH_ANALYSIS_STORAGE(i)
                    os.makedirs(i['storage_url'], exist_ok=True)
                elif i['status'] in self.skip_status:
                    skipped_tuples.append((i, i['status']))
                    continue

                try:
                    command = self.build_command_script(i)
                    command_tuples.append((i, command))
                except self.skip_exceptions as error:  # pragma: no cover
                    skipped_tuples.append((i, error))

        if commit:
            command_tuples = self.submit_analyses(command_tuples)

        return command_tuples, skipped_tuples

    def submit_analyses(self, command_tuples):
        """
        Submit pipelines as arrays grouped by the target methods.

        Arguments:
            command_tuples (list): of (analysis, command) tuples.

        Returns:
            list: analysis, status tuples.
        """
        # make sure analyses are in submitted status to avoid
        # merging project level analyses in every success
        identifiers = []
        api.patch_analyses_status([i for i, _ in command_tuples], 'SUBMITTED')
        label = f'Running {len(command_tuples)} analyses...'

        with progressbar(command_tuples, label=label) as bar:
            for i, j in bar:
                identifiers.append(i['pk'])
                log = join(i['storage_url'], 'head_job.log')
                err = join(i['storage_url'], 'head_job.err')

                with open(log, 'w') as stdout, open(err, 'w') as stderr:
                    subprocess.check_call(
                        args=['bash', j],
                        stdout=stdout,
                        stderr=stderr)

        analyses = api.get_instances('analyses', identifiers)
        return [(i, i['status']) for i in analyses]

    def build_command_script(self, analysis):
        """
        Get analysis command as a path to a bash script.

        Arguments:
            analysis (leuktools.Analysis): the analysis object to be run.

        Returns:
            str: analysis command as a path to a bash script.
        """
        template = '{{\n\n    {}\n\n}} | {{\n\n    {}\n\n}}'
        command = self.get_command(analysis, self.settings)
        outdir = analysis['storage_url']
        command_path = join(analysis['storage_url'], 'head_job.sh')

        # check status, if not run by admin will be marked as succeeded later
        status = self.get_status(analysis)
        assert status in {'IN_PROGRESS', 'FINISHED'}

        if system_settings.is_admin_user and status == 'FINISHED':
            status = 'SUCCEEDED'

        # build command
        failed = self.get_patch_status_command(analysis["pk"], "FAILED")
        started = self.get_patch_status_command(analysis["pk"], "STARTED")
        finished = self.get_patch_status_command(analysis["pk"], status)

        command = (
            f'umask g+wrx && date && cd {outdir} && '
            f'{started} && {command} && {finished}')

        if system_settings.is_admin_user:
            command += f" && chmod -R a-w {analysis['storage_url']}"

        # write command
        with open(command_path, "w") as f:
            f.write(template.format(command, failed))

        return command_path

    @staticmethod
    def echo_summary(command_tuples, skipped_tuples, invalid_tuples):
        """
        Echo errors for error tuples such as `invalid`, `cant_run`.

        Arguments:
            invalid_tuples (list): of (tuple, error) tuples.
            command_tuples (list): of (analysis, command) tuples.
            skipped_tuples (list): of (analysis, skip reason) tuples.
        """
        summary_keys = ['PK', 'PROJECTS', 'TARGETS', 'REFERENCES', 'MESSAGE']
        summary = ['\n']
        command_tuples.sort(key=lambda i: str(i[1]))
        invalid_tuples.sort(key=lambda i: str(i[1]))
        skipped_tuples.sort(key=lambda i: str(i[1]))

        def _style_msg(msg):
            color = 'blue'
            if isinstance(msg, Exception):
                color = 'red'
            elif msg == 'SUCCEEDED':
                color = 'green'
            elif msg == 'FAILED':
                color = 'red'
            return click.style(str(msg), fg=color)

        def _style_projects(targets):
            keys = set()
            for i in targets:
                for j in i['projects']:  # pragma: no cover
                    keys.add(f"{j['pk']} {i['technique']['method']}")
            return f"Projects {' '.join(keys)}" if keys else '0 projects'

        def _style_workflows(workflows, relation='targets'):
            if len(workflows) > 2 or not workflows:
                return f'{len(workflows)} {relation}'
            return ' '.join([i['system_id'] for i in workflows])

        for i, msg in command_tuples + skipped_tuples + invalid_tuples:
            if isinstance(i, dict):  # if skipped, command
                identifier = str(i['pk'])
                targets = i['targets']
                references = i['references']
            else:  # if invalid tuples
                identifier = 'INVALID'
                targets = i[0]
                references = i[1]

            row = {k: 'NA' for k in summary_keys}
            row['PK'] = identifier
            row['MESSAGE'] = _style_msg(msg)
            row['PROJECTS'] = _style_projects(targets)
            row['TARGETS'] = _style_workflows(targets)
            row['REFERENCES'] = _style_workflows(references, 'references')
            summary.append('\t'.join(row[k] for k in summary_keys))

        summary = '\n'.join(summary).expandtabs(12) + '\n'
        click.echo(f'{summary}\n')

    # -----------------------
    # ANALYSES CREATION LOGIC
    # -----------------------

    def get_or_create_analyses(self, tuples):
        """
        Get or create analysis for a pipeline and a list of tuples.

        Arguments:
            tuples (list): list of (targets, references, analyses) tuples.

        Returns:
            tuple: list of analyses, invalid tuples [(tuple, error), ...].
        """
        label = f'Getting analyses for {len(tuples)} tuples...\t\t'
        existing_analyses, tuples = self.get_existing_analyses(tuples)
        created_analyses = []
        invalid_tuples = []

        with click.progressbar(tuples, file=sys.stderr, label=label) as bar:
            for i in bar:
                try:
                    targets, references, analyses = i
                    self.validate_species(targets, references)
                    self.validate_tuple(targets, references, analyses)

                    analysis = api.create_instance(
                        endpoint='analyses',
                        pipeline=self.pipeline,
                        targets=targets,
                        references=references,
                        analyses=analyses)

                    created_analyses.append(
                        data.update_storage_url(
                            endpoint='analyses',
                            identifier=analysis['pk'],
                            use_hash=True))
                except exceptions.ValidationError as error:
                    invalid_tuples.append((i, error))

        return existing_analyses + created_analyses, invalid_tuples

    def get_existing_analyses(self, tuples):
        """
        Get list of existing analyses for a list of `tuples`.

        Check targets and references for an exact match. Check whether all the
        requested analyses are in `analysis.analyses`.

        Arguments:
            tuples (list): list of (targets, references, analyses) tuples.

        Returns:
            tuple: list of existing analyses, list of tuples without analysis
        """
        for i in tuples:
            if len(i) == 2:
                print(i)
        click.echo("Checking for existing analyses...", file=sys.stderr)
        projects = {j['pk'] for i, _, __ in tuples for j in i[0]['projects']}
        filters = dict(pipeline=self.pipeline['pk'], projects__pk__in=projects)
        cache = defaultdict(list)
        existing, missing = [], []

        def get_cache_key(targets, references):
            return (
                tuple(sorted(set(i['pk'] for i in targets))),
                tuple(sorted(set(i['pk'] for i in references))))

        for i in api.get_instances('analyses', **filters):
            cache[get_cache_key(i['targets'], i['references'])].append(i)

        for targets, references, analyses in tuples:
            expected = set(analyses)
            found_existing = False

            for i in cache[get_cache_key(targets, references)]:
                if expected.issubset(set(i['analyses'])):
                    found_existing = True
                    existing.append(i)
                    break

            if not found_existing:
                missing.append((targets, references, analyses))

        return existing, missing

    # ---------------------
    # DATA VALIDATION LOGIC
    # ---------------------

    def validate_is_file(self, path):
        """Validate path is file."""
        if not isfile(path):
            raise ValidationError(f'{path} is not a file.')

    def validate_is_dir(self, path):
        """Validate path is directory."""
        if not isdir(path):
            raise ValidationError(f'{path} is not a file.')

    def validate_reference_genome(self, reference):
        """Validate genome exists and has required indexes."""
        required = ".fai", ".amb", ".ann", ".bwt", ".pac", ".sa"

        if not all(isfile(reference + i) for i in required):
            raise ValidationError(
                f'Missing indexes please run:\n\n\t'
                f'bwa index {reference}\n\tsamtools faidx {reference}')

        if not isfile(reference + '.dict'):
            raise ValidationError(
                f'Please generate {reference + ".dict"}, e.g.\n\n\t'
                f'samtools dict -a {self.ASSEMBLY} -s {self.SPECIES} '
                f'{reference} > {reference + ".dict"}')

    def validate_has_raw_sequencing_data(self, targets, references):
        """Validate targets and references have sequencing data."""
        msg = []

        for i in targets + references:
            if not i['sequencing_data']:
                msg.append(f'{i["system_id"]} has no sequencing data...')

        if msg:
            raise ValidationError('\n'.join(msg))

    def validate_single_data_type(self, targets, references):
        """Validate targets and references have sequencing data."""
        self.validate_has_raw_sequencing_data(targets, references)
        types = set()

        for i in targets + references:
            for j in i['sequencing_data']:
                file_type = j['file_type']

                if file_type.startswith('FASTQ_'):
                    file_type = 'FASTQ'

                types.add(file_type)

        if len(types) > 1:
            raise ValidationError(f'Multiple data types not supported: {types}')

        return types.pop()

    def validate_fastq_only(self, targets, references):
        """Validate sequencing data is only fastq."""
        dtype = self.validate_single_data_type(targets, references)

        if dtype != 'FASTQ':
            raise ValidationError(f'Only {dtype} supported, found: {dtype}')

    def validate_tuple_is_pair(self, targets, references):
        """Validate targets, references tuple is a pair."""
        if len(targets) != 1 or len(references) != 1:
            raise ValidationError('Target, reference pairs required.')

    def validate_one_target_no_references(self, targets, references):
        """Test only one sample is passed targets and none on references."""
        if len(targets) != 1 or references:
            raise ValidationError('References not allowed.')

    def validate_at_least_one_target_one_reference(self, targets, references):
        """Validate that at least one reference and target are passed."""
        if not references or not targets:
            raise ValidationError('References and targets required.')

    def validate_targets_not_in_references(self, targets, references):
        """Make sure targets are not passed as references also."""
        refset = set(i['pk'] for i in references)
        msg = []

        for i in targets:
            if i['pk'] in refset:
                msg.append(f"{i['system_id']} was also used as reference.")

        if msg:
            raise ValidationError('\n'.join(msg))

    def validate_methods(self, targets, references, methods):
        """Make sure all targets and references are come from PDX samples."""
        msg = []

        for i in targets + references:
            if i['technique']['method'] not in methods:
                msg.append(
                    f"Only {methods} allowed, "
                    f"found {i['technique']['method']} for {i['system_id']}.")

        if msg:
            raise ValidationError('\n'.join(msg))

    def validate_pdx_only(self, targets, references):
        """Make sure all targets and references are come from PDX samples."""
        msg = []

        for i in targets + references:
            if not i['specimen']['is_pdx']:
                msg.append(f"{i['system_id']} specimen is not PDX derived")

        if msg:
            raise ValidationError('\n'.join(msg))

    def validate_dna_only(self, targets, references):
        """Make sure all targets and references are DNA data."""
        msg = []

        for i in targets + references:
            if i['technique']['analyte'] != 'DNA':
                msg.append(f"{i['system_id']} analyte is not DNA")

        if msg:
            raise ValidationError('\n'.join(msg))

    def validate_rna_only(self, targets, references):
        """Make sure all targets and references are DNA data."""
        msg = []

        for i in targets + references:
            if i['technique']['analyte'] != 'RNA':
                msg.append(f"{i['system_id']} analyte is not RNA")

        if msg:
            raise ValidationError('\n'.join(msg))

    def validate_dna_pairs(self, targets, references):
        """Validate targets, references tuples for base dna pipelines."""
        self.validate_tuple_is_pair(targets, references)
        self.validate_dna_only(targets, references)

    def validate_same_technique(self, targets, references):
        """Validate targets and references have same bedfile."""
        techniques = {i['technique']['slug'] for i in references}

        if len(techniques) > 1:
            raise ValidationError(
                f'Multiple references techniques: {techniques}')

        msg = []
        for i in targets:
            if i['technique']['slug'] not in techniques:
                msg.append(
                    f"References technique differ from {i['system_id']}: "
                    f"{i['technique']['slug']} =! "
                    f"{references[0]['technique']['slug']}.")

        if msg:
            raise ValidationError('\n'.join(msg))

    def validate_species(self, targets, references):
        """Validate targets and references have sequencing data."""
        msg = []

        for i in targets + references:
            if i['specimen']['individual']['species'] != self.SPECIES:
                msg.append(f'{i["system_id"]} species not supported')

        if msg:
            raise ValidationError('\n'.join(msg))
