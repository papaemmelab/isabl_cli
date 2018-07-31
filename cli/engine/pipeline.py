"""Abstract Pipeline."""
# pylint: disable=R0201,C0103

from collections import defaultdict
from os.path import join
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
from cli.settings import PipelineSettings
from cli.settings import system_settings

from .validator import Validator


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


class AbstractPipeline(Validator):

    """
    An Abstract pipeline.

    Whilst this project becomes more mature, the size of `AbstractPipeline` is
    expected to be reduced and split into components. As for now we will leave
    all logic together so that functionality is transparent.

    For instance, the analysis execution logic could be included in a backend,
    same as data validation and analyses retrieval and creation.

    **User Implementations**

    TODO

    **Project Level Merge**

    TODO

    **Command Line Interface**

    TODO

    **Analyses Execution**

    TODO

    **Analyses Creation**

    TODO

    **Data Validation**

    TODO

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

    # pipeline unique together definition
    NAME = None
    VERSION = None
    ASSEMBLY = None
    SPECIES = None

    # utils
    get_result = utils.get_result
    get_results = utils.get_results

    # pipeline configuration
    import_strings = {}
    engine_settings = {}
    pipeline_inputs = {}
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

    def process_cli_options(self, **cli_options):  # pylint: disable=W9008
        """
        Must return list of tuples given the parsed options.

        Arguments:
            cli_options (dict): parsed command line options.

        Returns:
            list: of (targets, references) tuples.
        """
        raise NotImplementedError()

    def validate_workflows(self, targets, references):  # pylint: disable=W9008
        """
        Must raise UsageError if tuple combination isnt valid else return True.

        Arguments:
            targets (list): list of targets dictionaries.
            references (list): list of references dictionaries.

        Raises:
            click.UsageError: if tuple is invalid.

        Returns:
            bool: True if (targets, references, analyses) combination is ok.
        """
        raise NotImplementedError()

    def get_dependencies(self, targets, references, settings):  # pylint: disable=W9008,W0613
        """
        Get dictionary of inputs, this function is run before `get_command`.

        Arguments:
            targets (list): created analysis instance.
            references (list): created analysis instance.
            settings (object): settings namespace.

        Returns:
            tuple: (list of analyses dependencies primary keys, inputs dict).
        """
        return [], {}

    def get_command(self, analysis, inputs, settings):  # pylint: disable=W9008
        """
        Must return command and final analysis status.

        Allowed final status are SUCCEEDED and IN_PROGRESS.

        Arguments:
            analysis (dict): an analysis object as retrieved from API.
            inputs (dict): as returned by `get_dependencies`.
            settings (dict): pipelines settings.

        Returns:
            tuple: command (str), final status (str)
        """
        raise NotImplementedError()

    # -----------------------------
    # USER OPTIONAL IMPLEMENTATIONS
    # -----------------------------

    @abc.abstractmethod  # add __isabstractmethod__ property to method
    def merge_project_analyses(self, storage_url, analyses):
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

    def get_analysis_results(self, analysis):  # pylint: disable=W9008,W0613
        """
        Get dictionary of analysis results.

        This function is run on completion.

        Arguments:
            analysis (dict): succeeded analysis instance.

        Returns:
            dict: a jsonable dictionary.
        """
        return {} # pragma: no cover

    def get_project_analysis_results(self, analysis):  # pylint: disable=W9008,W0613
        """
        Get dictionary of results for a project level analysis.

        This function is run on completion.

        Arguments:
            analysis (dict): succeeded analysis instance.

        Returns:
            dict: a jsonable dictionary.
        """
        return {} # pragma: no cover

    def get_after_completion_status(self, analysis):  # pylint: disable=W0613
        """Possible values are FINISHED and IN_PROGRESS."""
        return 'FINISHED'

    def validate_settings(self, settings):
        """Validate settings."""
        return

    # ----------------------
    # MERGE BY PROJECT LOGIC
    # ----------------------

    def run_project_merge(self, project):
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
                analysis = self.get_project_analysis(project)
                oldmask = os.umask(0o22)
                self.merge_project_analyses(analysis['storage_url'], analyses)
                os.umask(oldmask)
                api.patch_analysis_status(analysis, 'SUCCEEDED')
            except Exception as error:  # pragma: no cover pylint: disable=W0703
                api.patch_analysis_status(analysis, 'FAILED')
                raise error

    def submit_project_merge(self, project):
        """
        Directly call `merge_project_analyses`.

        Overwrite this method if the merge procedure should be submitted using
        a scheduler, see `commands.merge_project_analyses`.

        Arguments:
            project (dict): a project instance.
        """
        self.run_project_merge(project)

    def get_project_analysis(self, project):
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

    # ----------
    # PROPERTIES
    # ----------

    @property
    def key(self):
        """Space separated name, version and assembly."""
        return f'{self.NAME} {self.VERSION} {self.ASSEMBLY}'

    @cached_property
    def settings(self):
        """Return the pipeline settings."""
        defaults = self.pipeline_settings.copy()

        for i, j in self.engine_settings.items():
            if i in self.pipeline_settings:  # pragma: no cover
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

        if pipeline['pipeline_class'] != pipeline_class:  # pragma: no cover
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

            pipe.run(
                tuples=pipe.process_cli_options(**cli_options),
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

    def run(self, tuples, commit, force=False, verbose=True):
        """
        Run a list of targets, references tuples.

        Arguments:
            force (bool): if true, analyses are wiped before being submitted.
            commit (bool): if true, analyses are started (`force` overwrites).
            verbose (bool): whether or not verbose output should be printed.
            tuples (list): list of (targets, references) tuples.

        Returns:
            tuple: command_tuples, skipped_tuples, invalid_tuples
        """
        commit = True if force else commit
        tuples = [tuple(i) for i in tuples]  # coerce to tuple
        utils.echo_title(f'Running {len(tuples)} tuples for {self.key}')

        if not tuples:  # pragma: no cover
            return None

        # make sure required settings are set before building commands
        for i in self.pipeline_settings:
            getattr(self.settings, i)

        # run extra settings validation
        self.validate_settings(self.settings)

        # create and run analyses
        analyses, invalid_tuples = self.get_or_create_analyses(tuples)
        run_tuples, skipped_tuples = self.run_analyses(
            analyses=analyses,
            commit=commit,
            force=force)

        if verbose:
            self.echo_run_summary(run_tuples, skipped_tuples, invalid_tuples)

        click.echo(
            f'RAN {len(run_tuples)} | '
            f'SKIPPED {len(skipped_tuples)} | '
            f'INVALID {len(invalid_tuples)}\n')

        return run_tuples, skipped_tuples, invalid_tuples

    def run_analyses(self, analyses, commit, force):
        """
        Run a list of analyses.

        Returns:
            list: tuple of
        """
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
                    inputs = i.pop('pipeline_inputs')
                    command = self.get_command(i, inputs, self.settings)
                    command_tuples.append((i, command))
                except self.skip_exceptions as error:  # pragma: no cover
                    skipped_tuples.append((i, error))

        if commit:
            run_tuples = self.submit_analyses(command_tuples)
        else:
            run_tuples = [(i, 'STAGED') for i, _ in command_tuples]
            api.patch_analyses_status([i for i, _ in command_tuples], 'STAGED')

        return run_tuples, skipped_tuples

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
        ret = []
        api.patch_analyses_status([i for i, _ in command_tuples], 'SUBMITTED')
        label = f'Running {len(command_tuples)} analyses...'

        with progressbar(command_tuples, label=label) as bar:
            for i, j in bar:
                self.write_command_script(i, j)
                log = join(i['storage_url'], 'head_job.log')
                err = join(i['storage_url'], 'head_job.err')
                api.patch_analysis_status(i, 'STARTED')
                oldmask = os.umask(0o22)
                status = 'SUCCEEDED'

                with open(log, 'w') as stdout, open(err, 'w') as stderr:
                    try:
                        subprocess.check_call(
                            args=j,
                            cwd=i['storage_url'],
                            stdout=stdout,
                            stderr=stderr,
                            shell=True)

                    except subprocess.CalledProcessError:
                        status = 'FAILED'

                os.umask(oldmask)
                i = api.patch_analysis_status(i, status)
                ret.append((i, i['status']))

        return ret

    def write_command_script(self, analysis, command):
        """
        Write analysis command to bash script and get path to it.

        Arguments:
            analysis (dict): an analysis instance.
            command (str): command to be run.

        Returns:
            str: analysis command as a path to a bash script.
        """
        outdir = analysis['storage_url']
        command_path = join(analysis['storage_url'], 'head_job.sh')

        # check status, if not run by admin will be marked as succeeded later
        status = self.get_after_completion_status(analysis)
        assert status in {'IN_PROGRESS', 'FINISHED'}, 'Status not supported'

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
            template = '{{\n\n    {}\n\n}} | {{\n\n    {}\n\n}}'
            f.write(template.format(command, failed))

        return command_path

    @staticmethod
    def get_patch_status_command(key, status):
        """Return a command to patch the `status` of a given analysis `key`."""
        return f'cli patch_status --key {key} --status {status}'

    def _get_dependencies(self, targets, references):
        missing = []
        analyses, inputs = self.get_dependencies(
            targets, references, self.settings)

        for i, j in self.pipeline_inputs.items():  # pragma: no cover
            if i not in inputs and j is NotImplemented:
                missing.append(i)

        if missing:  # pragma: no cover
            missing = ', '.join(map(str, missing))
            raise exceptions.ConfigurationError(
                f'Required inputs missing from `get_dependencies`: {missing}')

        return analyses, inputs

    # --------------
    # PIPELINE UTILS
    # --------------

    @staticmethod
    def echo_run_summary(run_tuples, skipped_tuples, invalid_tuples):
        """
        Echo errors for error tuples such as `invalid`, `cant_run`.

        Arguments:
            invalid_tuples (list): of (tuple, error) tuples.
            run_tuples (list): of (analysis, status) tuples.
            skipped_tuples (list): of (analysis, skip reason) tuples.
        """
        cols = [
            'IDENTIFIER',
            'PROJECTS',
            'TARGETS',
            'REFERENCES',
            'MESSAGE']

        summary = ['\n', '\t'.join(cols)]
        run_tuples.sort(key=lambda i: str(i[1]))
        invalid_tuples.sort(key=lambda i: str(i[1]))
        skipped_tuples.sort(key=lambda i: str(i[1]))

        def _style_msg(msg):
            msg = str(msg)
            color = 'blue'
            blink = False

            if isinstance(msg, Exception):
                color = 'red'
            elif msg == 'SUCCEEDED':
                color = 'green'
            elif msg == 'FAILED':
                color = 'red'
            elif msg == 'STAGED':
                color = 'magenta'
                blink = True

            if len(msg) > 20:
                msg = msg[:100] + '...'

            return click.style(msg, fg=color, blink=blink)

        def _style_projects(targets):
            keys = set()
            for i in targets:
                for j in i['projects']:  # pragma: no cover
                    keys.add(f"{j['pk']} ({i['technique']['method']})")
            return f"{', '.join(keys)}" if keys else '0'

        def _style_workflows(workflows):
            if len(workflows) > 2 or not workflows:
                return f'{len(workflows)}'
            return ' '.join([i['system_id'] for i in workflows])

        for i, msg in run_tuples + skipped_tuples + invalid_tuples:
            extra_msg = ''

            if isinstance(i, dict):  # if skipped, succeeded or failed
                identifier = str(i['pk'])
                targets = i['targets']
                references = i['references']
                extra_msg = ' ' + i['storage_url']
            else:  # if invalid tuples
                identifier = 'INVALID'
                targets = i[0]
                references = i[1]

            row = {k: 'NA' for k in cols}
            row['IDENTIFIER'] = identifier
            row['MESSAGE'] = _style_msg(msg) + extra_msg
            row['PROJECTS'] = _style_projects(targets)
            row['TARGETS'] = _style_workflows(targets)
            row['REFERENCES'] = _style_workflows(references)
            summary.append('\t'.join(row[k] for k in cols))

        summary = '\n'.join(summary).expandtabs(20) + '\n'
        click.echo(f'{summary}\n')

    def get_bedfile(self, workflow, bedfile_type='targets'):
        """Get targets or baits bedfile for workflow."""
        return workflow['technique']['bed_files'][self.ASSEMBLY][bedfile_type]

    def get_bam(self, workflow):
        """Get workflow bam for pipeline assembly."""
        return workflow['bam_files'][self.ASSEMBLY]['url']

    def get_bams(self, workflows):
        """Get bams for multiple workflows."""
        return [self.get_bam(i) for i in workflows]


    def get_fastq(self, workflow):
        """
        Get workflow fastq R1 and R2 files.

        Arguments:
            workflow (dict): a workflow instance.

        Raises:
            MissingDataError: if number of R1 and R2 files is not the same.

        Returns:
            tuple: list of R1 fastq, list of R2.
        """
        read_1, read_2 = [], []

        for i in workflow['sequencing_data']:
            if i['file_type'] == 'FASTQ_R1':
                read_1.append(i['file_url'])
            elif i['file_type'] == 'FASTQ_R2':
                read_2.append(i['file_url'])

        if read_2 and len(read_1) != len(read_2):
            raise exceptions.MissingDataError(
                f'The # of read 1 files ({len(read_1)}) '
                f'and read 2 ({len(read_2)}) should be the same '
                f'for RNA paired-end sequencing, found: {read_1 + read_2}')

        return read_1, read_2

    def update_workflow_bam_file(self, workflow, bam_url, analysis_pk):
        """Update workflow default bam for assembly, ADMIN ONLY."""
        try:
            self.get_bam(workflow)
            return workflow
        except KeyError:
            pass

        return data.update_workflow_bam_file(
            workflow=workflow,
            assembly_name=self.ASSEMBLY,
            analysis_pk=analysis_pk,
            bam_url=bam_url)

    def validate_bams(self, workflows):
        """Raise error not all workflows have registered bams."""
        errors = []
        assembly = self.ASSEMBLY

        for i in workflows:
            try:
                self.get_bam(i)
            except KeyError:
                sample = i["system_id"]
                errors.append(f'{sample} has no registered bam for {assembly}')

        if errors:
            raise exceptions.ValidationError('\n'.join(errors))

    def validate_bedfiles(self, workflows, bedfile_type='targets'):
        """Raise error not all workflows have registered bams."""
        errors = []

        for i in workflows:
            try:
                self.get_bedfile(i, bedfile_type=bedfile_type)
            except KeyError:
                errors.append(f'{i["system_id"]} has no registered bedfile')

        if errors:
            raise exceptions.ValidationError('\n'.join(errors))

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

        # add dependencies and inputs
        tuples = [i + self._get_dependencies(*i) for i in tuples]
        existing_analyses, tuples = self.get_existing_analyses(tuples)
        invalid_tuples, created_analyses = [], []

        with click.progressbar(tuples, file=sys.stderr, label=label) as bar:
            for i in bar:
                try:
                    targets, references, analyses, inputs = i
                    self.validate_species(targets, references)
                    self.validate_workflows(targets, references)

                    analysis = api.create_instance(
                        endpoint='analyses',
                        pipeline=self.pipeline,
                        targets=targets,
                        references=references,
                        analyses=analyses)

                    analysis = data.update_storage_url(
                        endpoint='analyses',
                        identifier=analysis['pk'],
                        use_hash=True)

                    analysis['pipeline_inputs'] = inputs
                    created_analyses.append(analysis)
                except exceptions.ValidationError as error:
                    invalid_tuples.append((i, error))

        return existing_analyses + created_analyses, invalid_tuples

    def get_existing_analyses(self, tuples):
        """
        Get list of existing analyses for a list of `tuples`.

        Check targets and references for an exact match. Check whether all the
        requested analyses are in `analysis.analyses`.

        Arguments:
            tuples (list): (targets, references, analyses, inputs) tuples.

        Returns:
            tuple: list of existing analyses, list of tuples without analysis
        """
        click.echo("Checking for existing analyses...", file=sys.stderr)
        projects = {j['pk'] for i in tuples for j in i[0][0]['projects']}
        filters = dict(pipeline=self.pipeline['pk'], projects__pk__in=projects)
        cache = defaultdict(list)
        existing, missing = [], []

        def get_cache_key(targets, references):
            return (
                tuple(sorted(set(i['pk'] for i in targets))),
                tuple(sorted(set(i['pk'] for i in references))))

        for i in api.get_instances('analyses', **filters):
            cache[get_cache_key(i['targets'], i['references'])].append(i)

        for targets, references, analyses, inputs in tuples:
            expected = set(analyses)
            found_existing = False

            for i in cache[get_cache_key(targets, references)]:
                if expected.issubset(set(i['analyses'])):
                    found_existing = True
                    i['pipeline_inputs'] = inputs
                    existing.append(i)
                    break

            if not found_existing:
                missing.append((targets, references, analyses, inputs))

        return existing, missing
