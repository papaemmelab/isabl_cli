"""Abstract Application."""
# pylint: disable=R0201,C0103

from collections import defaultdict
from contextlib import redirect_stderr
from contextlib import redirect_stdout
from os.path import isdir
from os.path import isfile
from os.path import join
import abc
import os
import sys

from cached_property import cached_property
from click import progressbar
from slugify import slugify
import click

from isabl_cli import api
from isabl_cli import data
from isabl_cli import exceptions
from isabl_cli import utils
from isabl_cli.settings import get_application_settings
from isabl_cli.settings import system_settings


class AbstractApplication:  # pylint: disable=too-many-public-methods

    """An Abstract application."""

    # application unique together definition
    NAME = None
    VERSION = None
    ASSEMBLY = None
    SPECIES = None

    # optional application info
    URL = None

    # application configuration
    application_description = ""
    application_inputs = {}
    application_results = {}
    application_project_level_results = {}
    application_settings = {}
    application_import_strings = {}
    cli_help = ""
    cli_options = []
    skip_status = {"FAILED", "FINISHED", "STARTED", "SUBMITTED", "SUCCEEDED"}
    skip_exceptions = (
        click.UsageError,
        exceptions.MissingRequirementError,
        exceptions.ConfigurationError,
        exceptions.MissingOutputError,
    )

    # private result keys
    _command_script_key = "command_script"
    _command_log_key = "command_log"
    _command_err_key = "command_err"

    # base results
    _base_results = {
        _command_script_key: {
            "frontend_type": "text-file",
            "description": "Script used to execute the analysis.",
            "verbose_name": "Analysis Script",
        },
        _command_log_key: {
            "frontend_type": "text-file",
            "description": "Analysis standard output.",
            "verbose_name": "Standard Output",
        },
        _command_err_key: {
            "frontend_type": "text-file",
            "description": "Analysis standard error.",
            "verbose_name": "Standard Error",
        },
    }

    # -----------------------------
    # USER REQUIRED IMPLEMENTATIONS
    # -----------------------------

    def get_experiments_from_cli_options(self, **cli_options):  # pylint: disable=W9008
        """
        Must return list of target-reference experiment tuples given the parsed options.

        Arguments:
            cli_options (dict): parsed command line options.

        Returns:
            list: of (targets, references) tuples.
        """
        raise NotImplementedError()

    def validate_experiments(self, targets, references):  # pylint: disable=W9008
        """
        Must raise UsageError if tuple combination isn't valid else return True.

        Arguments:
            targets (list): list of targets dictionaries.
            references (list): list of references dictionaries.

        Raises:
            click.UsageError: if tuple is invalid.

        Returns:
            bool: True if (targets, references, analyses) combination is ok.
        """
        raise NotImplementedError()

    def get_dependencies(
        self, targets, references, settings
    ):  # pylint: disable=W9008,W0613
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
            settings (dict): applications settings.

        Returns:
            tuple: command (str), final status (str)
        """
        raise NotImplementedError()

    # -----------------------------
    # USER OPTIONAL IMPLEMENTATIONS
    # -----------------------------

    @abc.abstractmethod  # add __isabstractmethod__ property to method
    def merge_project_analyses(self, analysis, analyses):
        """
        Merge analyses on a project level basis.

        If implemented, a new project level analyses will be created. This
        function will only be called if no other analysis of the same
        application is currently running or submitted.

        Arguments:
            analysis (dict): the project level analysis.
            analyses (list): list of succeeded analyses instances.
        """
        raise NotImplementedError("No merging logic available!")

    def validate_project_analyses(self, project, analyses):
        """Raise AssertionError if project level analysis logic shouldn't happen."""
        return

    def get_analysis_results(self, analysis):  # pylint: disable=W9008,W0613
        """
        Get dictionary of analysis results.

        This function is run on completion.

        Arguments:
            analysis (dict): succeeded analysis instance.

        Returns:
            dict: a jsonable dictionary.
        """
        return {}  # pragma: no cover

    def get_project_analysis_results(self, analysis):  # pylint: disable=W9008,W0613
        """
        Get dictionary of results for a project level analysis.

        This function is run on completion.

        Arguments:
            analysis (dict): succeeded analysis instance.

        Returns:
            dict: a jsonable dictionary.
        """
        return {}  # pragma: no cover

    def get_after_completion_status(self, analysis):  # pylint: disable=W0613
        """Possible values are FINISHED and IN_PROGRESS."""
        return "FINISHED"

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
            endpoint="analyses",
            application=self.application["pk"],
            projects=project["pk"],
            status="SUCCEEDED",
        )

        if analyses:
            error = None

            try:
                self.validate_project_analyses(project, analyses)
                analysis = self.get_project_analysis(project, patch=True)
            except AssertionError as error:
                click.echo(f"Project analysis not created, validation failed: {error}")
                return

            with open(self.get_command_script_path(analysis), "w") as f:
                f.write(
                    f"isabl merge-project-analyses --project {project['pk']} "
                    f"--application {self.primary_key}\n"
                )

            stdout_path = self.get_command_log_path(analysis)
            stderr_path = self.get_command_err_path(analysis)
            oldmask = os.umask(0o22)

            with open(stdout_path, "w") as out, open(stderr_path, "w") as err:
                with redirect_stdout(out), redirect_stderr(err):
                    try:
                        # TODO: setting to submitted here temporarily, will fix later
                        api.patch_analysis_status(analysis, "SUBMITTED")
                        api.patch_analysis_status(analysis, "STARTED")
                        self.merge_project_analyses(analysis, analyses)
                        api.patch_analysis_status(analysis, "SUCCEEDED")
                    except Exception as e:  # pragma: no cover pylint: disable=W0703
                        click.echo(e, file=sys.stderr)
                        api.patch_analysis_status(analysis, "FAILED")
                        error = e

            os.umask(oldmask)
            if error is not None:
                raise Exception(error)

    def submit_project_merge(self, project):
        """
        Directly call `merge_project_analyses` if SUBMIT_PROJECT_LEVEL_MERGE is not set.

        Use setting SUBMIT_PROJECT_LEVEL_MERGE to submit the merge work with custom
        logic, for instance by using LSF with `bsub isabl merge-project-analyses`.

        Arguments:
            project (dict): a project instance.
        """
        if system_settings.SUBMIT_PROJECT_LEVEL_MERGE:
            system_settings.SUBMIT_PROJECT_LEVEL_MERGE(
                project=project, application=self
            )
        else:
            self.run_project_merge(project)

    def get_project_analysis(self, project, patch=False):
        """
        Get or create project level analysis.

        Arguments:
            project (dict): project instance.
            patch (dict): whether or not results and storage_url should be patched.

        Returns:
            dict: analysis instance.
        """
        analysis = api.create_instance(
            endpoint="analyses",
            project_level_analysis=project,
            application=self.project_level_application,
        )

        return self._patch_analysis(analysis) if patch else analysis

    # ----------
    # PROPERTIES
    # ----------

    @cached_property
    def primary_key(self):
        """Space separated name, version and assembly."""
        return self.application["pk"]

    @cached_property
    def settings(self):
        """Return the application settings."""
        defaults = self.application_settings.copy()
        import_strings = set(self.application_import_strings)
        import_strings.add("submit_analyses")

        if "submit_analyses" not in defaults:
            defaults["submit_analyses"] = "isabl_cli.batch_systems.submit_local"

        return get_application_settings(
            defaults=defaults,
            settings=self.application.settings or {},
            reference_data=self.application.assembly.reference_data or {},
            import_strings=import_strings,
        )

    @cached_property
    def application(self):
        """Get application database object."""
        if not all([self.NAME, self.VERSION, self.ASSEMBLY, self.SPECIES]):
            raise NotImplementedError(
                f"NAME must be set: {self.NAME}\n"
                f"VERSION must be set: {self.VERSION}\n"
                f"ASSEMBLY must be set: {self.ASSEMBLY}\n"
                f"SPECIES must be set: {self.SPECIES}\n"
            )

        application = api.create_instance(
            endpoint="applications",
            name=self.NAME,
            version=self.VERSION,
            assembly={"name": self.ASSEMBLY, "species": self.SPECIES},
        )

        api.patch_instance(  # pragma: no cover
            description=self.application_description,
            endpoint="applications",
            identifier=application["pk"],
            application_class=f"{self.__module__}.{self.__class__.__name__}",
            results=self._application_results,
            url=self.URL,
        )

        return application

    @cached_property
    def project_level_application(self):
        """Get or create a project level application database object."""
        application = api.create_instance(
            endpoint="applications",
            name=f"{self.NAME} Project Application",
            version=self.VERSION,
            assembly=self.application["assembly"],
        )

        api.patch_instance(  # pragma: no cover
            description=f"{self.NAME} {self.VERSION} Project Level Application.",
            endpoint="applications",
            identifier=application["pk"],
            results=self._application_project_level_results,
            application_class=self.application["application_class"],
            url=self.URL,
        )

        return application

    @cached_property
    def assembly(self):
        """Get assembly database object."""
        return self.application["assembly"]

    @property
    def _application_results(self):
        ret = self.application_results
        ret.update(self._base_results)
        return ret

    @property
    def _application_project_level_results(self):
        ret = self.application_project_level_results
        ret.update(self._base_results)
        return ret

    # ----------------------------
    # COMMAND LINE INTERFACE LOGIC
    # ----------------------------

    @classmethod
    def as_cli_command(cls):
        """Get application as click command line interface."""
        pipe = cls()

        commit = click.option(
            "--commit",
            show_default=True,
            is_flag=True,
            help="Submit application analyses.",
        )

        verbose = click.option(
            "--verbose",
            show_default=True,
            is_flag=True,
            help="Print verbose output of the operation.",
        )

        force = click.option(
            "--force",
            help="Wipe unfinished analyses and start from scratch.",
            is_flag=True,
        )

        restart = click.option(
            "--restart",
            help="Attempt restarting failed analyses from previous checkpoint.",
            is_flag=True,
        )

        def print_url(ctx, _, value):
            """Print the application url url."""
            if value:
                url = pipe.application["url"]
                if url:
                    click.secho("\nApplication URL:\n", fg="green")
                    click.secho(url)
                else:
                    click.secho("\nApplication has no url defined.\n")
                ctx.exit()

        @click.command(name=pipe.get_cli_command_name(), help=pipe.cli_help)
        @click.option(
            "--url",
            help="Show the url or main repo of the app.",
            is_eager=True,
            is_flag=True,
            expose_value=False,
            callback=print_url,
        )
        @utils.apply_decorators(pipe.cli_options + [commit, force, verbose, restart])
        def command(commit, force, verbose, restart, **cli_options):
            """Click command to be used in the CLI."""
            if commit and force:
                raise click.UsageError("--commit not required when using --force")

            if commit and restart:
                raise click.UsageError("--restart not required when using --force")

            if force and restart:
                raise click.UsageError("cant use --force and --restart together")

            pipe.run(
                tuples=pipe.get_experiments_from_cli_options(**cli_options),
                commit=commit,
                force=force,
                verbose=verbose,
                restart=restart,
                run_args=cli_options,
            )

            if not (commit or force or restart):
                utils.echo_add_commit_message()

        return command

    def get_cli_command_name(self):
        """Get name for isabl_cli command."""
        return f"{slugify(self.NAME)}-{slugify(self.VERSION, separator='.')}"

    # ------------------------
    # ANALYSES EXECUTION LOGIC
    # ------------------------

    def run(
        self, tuples, commit, force=False, restart=False, verbose=True, run_args=None
    ):
        """
        Run a list of targets, references tuples.

        Arguments:
            restart (bool): set settings.restart = True.
            force (bool): if true, analyses are wiped before being submitted.
            commit (bool): if true, analyses are started (`force` overwrites).
            verbose (bool): whether or not verbose output should be printed.
            tuples (list): list of (targets, references) tuples.
            run_args (dict): dictionary of extra arguments required to run the app.

        Returns:
            tuple: command_tuples, skipped_tuples, invalid_tuples
        """
        key = f"{self.NAME} {self.VERSION} {self.ASSEMBLY}"
        commit = True if force or restart else commit
        tuples = [tuple(i) for i in tuples]  # coerce to tuple
        utils.echo_title(f"Running {len(tuples)} tuples for {key}")

        if not tuples:  # pragma: no cover
            return None

        # make sure required settings are set before building commands
        for i in self.application_settings:
            getattr(self.settings, i)

        # set restart attribute
        self.settings.restart = restart

        # update run arguments attribute
        self.settings.run_args = run_args or {}

        # run extra settings validation
        self.validate_settings(self.settings)

        # create and run analyses
        analyses, invalid_tuples = self.get_or_create_analyses(tuples)
        run_tuples, skipped_tuples = self.run_analyses(
            analyses=analyses, commit=commit, force=force, restart=restart
        )

        if verbose:
            self.echo_run_summary(run_tuples, skipped_tuples, invalid_tuples)

        click.echo(
            f"RAN {len(run_tuples)} | "
            f"SKIPPED {len(skipped_tuples)} | "
            f"INVALID {len(invalid_tuples)}\n"
        )

        return run_tuples, skipped_tuples, invalid_tuples

    def run_analyses(self, analyses, commit, force, restart):
        """
        Run a list of analyses.

        Returns:
            list: tuple of
        """
        skipped_tuples = []
        command_tuples = []

        with progressbar(analyses, label="Building commands...\t\t") as bar:
            for i in bar:
                if force and i["status"] not in {"SUCCEEDED", "FINISHED"}:
                    system_settings.TRASH_ANALYSIS_STORAGE(i)
                    os.makedirs(i["storage_url"], exist_ok=True)
                elif (
                    not (restart and i["status"] == "FAILED")
                    and i["status"] in self.skip_status
                ):
                    skipped_tuples.append((i, i["status"]))
                    continue

                try:
                    inputs = i.pop("application_inputs")
                    command = self.get_command(i, inputs, self.settings)
                    command_tuples.append((i, command))
                    self.write_command_script(i, command)
                except self.skip_exceptions as error:  # pragma: no cover
                    skipped_tuples.append((i, error))

        if commit:
            click.echo("Running analyses...")
            run_tuples = self.settings.submit_analyses(self, command_tuples)
        else:
            run_tuples = [(i, "STAGED") for i, _ in command_tuples]
            api.patch_analyses_status([i for i, _ in command_tuples], "STAGED")

        return run_tuples, skipped_tuples

    def get_command_script_path(self, analysis):
        """Get path to analysis command script."""
        return join(analysis["storage_url"], "head_job.sh")

    def get_command_log_path(self, analysis):
        """Get path to analysis log file."""
        return join(analysis["storage_url"], "head_job.log")

    def get_command_err_path(self, analysis):
        """Get path to analysis err file."""
        return join(analysis["storage_url"], "head_job.err")

    def write_command_script(self, analysis, command):
        """
        Write analysis command to bash script and get path to it.

        Arguments:
            analysis (dict): an analysis instance.
            command (str): command to be run.
        """
        # check status, if not run by admin will be marked as succeeded later
        outdir = analysis["storage_url"]
        status = self.get_after_completion_status(analysis)
        assert status in {"IN_PROGRESS", "FINISHED"}, "Status not supported"

        if system_settings.is_admin_user and status == "FINISHED":
            status = "SUCCEEDED"

        # build and write command
        failed = self.get_patch_status_command(analysis["pk"], "FAILED")
        started = self.get_patch_status_command(analysis["pk"], "STARTED")
        finished = self.get_patch_status_command(analysis["pk"], status)

        command = (
            f"umask g+wrx && date && cd {outdir} && "
            f"{started} && {command} && {finished}"
        )

        if system_settings.is_admin_user:
            command += f" && chmod -R a-w {analysis['storage_url']}"

        with open(self.get_command_script_path(analysis), "w") as f:
            template = "{{\n\n    {}\n\n}} | {{\n\n    {}\n\n}}"
            f.write(template.format(command, failed))

    @staticmethod
    def get_patch_status_command(key, status):
        """Return a command to patch the `status` of a given analysis `key`."""
        return f"isabl_cli patch_status --key {key} --status {status}"

    def _get_dependencies(self, targets, references):
        missing = []
        analyses, inputs = self.get_dependencies(targets, references, self.settings)

        for i, j in self.application_inputs.items():  # pragma: no cover
            if i not in inputs and j is NotImplemented:
                missing.append(i)

        if missing:  # pragma: no cover
            missing = ", ".join(map(str, missing))
            raise exceptions.ConfigurationError(
                f"Required inputs missing from `get_dependencies`: {missing}"
            )

        return analyses, inputs

    def _get_analysis_results(self, analysis, created=False):
        """Get results dictionary and append head job script, logs and error files."""
        results = {
            self._command_log_key: self.get_command_log_path(analysis),
            self._command_err_key: self.get_command_err_path(analysis),
            self._command_script_key: self.get_command_script_path(analysis),
        }

        if created:
            return results

        if analysis["project_level_analysis"]:
            specification = self.application_project_level_results
            get_results = self.get_project_analysis_results
        else:
            specification = self.application_results
            get_results = self.get_analysis_results

        results.update(get_results(analysis))

        for i in specification:
            assert i in results, f"Missing expected result {i} in: {results}"

        return results

    # -----------------
    # APPLICATION UTILS
    # -----------------

    @staticmethod
    def get_result(*args, **kwargs):
        """Get an application result."""
        return utils.get_result(*args, **kwargs)

    @staticmethod
    def get_results(*args, **kwargs):
        """Get application results."""
        return utils.get_results(*args, **kwargs)

    def patch_application_settings(self, **settings):
        """Patch application settings if necessary."""
        click.echo(f"Patching settings for {self.NAME} {self.VERSION} {self.ASSEMBLY}")

        try:
            assert self.application["settings"] == settings
            click.secho(f"\n\tno changes detected, skipping patch.\n", fg="yellow")
        except AssertionError:
            try:
                api.patch_instance("applications", self.primary_key, settings=settings)
                click.secho(f"\n\tSuccessfully patched settings.\n", fg="green")
            except TypeError as error:
                click.secho(f"\n\tPatched failed with error: {error}.\n", fg="red")

    @staticmethod
    def get_job_name(analysis):
        """Get job name given an analysis instance."""
        targets = analysis["targets"]
        references = analysis["references"]
        methods = {i["technique"]["method"] for i in targets}
        projects = {str(j["pk"]) for i in targets for j in i["projects"]}

        if len(targets) > 2 or not targets:
            targets = f"{len(targets)} samples."
        else:
            targets = " ".join([i["system_id"] for i in targets])

        if len(references) > 2 or not references:
            references = f"{len(references)} samples."
        else:
            references = " ".join([i["system_id"] for i in references])

        return " | ".join(
            [
                f'analysis: {analysis["pk"]}',
                f"targets: {targets}",
                f"references: {references} |",
                f'methods: {" ".join(methods)}',
                f'projects: {" ".join(projects)}',
                f'rundir: {analysis["storage_url"]}',
                f'application: {analysis["application"]["pk"]}',
            ]
        )

    @staticmethod
    def echo_run_summary(run_tuples, skipped_tuples, invalid_tuples):
        """
        Echo errors for error tuples such as `invalid`, `cant_run`.

        Arguments:
            invalid_tuples (list): of (tuple, error) tuples.
            run_tuples (list): of (analysis, status) tuples.
            skipped_tuples (list): of (analysis, skip reason) tuples.
        """
        cols = ["IDENTIFIER", "PROJECTS", "TARGETS", "REFERENCES", "MESSAGE"]

        summary = ["\n", "\t".join(cols)]
        run_tuples.sort(key=lambda i: str(i[1]))
        invalid_tuples.sort(key=lambda i: str(i[1]))
        skipped_tuples.sort(key=lambda i: str(i[1]))

        def _style_msg(msg):
            msg = str(msg)
            color = "blue"
            blink = False

            if isinstance(msg, Exception):
                color = "red"
            elif msg == "SUCCEEDED":
                color = "green"
            elif msg == "FAILED":
                color = "red"
            elif msg == "STAGED":
                color = "magenta"
                blink = True

            if len(msg) > 20:
                msg = msg[:100] + "..."

            return click.style(msg, fg=color, blink=blink)

        def _style_projects(targets):
            keys = set()
            for i in targets:
                for j in i["projects"]:  # pragma: no cover
                    keys.add(f"{j['pk']} ({i['technique']['method']})")
            return f"{', '.join(keys)}" if keys else "0"

        def _style_experiments(experiments):
            if len(experiments) > 2 or not experiments:
                return f"{len(experiments)}"
            return " ".join([i["system_id"] for i in experiments])

        for i, msg in run_tuples + skipped_tuples + invalid_tuples:
            extra_msg = ""

            if isinstance(i, dict):  # if skipped, succeeded or failed
                identifier = str(i["pk"])
                targets = i["targets"]
                references = i["references"]
                extra_msg = " " + i["storage_url"]
            else:  # if invalid tuples
                identifier = "INVALID"
                targets = i[0]
                references = i[1]

            row = {k: "NA" for k in cols}
            row["IDENTIFIER"] = identifier
            row["MESSAGE"] = _style_msg(msg) + extra_msg
            row["PROJECTS"] = _style_projects(targets)
            row["TARGETS"] = _style_experiments(targets)
            row["REFERENCES"] = _style_experiments(references)
            summary.append("\t".join(row[k] for k in cols))

        summary = "\n".join(summary).expandtabs(20) + "\n"
        click.echo(f"{summary}\n")

    def get_bedfile(self, experiment, bedfile_type="targets"):
        """Get targets or baits bedfile for experiment."""
        return experiment["technique"]["bed_files"][self.ASSEMBLY][bedfile_type]

    def get_bam(self, experiment):
        """Get experiment bam for application assembly."""
        return experiment["bam_files"][self.ASSEMBLY]["url"]

    def get_bams(self, experiments):
        """Get bams for multiple experiments."""
        return [self.get_bam(i) for i in experiments]

    def get_fastq(self, experiment):
        """
        Get experiment fastq R1 and R2 files.

        Arguments:
            experiment (dict): a experiment instance.

        Raises:
            MissingDataError: if number of R1 and R2 files is not the same.

        Returns:
            tuple: list of R1 fastq, list of R2.
        """
        read_1, read_2 = [], []

        for i in experiment["sequencing_data"]:
            if i["file_type"] == "FASTQ_R1":
                read_1.append(i["file_url"])
            elif i["file_type"] == "FASTQ_R2":
                read_2.append(i["file_url"])

        if read_2 and len(read_1) != len(read_2):
            raise exceptions.MissingDataError(
                f"The # of read 1 files ({len(read_1)}) "
                f"and read 2 ({len(read_2)}) should be the same "
                f"for RNA paired-end sequencing, found: {read_1 + read_2}"
            )

        return read_1, read_2

    def update_experiment_bam_file(self, experiment, bam_url, analysis_pk):
        """Update experiment default bam for assembly, ADMIN ONLY."""
        try:
            self.get_bam(experiment)
            return experiment
        except KeyError:
            pass

        return data.update_experiment_bam_file(
            experiment=experiment,
            assembly_name=self.ASSEMBLY,
            analysis_pk=analysis_pk,
            bam_url=bam_url,
        )

    def validate_bams(self, experiments):
        """Raise error not all experiments have registered bams."""
        errors = []
        assembly = self.ASSEMBLY

        for i in experiments:
            try:
                self.get_bam(i)
            except KeyError:
                sample = i["system_id"]
                errors.append(f"{sample} has no registered bam for {assembly}")

        if errors:
            raise exceptions.ValidationError("\n".join(errors))

    def validate_bedfiles(self, experiments, bedfile_type="targets"):
        """Raise error not all experiments have registered bedfile."""
        errors = []

        for i in experiments:
            try:
                bed = self.get_bedfile(i, bedfile_type=bedfile_type)
                assert isfile(bed), f"BED file does not exist: {bed}"
            except (KeyError, AssertionError):
                errors.append(f'{i["system_id"]} has no registered bedfile')

        if errors:
            raise exceptions.ValidationError("\n".join(errors))

    # -----------------------
    # ANALYSES CREATION LOGIC
    # -----------------------

    def get_or_create_analyses(self, tuples):
        """
        Get or create analysis for a application and a list of tuples.

        Arguments:
            tuples (list): list of (targets, references, analyses) tuples.

        Returns:
            tuple: list of analyses, invalid tuples [(tuple, error), ...].
        """
        label = f"Getting analyses for {len(tuples)} tuples...\t\t"

        # add dependencies and inputs
        tuples = [i + self._get_dependencies(*i) for i in tuples]
        existing_analyses, tuples = self.get_existing_analyses(tuples)
        invalid_tuples, created_analyses = [], []

        with click.progressbar(tuples, file=sys.stderr, label=label) as bar:
            for i in bar:
                try:
                    targets, references, analyses, inputs = i
                    self.validate_species(targets + references)
                    self.validate_experiments(targets, references)
                    analysis = self._patch_analysis(
                        api.create_instance(
                            endpoint="analyses",
                            application=self.application,
                            targets=targets,
                            references=references,
                            analyses=analyses,
                        )
                    )

                    analysis["application_inputs"] = inputs
                    created_analyses.append(analysis)
                except (exceptions.ValidationError, AssertionError) as error:
                    invalid_tuples.append((i, exceptions.ValidationError(*error.args)))

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
        projects = {j["pk"] for i in tuples for j in i[0][0]["projects"]}
        cache = defaultdict(list)
        existing, missing = [], []
        filters = dict(application=self.application["pk"], projects__pk__in=projects)

        def get_cache_key(targets, references):
            return (
                tuple(sorted(set(i["pk"] for i in targets))),
                tuple(sorted(set(i["pk"] for i in references))),
            )

        for i in api.get_instances("analyses", **filters):
            cache[get_cache_key(i["targets"], i["references"])].append(i)

        for targets, references, analyses, inputs in tuples:
            expected = set(analyses)
            found_existing = False

            for i in cache[get_cache_key(targets, references)]:
                if expected.issubset(set(i["analyses"])):
                    found_existing = True
                    i["application_inputs"] = inputs
                    existing.append(i)
                    break

            if not found_existing:
                missing.append((targets, references, analyses, inputs))

        return existing, missing

    def _patch_analysis(self, analysis):
        analysis["storage_url"] = data.get_storage_url(
            endpoint="analyses", identifier=analysis["pk"], use_hash=True
        )

        return api.patch_instance(
            "analyses",
            analysis["pk"],
            results=self._get_analysis_results(analysis, created=True),
            storage_url=analysis["storage_url"],
        )

    # -------------------------
    # ANALYSES VALIDATION UTILS
    # -------------------------

    def validate_is_file(self, path):  # pragma: no cover
        """Validate path is file."""
        assert isfile(path), f"{path} is not a file."

    def validate_is_dir(self, path):  # pragma: no cover
        """Validate path is directory."""
        assert isdir(path), f"{path} is not a file."

    def validate_reference_genome(self, reference):
        """Validate genome exists and has required indexes."""
        required = ".fai", ".amb", ".ann", ".bwt", ".pac", ".sa"

        assert all(isfile(reference + i) for i in required), (
            f"Missing indexes please run:\n\n\t"
            f"bwa index {reference}\n\tsamtools faidx {reference}"
        )

        assert isfile(reference + ".dict"), (
            f'Please generate {reference + ".dict"}, e.g.\n\n\t'
            f"samtools dict -a {self.ASSEMBLY} -s {self.SPECIES} "
            f'{reference} > {reference + ".dict"}'
        )

    def validate_has_raw_sequencing_data(self, experiments):
        """Validate experiments have sequencing data."""
        msg = []

        for i in experiments:
            if not i["sequencing_data"]:
                msg.append(f'{i["system_id"]} has no sequencing data...')

        assert not msg, "\n".join(msg)

    def validate_single_data_type(self, experiments):
        """Validate experiments have only one type of sequencing data."""
        self.validate_has_raw_sequencing_data(experiments)
        types = defaultdict(list)

        for i in experiments:
            for j in i["sequencing_data"]:
                file_type = j["file_type"]

                if file_type.startswith("FASTQ_"):
                    file_type = "FASTQ"

                types[file_type].append(i["system_id"])

        assert len(types) == 1, f"Multiple types not supported: {dict(types)}"
        return list(types.keys())[0]

    def validate_fastq_only(self, experiments):
        """Validate sequencing data is only fastq."""
        dtype = self.validate_single_data_type(experiments)
        assert dtype == "FASTQ", f"Only FASTQ supported, found: {dtype}"

    def validate_is_pair(self, targets, references):
        """Validate targets, references tuple is a pair."""
        assert len(targets) == 1 and len(references) == 1, "Pairs only."
        assert (
            targets[0]["pk"] != references[0]["pk"]
        ), "Target can't be same as reference."

    def validate_one_target(self, targets):
        """Validate only one target."""
        assert len(targets) == 1, f"Only 1 target allowed, got: {len(targets)}"

    def validate_one_target_no_references(self, targets, references):
        """Validate only one sample is passed targets and none on references."""
        self.validate_one_target(targets)
        assert not references, f"No reference experiments, got: {len(references)}"

    def validate_at_least_one_target_one_reference(
        self, targets, references
    ):  # pylint: disable=C0103
        """Validate that at least one reference and target are passed."""
        assert targets and references, "References and targets required."

    def validate_targets_not_in_references(self, targets, references):
        """Make sure targets are not passed as references also."""
        refset = set(i["pk"] for i in references)
        template = "%s was also used as reference."
        msg = [template % i["system_id"] for i in targets if i["pk"] in refset]
        assert not msg, "\n".join(msg)

    def validate_methods(self, experiments, methods):
        """Make sure all experiments methods are those expected."""
        msg = []

        for i in experiments:
            if i["technique"]["method"] not in methods:
                msg.append(
                    f"Only '{methods}' sequencing method allowed, "
                    f"found {i['technique']['method']} for {i['system_id']}."
                )

        assert not msg, "\n".join(msg)

    def validate_pdx_only(self, experiments):
        """Make sure experiments come from PDX samples."""
        msg = []

        for i in experiments:
            if not i["is_pdx"]:
                msg.append(f"{i['system_id']} sample is not PDX derived")

        assert not msg, "\n".join(msg)

    def validate_dna_only(self, experiments):
        """Make sure experiments are DNA data."""
        msg = []

        for i in experiments:
            if i["technique"]["analyte"] != "DNA":
                msg.append(f"{i['system_id']} analyte is not DNA")

        assert not msg, "\n".join(msg)

    def validate_rna_only(self, experiments):
        """Make sure experiments are RNA data."""
        msg = []

        for i in experiments:
            if i["technique"]["analyte"] != "RNA":
                msg.append(f"{i['system_id']} analyte is not RNA")

        assert not msg, "\n".join(msg)

    def validate_dna_pairs(self, targets, references):
        """Validate targets, references tuples for base dna applications."""
        self.validate_is_pair(targets, references)
        self.validate_dna_only(targets + references)

    def validate_same_technique(self, targets, references):
        """Validate targets and references have same bedfile."""
        ttec = {i["technique"]["slug"] for i in targets}
        rtec = {i["technique"]["slug"] for i in references}
        assert len(rtec) == 1, f"Expected one technique, got: {rtec}"
        assert len(ttec) == 1, f"Expected one technique, got: {ttec}"
        assert rtec == ttec, f"Same techniques required: {ttec}, {rtec}"

    def validate_species(self, experiments):
        """Validate experiments's species is same as application's setting."""
        msg = []

        for i in experiments:
            if i["sample"]["individual"]["species"] != self.SPECIES:
                msg.append(f'{i["system_id"]} species not supported')

        assert not msg, "\n".join(msg)
