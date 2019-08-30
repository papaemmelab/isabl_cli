"""Abstract Application."""
# pylint: disable=R0201,C0103,too-many-lines

from collections import defaultdict
from contextlib import redirect_stderr
from contextlib import redirect_stdout
from os.path import isdir
from os.path import isfile
from os.path import join
import abc
import os
import sys
import traceback

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
    STAGED_MSG = "READY FOR SUBMISSION"

    # application configuration
    application_protect_results = True
    application_description = ""
    application_inputs = {}
    application_results = {}
    application_settings = {}
    application_import_strings = {}

    # auto merge configurations
    application_project_level_results = {}
    application_individual_level_results = {}

    # individual level analyses configs
    unique_analysis_per_individual = False

    cli_help = ""
    cli_options = []
    cli_allow_force = True
    cli_allow_restart = True
    skip_status = {"FAILED", "FINISHED", "STARTED", "SUBMITTED", "SUCCEEDED"}
    skip_exceptions = (
        AssertionError,
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

    # ------------------------
    # OPTIONAL IMPLEMENTATIONS
    # ------------------------

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

    def get_after_completion_status(self, analysis):  # pylint: disable=W0613
        """Possible values are FINISHED and IN_PROGRESS."""
        return "FINISHED"

    def validate_settings(self, settings):
        """Validate settings."""
        return

    # --------------------------------------
    # PROJECT LEVEL OPTIONAL IMPLEMENTATIONS
    # --------------------------------------

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

    # --------------------------------------
    # INDIVIDUAL LEVEL OPTIONAL IMPLEMENTATIONS
    # --------------------------------------

    def get_individual_analysis_results(self, analysis):  # pylint: disable=W9008,W0613
        """
        Get dictionary of results for a individual level analysis.

        This function is run on completion.

        Arguments:
            analysis (dict): succeeded analysis instance.

        Returns:
            dict: a jsonable dictionary.
        """
        return {}  # pragma: no cover

    @abc.abstractmethod  # add __isabstractmethod__ property to method
    def merge_individual_analyses(self, analysis, analyses):
        """
        Merge analyses on a individual level basis.

        If implemented, a new individual level analyses will be created. This
        function will only be called if no other analysis of the same
        application is currently running or submitted.

        Arguments:
            analysis (dict): the individual level analysis.
            analyses (list): list of succeeded analyses instances.
        """
        raise NotImplementedError("No merging logic available!")

    def validate_individual_analyses(self, individual, analyses):
        """Raise AssertionError if individual level analysis logic shouldn't happen."""
        return

    # --------------------
    # MERGE ANALYSES LOGIC
    # --------------------

    def submit_merge_analysis(self, instance):
        """
        Directly call merge logic unless SUBMIT_MERGE_ANALYSIS.

        Use setting SUBMIT_MERGE_ANALYSIS to submit the merge work with custom logic.

        Arguments:
            unused-argument (dict): a project or individual instance.
        """
        if system_settings.SUBMIT_MERGE_ANALYSIS:
            system_settings.SUBMIT_MERGE_ANALYSIS(
                instance=instance,
                application=self,
                command=self._get_cli_merge_command(instance),
            )
        elif "species" in instance:
            self.run_individual_merge(instance)
        else:
            self.run_project_merge(instance)

    def _run_analyses_merge(self, instance, analyses):
        merge_analyses = self.merge_project_analyses
        validate_analyses = self.validate_project_analyses
        get_analysis = self.get_project_level_auto_merge_analysis

        if "species" in instance:
            merge_analyses = self.merge_individual_analyses
            validate_analyses = self.validate_individual_analyses
            get_analysis = self.get_individual_level_auto_merge_analysis

        if not analyses or len(analyses) < 2:
            return

        try:
            validate_analyses(instance, analyses)
            analysis = get_analysis(instance)
        except AssertionError as error:  # pragma: no cover
            click.echo(f"Analysis not created, validation failed: {error}")
            return

        if analysis.status == "STARTED":  # pragma: no cover
            click.secho(f"analysis {analysis} is 'STARTED', exiting...", fg="blue")
            return

        if not analysis.storage_url:
            self._patch_analysis(analysis)

        try:
            # raise error if can't write in directory
            with open(self.get_command_script_path(analysis), "w") as f:
                f.write(self._get_cli_merge_command(instance))
        except PermissionError as error:  # pragma: no cover
            api.patch_analysis_status(analysis, "FAILED")
            raise error

        # make sure merge level analyses are group writable
        oldmask = os.umask(0o07)
        stdout_path = self.get_command_log_path(analysis)
        stderr_path = self.get_command_err_path(analysis)
        error = None

        with open(stdout_path, "w") as out, open(stderr_path, "w") as err:
            with redirect_stdout(out), redirect_stderr(err):
                try:
                    # TODO: setting to submitted here temporarily, will fix later
                    api.patch_analysis_status(analysis, "SUBMITTED")
                    api.patch_analysis_status(analysis, "STARTED")
                    merge_analyses(analysis, analyses)
                    api.patch_analysis_status(analysis, "SUCCEEDED")
                except Exception as e:  # pragma: no cover pylint: disable=W0703
                    print(traceback.format_exc())
                    click.echo(traceback.format_exc(), file=sys.stderr)
                    click.echo(e, file=sys.stderr)
                    api.patch_analysis_status(analysis, "FAILED")
                    error = e

        api.patch_instance("analyses", analysis.pk, data=analysis.data)
        os.umask(oldmask)
        if error is not None:
            raise Exception(error)

    def _get_cli_merge_command(self, instance):
        if "species" in instance:
            return (
                f"isabl merge-individual-analyses --individual {instance['pk']} "
                f"--application {self.primary_key}\n"
            )

        return (
            f"isabl merge-project-analyses --project {instance['pk']} "
            f"--application {self.primary_key}\n"
        )

    # ----------------------
    # MERGE BY PROJECT LOGIC
    # ----------------------

    def run_project_merge(self, project):
        """
        Merge analyses on a project level basis.

        Arguments:
            project (dict): project instance.
        """
        self._run_analyses_merge(
            project,
            api.get_instances(
                endpoint="analyses",
                application=self.application["pk"],
                projects=project["pk"],
                status="SUCCEEDED",
            ),
        )

    def get_project_level_auto_merge_analysis(self, project):
        """
        Get or create project level analysis.

        Arguments:
            project (dict): project instance.

        Returns:
            dict: analysis instance.
        """
        return api.create_instance(
            endpoint="analyses",
            project_level_analysis=project,
            application=self.project_level_auto_merge_application,
        )

    # -------------------------
    # MERGE BY INDIVIDUAL LOGIC
    # -------------------------

    def run_individual_merge(self, individual):
        """
        Merge analyses on a individual level basis.

        Arguments:
            individual (dict): individual instance.
        """
        assert not self.unique_analysis_per_individual, (
            "Applications that require a unique analysis per individual "
            "don't support individual"
        )

        self._run_analyses_merge(
            individual,
            api.get_instances(
                endpoint="analyses",
                application=self.application["pk"],
                targets__sample__individual__pk=individual["pk"],
                status="SUCCEEDED",
            ),
        )

    def get_individual_level_auto_merge_analysis(self, individual):
        """
        Get or create individual level analysis.

        Arguments:
            individual (dict): individual instance.

        Returns:
            dict: analysis instance.
        """
        assert not self.unique_analysis_per_individual, (
            "Applications that require a unique analysis per individual "
            "don't support individual auto merge"
        )

        return api.create_instance(
            endpoint="analyses",
            individual_level_analysis=individual,
            application=self.individual_level_auto_merge_application,
        )

    # ----------
    # PROPERTIES
    # ----------

    @cached_property
    def client_id(self):
        """Get current client ID."""
        return system_settings.client.get("pk", "default_client")

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
            defaults["submit_analyses"] = os.getenv(
                "ISABL_DEFAULT_SUBMIT_ANALYSES", "isabl_cli.batch_systems.submit_local"
            )

        return get_application_settings(
            defaults=defaults,
            settings=self._settings_for_client,
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

        if application.settings.get("default_client") is None:
            application.settings["default_client"] = {}

        return api.patch_instance(
            description=self.application_description,
            endpoint="applications",
            instance_id=application["pk"],
            application_class=f"{self.__module__}.{self.__class__.__name__}",
            results=self._application_results,
            url=self.URL,
            settings=application.settings,
        )

    @cached_property
    def individual_level_auto_merge_application(self):
        """Get or create an individual level application database object."""
        assert not hasattr(
            self.merge_individual_analyses, "__isabstractmethod__"
        ), "No logic implemented to merge analyses for an individual..."

        application = api.create_instance(
            endpoint="applications",
            name=f"{self.NAME} Individual Application",
            version=self.VERSION,
            assembly=self.application["assembly"],
        )

        return api.patch_instance(
            description=f"{self.NAME} {self.VERSION} Individual Level Application.",
            endpoint="applications",
            instance_id=application["pk"],
            results=self._application_individual_level_results,
            application_class=self.application["application_class"],
            url=self.URL,
        )

    @cached_property
    def project_level_auto_merge_application(self):
        """Get or create a project level application database object."""
        assert not hasattr(
            self.merge_project_analyses, "__isabstractmethod__"
        ), "No logic implemented to merge project analyses..."

        application = api.create_instance(
            endpoint="applications",
            name=f"{self.NAME} Project Application",
            version=self.VERSION,
            assembly=self.application["assembly"],
        )

        return api.patch_instance(
            description=f"{self.NAME} {self.VERSION} Project Level Application.",
            endpoint="applications",
            instance_id=application["pk"],
            results=self._application_project_level_results,
            application_class=self.application["application_class"],
            url=self.URL,
        )

    @cached_property
    def assembly(self):
        """Get assembly database object."""
        return self.application["assembly"]

    @property
    def has_project_auto_merge(self):
        """Return True if project level auto merge logic is defined."""
        return not hasattr(self.merge_project_analyses, "__isabstractmethod__")

    @property
    def has_individual_auto_merge(self):
        """Return True if individual level auto merge logic is defined."""
        return not hasattr(self.merge_individual_analyses, "__isabstractmethod__")

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

    @property
    def _application_individual_level_results(self):
        ret = self.application_individual_level_results
        ret.update(self._base_results)
        return ret

    @property
    def _settings_for_client(self):
        return self.application.settings.get(self.client_id) or {}

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

        quiet = click.option(
            "--quiet",
            show_default=True,
            is_flag=True,
            help="Don't print verbose output of the operation.",
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

        cli_options = [  # pylint: disable=unused-variable
            (quiet, True),
            (commit, True),
            (force, cls.cli_allow_force),
            (restart, cls.cli_allow_restart),
        ]

        def print_url(ctx, _, value):
            """Print the application url url."""
            if value:
                url = pipe.application["url"]
                if url:
                    click.secho("\nApplication URL:\n", fg="green")
                    click.secho(url)
                else:  # pragma: no cover
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
        @utils.apply_decorators(pipe.cli_options + [i for i, j in cli_options if j])
        def command(commit, quiet, **cli_options):
            """Click command to be used in the CLI."""
            force = cli_options.pop("force", False)
            restart = cli_options.pop("restart", False)

            if commit and force:
                raise click.UsageError("--commit not required when using --force")

            if commit and restart:  # pragma: no cover
                raise click.UsageError("--restart not required when using --force")

            if force and restart:
                raise click.UsageError("cant use --force and --restart together")

            pipe.run(
                tuples=pipe.get_experiments_from_cli_options(**cli_options),
                commit=commit,
                force=force,
                verbose=not quiet,
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

        # set force attribute
        self.settings.force = force

        # update run arguments attribute
        self.settings.run_args = run_args or {}

        # run extra settings validation
        self.validate_settings(self.settings)

        # create analyses
        analyses, invalid_tuples = self.get_or_create_analyses(tuples)

        # make sure outdir is set
        for i in analyses:
            if not i.storage_url:
                self._patch_analysis(i)

        # run analyses
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
                # if not protect results we need to change the status to staged
                if (
                    not self.application_protect_results
                    and commit
                    and i.status == "SUCCEEDED"
                ):
                    api.patch_analysis_status(i, "STAGED")

                # trash analysis and create outdir again
                elif force and i["status"] not in {"SUCCEEDED", "FINISHED"}:
                    system_settings.TRASH_ANALYSIS_STORAGE(i)
                    os.makedirs(i["storage_url"], exist_ok=True)

                # only restart failed analyses
                elif (
                    not (restart and i["status"] == "FAILED")
                    and i["status"] in self.skip_status
                ):
                    skipped_tuples.append((i, i["status"]))
                    continue

                try:
                    inputs = i.pop(
                        "application_inputs",
                        self._get_dependencies(i.targets, i.references),
                    )

                    command = self.get_command(i, inputs, self.settings)
                    command_tuples.append((i, command))
                    self.write_command_script(i, command)
                except self.skip_exceptions as error:  # pragma: no cover
                    skipped_tuples.append((i, error))

        if commit:
            click.echo("Running analyses...")
            run_tuples = self.settings.submit_analyses(self, command_tuples)
        else:
            run_tuples = [(i, self.STAGED_MSG) for i, _ in command_tuples]

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
        os.makedirs(outdir, exist_ok=True)
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

        with open(self.get_command_script_path(analysis), "w") as f:
            f.write(f"{{\n\n    {command}\n\n}} || {{\n\n    {failed}\n\n}}")

    @staticmethod
    def get_patch_status_command(key, status):
        """Return a command to patch the `status` of a given analysis `key`."""
        return f"isabl patch-status --key {key} --status {status}"

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

        if analysis.project_level_analysis:
            specification = self.application_project_level_results
            get_results = self.get_project_analysis_results
        elif analysis.individual_level_analysis and self.has_individual_auto_merge:
            specification = self.application_individual_level_results
            get_results = self.get_individual_analysis_results
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
        assert system_settings.is_admin_user, "Apps can be patched only by admin user."
        click.echo(f"Patching settings: {self.NAME} {self.VERSION} {self.ASSEMBLY}\n")

        try:
            assert self._settings_for_client == settings
            click.secho(f"\tNo changes detected, skipping patch.\n", fg="yellow")
        except AssertionError:
            try:
                api.patch_instance(
                    "applications",
                    self.primary_key,
                    settings={**self.application.settings, self.client_id: settings},
                )

                del self.settings  # make sure cached settings are re-computed
                del self.application  # make sure application is refetched
                click.secho("\tSuccessfully patched settings.\n", fg="green")
            except TypeError as error:  # pragma: no cover
                click.secho(f"\tPatched failed with error: {error}.\n", fg="red")

        # create or update project level application
        if self.has_project_auto_merge:
            assert self.project_level_auto_merge_application
            click.secho(
                "\tSuccessfully updated project level application.\n", fg="magenta"
            )

        # create or update individual level application
        if self.has_individual_auto_merge:
            assert self.individual_level_auto_merge_application
            click.secho(
                "\tSuccessfully updated individual level application.\n", fg="magenta"
            )

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
                f"references: {references}",
                f'methods: {" ".join(methods)}',
                f'projects: {" ".join(projects)}',
                f'rundir: {analysis["storage_url"]}',
                f'application: {analysis["application"]["pk"]}',
            ]
        )

    def echo_run_summary(self, run_tuples, skipped_tuples, invalid_tuples):
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
            elif msg == "INVALID":
                color = "yellow"
            elif msg == "SUBMITTED":
                color = "cyan"
            elif msg == self.STAGED_MSG:
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
                extra_msg = " " + str(msg)
                msg = identifier = "INVALID"
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
        try:
            return experiment["bam_files"][self.ASSEMBLY]["url"]
        except KeyError:
            raise exceptions.ValidationError(
                f"{experiment.system_id} has no registered bam for {self.ASSEMBLY}"
            )

    def get_bams(self, experiments):
        """Get bams for multiple experiments."""
        return [self.get_bam(i) for i in experiments]

    def update_experiment_bam_file(self, experiment, bam_url, analysis_pk):
        """Update experiment default bam for assembly, ADMIN ONLY."""
        try:
            self.get_bam(experiment)
            return experiment
        except exceptions.ValidationError:
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

        for i in experiments:
            try:
                self.get_bam(i)
            except exceptions.ValidationError as error:
                errors.append(error)

        if errors:
            raise exceptions.ValidationError("\n".join(map(str, errors)))

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

    def validate_are_normals(self, experiments):
        """Raise error not all experiments come from NORMAL sample."""
        for i in experiments:
            msg = f"Experiment Sample {i.sample.system_id} is not NORMAL."
            assert i.sample.sample_class == "NORMAL", msg

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
        # add dependencies and inputs
        invalid_tuples, valid_tuples, created_analyses = [], [], []
        unique_individuals = set()

        for i in tuples:
            try:
                valid_tuples.append(i + self._get_dependencies(*i))

                if self.unique_analysis_per_individual:
                    individual = self._get_individual_from_tuple(i[0], i[1])

                    if individual.pk not in unique_individuals:
                        unique_individuals.add(individual.pk)
                    else:  # pragma: no cover
                        raise exceptions.ValidationError(
                            "Another tuple with the same individual has been passed."
                        )
            except (
                exceptions.ValidationError,
                exceptions.ConfigurationError,
                AssertionError,
            ) as error:
                invalid_tuples.append((i, exceptions.ValidationError(*error.args)))

        # get existing analyses from valid tuples
        if self.unique_analysis_per_individual:
            get_existing_analyses = self.get_individual_level_analyses
        else:
            get_existing_analyses = self.get_existing_analyses

        existing_analyses, valid_tuples = get_existing_analyses(valid_tuples)

        if len(valid_tuples) > 1000:  # pragma: no cover
            click.secho(
                f"Attempting to create {len(valid_tuples)} analyses, "
                f"this process might take some time...",
                fg="yellow",
            )

        with click.progressbar(
            valid_tuples,
            file=sys.stderr,
            label=f"Creating analyses for {len(valid_tuples)} tuples...\t\t",
        ) as bar:
            for i in bar:
                try:
                    targets, references, analyses, inputs, individual = i
                    self.validate_species(targets + references)
                    self.validate_experiments(targets, references)
                    analysis = self._patch_analysis(
                        api.create_instance(
                            endpoint="analyses",
                            application=self.application,
                            targets=targets,
                            references=references,
                            analyses=analyses,
                            individual_level_analysis=individual,
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
        if not tuples:  # pragma: no cover
            return [], []

        click.echo("Checking for existing analyses...", file=sys.stderr)
        projects = {j["pk"] for i in tuples for j in i[0][0]["projects"]}
        cache = defaultdict(list)
        existing, missing = [], []
        targets_pks = ",".join(str(j.pk) for i in tuples for j in i[0])
        filters = dict(application=self.application["pk"], projects__pk__in=projects)

        if targets_pks:
            filters["targets__pk__in"] = targets_pks

        def get_cache_key(targets, references):
            return (
                tuple(sorted(set(i["pk"] for i in targets))),
                tuple(sorted(set(i["pk"] for i in references))),
            )

        for i in api.get_instances("analyses", limit=2000, **filters):
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
                missing.append((targets, references, analyses, inputs, None))

        return existing, missing

    def get_individual_level_analyses(self, tuples):
        """Get existing individual level analyses."""
        # add individual to tuples
        tuples = [i + (self._get_individual_from_tuple(*i[:2]),) for i in tuples]
        tuples_map = {i[-1].pk: i for i in tuples}
        existing = {}

        for i in api.get_analyses(
            application=self.application.pk,
            individual_level_analysis__pk__in=",".join(map(str, tuples_map)),
        ):
            individual = i.individual_level_analysis
            current_tuple = tuples_map[individual.pk]

            # make sure we only have one analysis per individual
            assert individual.pk not in existing, f"Multiple analyses for {individual}"
            existing[individual.pk] = i

            # patch analysis if tuples differ
            for ix, key in enumerate(["targets", "references", "analyses"]):
                if {getattr(j, "pk", j) for j in i[key]} != {
                    getattr(j, "pk", j) for j in current_tuple[ix]
                }:
                    existing[individual.pk] = api.patch_instance(
                        "analyses",
                        i.pk,
                        targets=current_tuple[0],
                        references=current_tuple[1],
                        analyses=current_tuple[2],
                    )
                    break

        return list(existing.values()), [i for i in tuples if i[-1].pk not in existing]

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

    def _get_individual_from_tuple(self, targets, references):
        individual = {
            i.sample.individual.pk: i.sample.individual for i in targets + references
        }

        assert len(individual) == 1, "More than one individual passed!"
        return list(individual.values())[0]

    # ---------------
    # SPECIAL METHODS
    # ---------------

    def __repr__(self):
        """Print name, version and assembly."""
        return f"{self.NAME} {self.VERSION} ({self.ASSEMBLY})"

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

    def validate_has_raw_data(self, experiments):
        """Validate experiments have raw data."""
        msg = []

        for i in experiments:
            if not i["raw_data"]:
                msg.append(f'{i["system_id"]} has no raw data...')

        assert not msg, "\n".join(msg)

    def validate_single_data_type(self, experiments):
        """Validate experiments have only one type of raw data."""
        self.validate_has_raw_data(experiments)
        types = defaultdict(list)

        for i in experiments:
            for j in i["raw_data"]:
                file_type = j["file_type"]

                if file_type.startswith("FASTQ_"):
                    file_type = "FASTQ"

                types[file_type].append(i["system_id"])

        assert len(types) == 1, f"Multiple types not supported: {dict(types)}"
        return list(types.keys())[0]

    def validate_fastq_only(self, experiments):
        """Validate raw data is only fastq."""
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
        """Make sure all experimental methods are those expected."""
        msg = []

        for i in experiments:
            if i["technique"]["method"] not in methods:
                msg.append(
                    f"Only '{methods}' method(s) allowed, "
                    f"found {i['technique']['method']} for {i['system_id']}."
                )

        assert not msg, "\n".join(msg)

    def validate_pdx_only(self, experiments):
        """Make sure experiments come from PDX samples."""
        msg = []

        for i in experiments:
            if not i.get("is_pdx"):
                msg.append(f"{i['system_id']} sample is not PDX derived")

        assert not msg, "\n".join(msg)

    def validate_dna_only(self, experiments):
        """Make sure experiments are DNA data."""
        msg = []

        for i in experiments:
            if i["technique"]["category"] != "DNA":
                msg.append(f"{i['system_id']} category is not DNA")

        assert not msg, "\n".join(msg)

    def validate_rna_only(self, experiments):
        """Make sure experiments are RNA data."""
        msg = []

        for i in experiments:
            if i["technique"]["category"] != "RNA":
                msg.append(f"{i['system_id']} category is not RNA")

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
        assert rtec == ttec, f"Same techniques required: {ttec}, {rtec}."

    def validate_species(self, experiments):
        """Validate experiments's species is same as application's setting."""
        msg = []

        for i in experiments:
            if i["sample"]["individual"]["species"] != self.SPECIES:
                msg.append(f'{i["system_id"]} species not supported')

        assert not msg, "\n".join(msg)
