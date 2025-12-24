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
import analytics
import click

from isabl_cli import api
from isabl_cli import data
from isabl_cli import exceptions
from isabl_cli import utils
from isabl_cli import options
from isabl_cli.batch_systems import submit_local
from isabl_cli.settings import get_application_settings
from isabl_cli.settings import system_settings


class AbstractApplication:  # pylint: disable=too-many-public-methods

    """An Abstract Isabl application."""

    # The uniqueness of an application is determined by it's name and version.
    # A good strategy to applications is to ask are results still comparable?
    # An optimization that doesn't change outputs might not require a version change.
    NAME = None
    VERSION = None

    # Optionally set ASSEMBLY and SPECIES to version as a function of genome build.
    # This is particularly useful for NGS applications as often results are only
    # comparable if data was analyzed against the same version of the genome.
    ASSEMBLY = None
    SPECIES = None

    # URL (or comma separated URLs) to be stored in the application database object.
    application_url = None
    application_description = ""

    # Applications can depend on multiple configurations such as paths to executables,
    # references files, compute requirements, etc. These settings are explicitly
    # defined using the application_settings dictionary. Learn more at:
    # https://docs.isabl.io/writing-applications#application-settings
    application_settings = dict()
    application_import_strings = set()

    # Applications can be launched from the command line. To support this capability
    # you have to tell the application how to link analyses to different experiments.
    # Learn more: https://docs.isabl.io/writing-applications#command-line-configuration
    cli_help = ""
    cli_options = []
    cli_allow_force = True
    cli_allow_restart = True
    cli_allow_local = True

    # You can provide an specification your application results using this attribute.
    # Each key is a result id and the value is a dictionary with specs of the result.
    # By default, analysis results are protected upon completion (i.e. permissions are
    # set to read only).
    application_results = dict()

    # Disabling application_protect_results makes an app re-runnable. Meaning that
    # analyses in SUCCEEDED status can be re-run and overwritten.
    application_protect_results = True

    # Isabl applications can produce auto-merge analyses at a project and individual
    # level. For example, you may want to merge variants whenever new results are
    # available for a given project, or update quality control reports when a new
    # sample is added to an individual. A newly versioned analysis will be created for
    # each type of auto-merge and your role is to take a list of succeeded analysis
    # and implement the merge logic.
    application_project_level_results = dict()
    application_individual_level_results = dict()

    # application_inputs are analysis-specific settings (settings are the same for all
    # analyses, yet inputs are potentially different for each analysis). Each input set
    # to NotImplemented is considered required and must be resolved at get_dependencies
    application_inputs = dict()

    dependencies_results = []

    # It is possible to create applications that are unique at the individual level.
    # A good example of a unique per individual application could be a patient centric
    # report that aggregates results across all samples. Applications that require a
    # unique analysis per individual don't support individual level auto-merge.
    unique_analysis_per_individual = False

    # Analyses in these status won't be prepared for submission. To re-rerun SUCCEEDED
    # analyses see unique_analysis_per_individual. To re-rerun failed analyses use
    # either --force or --restart.
    skip_status = {
        "FAILED",
        "FINISHED",
        "STARTED",
        "SUBMITTED",
        "SUCCEEDED",
        "REJECTED",
    }

    # If any of these errors is raised during the command generation process, the
    # submission will continue. Errors or valdation messages are presented at the end.
    skip_exceptions = (
        AssertionError,
        click.UsageError,
        exceptions.MissingRequirementError,
        exceptions.ConfigurationError,
        exceptions.MissingOutputError,
    )

    # private variables
    _staged_message = "READY FOR SUBMISSION"

    # private result keys
    _command_script_key = "command_script"
    _command_log_key = "command_log"
    _command_err_key = "command_err"
    _base_results = {
        _command_script_key: {
            "frontend_type": "ansi",
            "description": "Script used to execute the analysis.",
            "verbose_name": "Analysis Script",
        },
        _command_log_key: {
            "frontend_type": "ansi",
            "description": "Analysis standard output.",
            "verbose_name": "Standard Output",
        },
        _command_err_key: {
            "frontend_type": "ansi",
            "description": "Analysis standard error.",
            "verbose_name": "Standard Error",
        },
    }

    # -----------------------------
    # USER REQUIRED IMPLEMENTATIONS
    # -----------------------------

    def get_command(self, analysis, inputs, settings):  # pylint: disable=W9008
        """
        Must return a shell command for the analysis as a string.

        Arguments:
            analysis (dict): an analysis object as retrieved from API.
            inputs (dict): as returned by `get_dependencies`.
            settings (dict): applications settings.

        Returns:
            str: command
        """
        raise NotImplementedError()

    # ------------------------
    # OPTIONAL IMPLEMENTATIONS
    # ------------------------

    @abc.abstractmethod
    def get_experiments_from_cli_options(self, **cli_options):  # pylint: disable=W9008
        """
        Must return list of target-reference experiment tuples given the parsed options.

        Arguments:
            cli_options (dict): parsed command line options.

        Returns:
            list: of (targets, references) tuples.
        """
        raise NotImplementedError()

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
        """
        Validate settings.

        Make sure application settings are valid by raising an AssertionError if
        something is properly configured.

        Arguments:
            settings (isabl_cli.settings.ApplicationSettings): an application settings
            object, which is also a Munch-like dictionary.
        """
        return

    def validate_experiments(self, targets, references):  # pylint: disable=W9008
        """
        Raise AssertionError if tuple combination isn't valid.

        Some of the advantages of metadata-driven applications is that we can prevent
        analyses that don't make sense, for example running a variant calling
        application on imaging data.

        Arguments:
            targets (list): list of targets dictionaries.
            references (list): list of references dictionaries.

        Raises:
            AssertionError: if tuple is invalid.
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
            instance (obj): a project or individual instance.
        """
        submit_merge = system_settings.SUBMIT_MERGE_ANALYSIS

        if submit_merge:
            click.secho(
                f"Submitting merge analyses for {instance} "
                f"using: {submit_merge.__module__}.{submit_merge.__name__}"
            )

            submit_merge(
                instance=instance,
                application=self,
                command=self._get_cli_merge_command(instance),
            )
        elif "species" in instance:
            self.run_individual_merge(instance)
        else:
            self.run_project_merge(instance)

    def _run_analyses_merge(self, instance, analyses):
        merge_type = "project"
        merge_analyses = self.merge_project_analyses
        validate_analyses = self.validate_project_analyses
        get_analysis = self.get_project_level_auto_merge_analysis

        if "species" in instance:
            merge_type = "individual"
            merge_analyses = self.merge_individual_analyses
            validate_analyses = self.validate_individual_analyses
            get_analysis = self.get_individual_level_auto_merge_analysis

        if not analyses or len(analyses) < 2:
            click.secho(
                f"Not enough analyses for {instance} merge, "
                f"at least 2 required but got: {len(analyses)}",
                err=True,
                fg="yellow",
            )
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
        error_msg = ""

        with open(stdout_path, "w") as out, open(stderr_path, "w") as err:
            with redirect_stdout(out), redirect_stderr(err):
                try:
                    # TODO: setting to submitted here temporarily, will fix later
                    api.patch_analysis_status(analysis, "SUBMITTED")
                    api.patch_analysis_status(analysis, "STARTED")
                    merge_analyses(analysis, analyses)
                    api.patch_analysis_status(analysis, "SUCCEEDED")
                    analysis.data["merged_analyses"] = len(analyses)
                except Exception as e:  # pragma: no cover pylint: disable=W0703
                    error_msg = traceback.format_exc()
                    click.echo(error_msg, file=sys.stderr)
                    click.echo(e, file=sys.stderr)
                    api.patch_analysis_status(analysis, "FAILED")
                    error = e
                    print(error_msg)

        api.patch_instance("analyses", analysis.pk, data=analysis.data)
        os.umask(oldmask)
        if error is not None:
            raise Exception(
                f"Merge analysis ({analysis}) failed for {instance}. You can "
                f"retry this operation with `isabl merge-{merge_type}-analyses "
                f"--{merge_type} {instance.pk} "
                f"--application {analyses[0].application.pk}`. "
                f"The error traceback, if any, was: {error} {error_msg}"
            )

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
            "don't support individual level auto-merge."
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
            "don't support individual level auto-merge."
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
        return str(system_settings.client.get("pk", "default_client"))

    @cached_property
    def primary_key(self):
        """Space separated name, version and assembly."""
        return self.application["pk"]

    @cached_property
    def settings(self):
        """Return the application settings."""
        reference_data = dict()
        import_strings = set(self.application_import_strings)
        import_strings.add("submit_analyses")
        defaults = self.application_settings.copy()
        defaults["submit_analyses"] = defaults.get("submit_analyses", None)

        if self.application.assembly:
            reference_data = self.application.assembly.reference_data or {}

        return get_application_settings(
            defaults=defaults,
            settings=self._settings_for_client,
            reference_data=reference_data,
            import_strings=import_strings,
        )

    @cached_property
    def application(self):
        """Get application database object."""
        assert all([self.NAME, self.VERSION]), (
            f"NAME must be set: {self.NAME}\n" f"VERSION must be set: {self.VERSION}\n"
        )

        if any([self.ASSEMBLY, self.SPECIES]):
            assert all([self.ASSEMBLY, self.SPECIES]), (
                f"ASSEMBLY must be set: {self.ASSEMBLY}\n"
                f"SPECIES must be set: {self.SPECIES}\n"
            )

        application = api.create_instance(
            endpoint="applications",
            name=self.NAME,
            version=self.VERSION,
            assembly=(
                {"name": self.ASSEMBLY, "species": self.SPECIES}
                if self.ASSEMBLY
                else None
            ),
        )

        if application.settings.get("default_client") is None:
            application.settings["default_client"] = {}

        return api.patch_instance(
            description=self.application_description,
            endpoint="applications",
            instance_id=application["pk"],
            application_class=f"{self.__module__}.{self.__class__.__name__}",
            results=self._application_results,
            url=self.application_url,
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
            url=self.application_url,
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
            url=self.application_url,
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

        local = click.option(
            "--local",
            help="Run analyses' head jobs locally, one after the other one.",
            is_flag=True,
        )

        cli_options = [  # pylint: disable=unused-variable
            (quiet, True),
            (commit, True),
            (force, cls.cli_allow_force),
            (restart, cls.cli_allow_restart),
            (local, cls.cli_allow_local),
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
            local = cli_options.pop("local", False)
            tuples = []

            if commit and force:
                raise click.UsageError("--commit not required when using --force")

            if commit and restart:  # pragma: no cover
                raise click.UsageError("--restart not required when using --force")

            if force and restart:
                raise click.UsageError("cant use --force and --restart together")

            if not hasattr(
                cls.get_experiments_from_cli_options, "__isabstractmethod__"
            ):
                tuples = pipe.get_experiments_from_cli_options(**cli_options)
            else:
                tuples = cls.get_experiments_from_default_cli_options(cli_options)

            pipe.run(
                tuples=tuples,
                commit=commit,
                force=force,
                verbose=not quiet,
                restart=restart,
                local=local,
                run_args=cli_options,
            )

            if not (commit or force or restart):
                utils.echo_add_commit_message()

        return command

    def get_cli_command_name(self):
        """Get name for isabl_cli command."""
        return f"{slugify(self.NAME)}-{slugify(self.VERSION, separator='.')}"

    @classmethod
    def get_experiments_from_default_cli_options(cls, cli_options):
        """Get experiments from default CLI options."""
        tuples = []
        references = []
        supported_cli_options = [
            options.ANALYSES,
            options.NORMAL_TARGETS,
            options.TARGETS,
            options.REFERENCES,
            options.NULLABLE_REFERENCES,
            options.PAIR,
            options.PAIRS,
            options.PAIRS_FROM_FILE,
        ]

        assert any(i in cls.cli_options for i in supported_cli_options), (
            f"'{cls.__name__}.cli_options' must include at least one of the "
            f"supported default options to get experiments: {supported_cli_options}"
        )

        if options.PAIR in cls.cli_options:
            tuples += cli_options.get("pair", [])

        if options.PAIRS in cls.cli_options:
            tuples += cli_options.get("pairs", [])

        if options.PAIRS_FROM_FILE in cls.cli_options:
            tuples += cli_options.get("pairs_from_file", [])

        if options.ANALYSES in cls.cli_options:
            for i in cli_options.get("analyses_filters", []):
                tuples.append((i.targets, i.references))

        if (
            options.REFERENCES in cls.cli_options
            or options.NULLABLE_REFERENCES in cls.cli_options
        ):
            references += cli_options.get("references", [])

        if any(i in cls.cli_options for i in [options.NORMAL_TARGETS, options.TARGETS]):
            for i in cli_options.get("targets", []):
                tuples.append(([i], references))

        return tuples

    # ------------------------
    # ANALYSES EXECUTION LOGIC
    # ------------------------

    def run(
        self,
        tuples,
        commit,
        force=False,
        restart=False,
        verbose=True,
        run_args=None,
        local=False,
    ):
        """
        Run a list of targets, references tuples.

        Arguments:
            restart (bool): set settings.restart = True.
            force (bool): if true, analyses are wiped before being submitted.
            local (bool): if true, analyses will be run locally one by one.
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
        self.settings.run_args = api.isablfy(run_args or {})

        # run extra settings validation
        self.validate_settings(self.settings)

        # create analyses
        analyses, invalid_tuples = self.get_or_create_analyses(tuples)

        # make sure outdir is set
        for i in analyses:
            if not i.storage_url:  # pragma: no cover
                self._patch_analysis(i)

        # run analyses
        run_tuples, skipped_tuples, invalid_run_tuples = self.run_analyses(
            analyses=analyses, commit=commit, force=force, restart=restart, local=local
        )

        invalid_tuples.extend(invalid_run_tuples)

        if verbose:
            self.echo_run_summary(run_tuples, skipped_tuples, invalid_tuples)

        if commit:
            click.echo(
                f"RAN {len(run_tuples)} | "
                f"SKIPPED {len(skipped_tuples)} | "
                f"INVALID {len(invalid_tuples)}\n"
            )

        else:
            click.echo(
                f"STAGED {len(run_tuples)} | "
                f"SKIPPED {len(skipped_tuples)} | "
                f"INVALID {len(invalid_tuples)}\n"
            )

            num_run_on_commit = len(run_tuples)
            num_succeeded = 0
            num_failed = 0

            for i in skipped_tuples:
                if i[1] == "SUCCEEDED":
                    num_succeeded += 1
                    if not self.application_protect_results:
                        num_run_on_commit += 1
                else:
                    num_failed += 1

            if not self.application_protect_results:
                if num_run_on_commit == 1:
                    click.echo(f"{num_run_on_commit} analysis available to run:")
                else:
                    click.echo(f"{num_run_on_commit} analyses available to run:")

                click.echo(f"\t{len(run_tuples)} STAGED")
                click.echo(f"\t{num_succeeded} SUCCEEDED (Unprotected)")

        return run_tuples, skipped_tuples, invalid_tuples

    def run_analyses(self, analyses, commit, force, restart, local):
        """
        Run a list of analyses.

        Returns:
            list: tuple of
        """
        skipped_tuples = []
        invalid_tuples = []
        command_tuples = []
        submit_analyses = (
            submit_local
            if local
            else (self.settings.submit_analyses or system_settings.SUBMIT_ANALYSES)
        )

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
                elif force and i["status"] not in {"SUCCEEDED", "FINISHED", "REJECTED"}:
                    system_settings.TRASH_ANALYSIS_STORAGE(i)
                    utils.makedirs(i["storage_url"])

                # only restart failed analyses
                elif (
                    not (restart and i["status"] == "FAILED")
                    and i["status"] in self.skip_status
                ):
                    skipped_tuples.append((i, i["status"]))
                    continue

                elif restart and i.ran_by != system_settings.api_username:
                    invalid_tuples.append(
                        (
                            i,
                            "Can't restart: started by different user. Consider --force",
                        )
                    )
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

        # track user app usage from cli
        self.send_analytics(
            command_tuples=command_tuples,
            skipped_tuples=skipped_tuples,
            invalid_tuples=invalid_tuples,
            commit=commit,
            force=force,
            restart=restart,
            submitter=submit_analyses.__name__,
        )

        if commit:
            click.echo(f"Running analyses with {submit_analyses.__name__}...")
            run_tuples = submit_analyses(self, command_tuples)
        else:
            run_tuples = [(i, self._staged_message) for i, _ in command_tuples]

        return run_tuples, skipped_tuples, invalid_tuples

    def send_analytics(
        self,
        command_tuples,
        skipped_tuples,
        invalid_tuples,
        commit,
        force,
        restart,
        submitter,
    ):
        """Send analytics event of analyses ran from cli."""
        analyses_tuples = []
        for i, _ in command_tuples:
            status = (
                "FORCED"
                if force
                else (
                    "RESTARTED"
                    if restart
                    else "SUBMITTED" if commit else self._staged_message
                )
            )
            analyses_tuples.append((i, status))
        for i, _ in skipped_tuples + invalid_tuples:
            analyses_tuples.append((i, "INVALID"))

        analyses = []
        for i, status in analyses_tuples:
            analysis = {
                "analysis": i.pk,
                "application": {
                    "name": i.application.name,
                    "pk": i.application.pk,
                    "version": i.application.version,
                },
                "status": status,
            }
            analyses.append(analysis)

        analytics.track(
            system_settings.api_username,
            "Ran app",
            {
                "analyses": analyses[:200],
                "total": len(analyses),
                "submitter": submitter,
                "valid": len(command_tuples),
                "invalid": len(skipped_tuples) + len(invalid_tuples),
            },
        )

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
        utils.makedirs(outdir)
        status = self._get_after_completion_status(analysis)

        # build and write command
        tmpdir = os.getenv("TMP", "/tmp")
        tmpdir = " && ".join(
            f"export {i}={tmpdir}" for i in ["TMP", "TMPDIR", "TMP_DIR"]
        )
        failed = self.get_patch_status_command(analysis["pk"], "FAILED")
        started = self.get_patch_status_command(analysis["pk"], "STARTED")
        finished = self.get_patch_status_command(analysis["pk"], status)
        command = (
            f"umask g+wrx && date && cd {outdir} && {tmpdir} && "
            f"{started} && {command} && {finished} && date"
        )

        with open(self.get_command_script_path(analysis), "w") as f:
            f.write(f"{{\n\n    {command}\n\n}} || {{\n\n    {failed}\n\n}}")

    @staticmethod
    def get_patch_status_command(key, status):
        """Return a command to patch the `status` of a given analysis `key`."""
        return f"isabl patch-status --key {key} --status {status}"

    def _get_after_completion_status(self, analysis):
        status = self.get_after_completion_status(analysis)
        assert status in {"IN_PROGRESS", "FINISHED"}, "Status not supported"

        if system_settings.is_admin_user and status == "FINISHED":
            status = "SUCCEEDED"

        return status

    def _get_dependencies(self, targets, references):
        missing = []
        if self.dependencies_results:
            analyses, inputs = self._get_dependencies_results(targets, references)
        else:
            analyses, inputs = self.get_dependencies(targets, references, self.settings)

        for i, j in self.application_inputs.items():  # pragma: no cover
            if i not in inputs:
                if j is NotImplemented:
                    missing.append(i)
                else:
                    inputs[i] = j

        if missing:  # pragma: no cover
            missing = ", ".join(map(str, missing))
            raise exceptions.ConfigurationError(
                f"Required inputs missing from `get_dependencies`: {missing}"
            )

        return analyses, inputs

    def _get_dependencies_results(self, targets, references):
        """
        Get dependencies' results from a defined application version.

        It's called when `self.dependencies_results` is an array containing an object:
            {
                result (str): Result key `Application.results.result_key`.
                name (str): Name the result will have in the inputs objects.
                app (obj): `Application` instance.
                app_name (str): `Application.name`.
                app_version (str): `Application.version`. If not defined, use
                    any available. Use `any` if you want to use the latest available.
                app_assembly (str): `Application.assembly.name`. i.e: GRCh37
                linked (bool): if False, the analysis is not linked as a dependencie of
                    the analysis. As adding new dependencies to an analysis forces isabl
                    to not recognize existing ones (Default: True).
            }
        """
        inputs = {}
        analyses = []
        for dependency in self.dependencies_results:
            result_args = {
                "result_key": dependency["result"],
                "targets": targets,
                "references": references,
            }
            # Match app by name, and optionally by version
            if dependency.get("app_name"):
                result_args["application_name"] = dependency.get("app_name")
                if "app_version" in dependency:
                    result_args["application_version"] = dependency.get("app_version")
            else:
                # Match app by primary key
                result_args["application_key"] = dependency.get("app").primary_key
                result_args["application_name"] = dependency.get("app").NAME

            if "app_assembly" in dependency:
                result_args["application_assembly"] = dependency.get("app_assembly")

            if "status" in dependency:
                result_args["status"] = dependency["status"]

            input_name = dependency.get("name")
            inputs[input_name], key = self.get_result(targets[0], **result_args)

            if not "linked" in dependency or dependency["linked"]:
                # Avoid linking the analysis as dependency to avoid creating new runs.
                analyses.append(key)

        return analyses, inputs

    def _get_analysis_results_from_patterns(self, analysis, specification):
        """Get first matching file, if "pattern" is in result specification."""
        results = {}
        for name, result_specification in specification.items():
            pattern = result_specification.get("pattern")
            exclude = result_specification.get("exclude")
            optional = result_specification.get("optional")
            if pattern:
                results[name] = utils.first_matching_file(
                    directory=analysis.storage_url,
                    pattern=pattern,
                    exclude=exclude,
                    optional=optional,
                )
        return results

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

        results.update(
            self._get_analysis_results_from_patterns(analysis, specification)
        )
        results.update(get_results(analysis))

        for i in specification:
            assert i in results, f"Missing expected result {i} in: {results}"

        # update bam with result if specified
        for key, attributes in specification.items():
            if attributes.get("store_as_bam") and results[key]:
                self.update_experiment_bam_file(
                    experiment=analysis["targets"][0],
                    bam_url=results[key],
                    analysis_pk=analysis["pk"],
                )

        return results

    # -----------------
    # APPLICATION UTILS
    # -----------------

    @staticmethod
    def get_result(*args, **kwargs):  # pragma: no cover
        """Get an application result."""
        return utils.get_result(*args, **kwargs)

    @staticmethod
    def get_results(*args, **kwargs):
        """Get application results."""
        return utils.get_results(*args, **kwargs)

    def patch_application_settings(self, client_id=None, **settings):
        """Patch application settings if necessary."""
        assert system_settings.is_admin_user, "Apps can be patched only by admin user."
        click.echo(f"Patching settings: {self.NAME} {self.VERSION} {self.ASSEMBLY}\n")
        client_id = str(client_id or self.client_id)  # must be a string

        try:
            assert self.application.settings.get(client_id) == settings
            click.secho(
                f"\tNo changes detected, skipping patch.\n", err=True, fg="yellow"
            )
        except AssertionError:
            try:
                api.patch_instance(
                    "applications",
                    self.primary_key,
                    settings={**self.application.settings, client_id: settings},
                )

                try:
                    del self.settings  # make sure cached settings are re-computed
                except AttributeError:
                    pass

                try:
                    del self.application  # make sure application is refetched
                except AttributeError:
                    pass

                click.secho("\tSuccessfully patched settings.\n", fg="green")
            except TypeError as error:  # pragma: no cover
                click.secho(
                    f"\tPatched failed with error: {error}.\n", err=True, fg="red"
                )

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
        else:  # pragma: no cover
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

            if isinstance(msg, Exception):  # pragma: no cover
                color = "red"
            elif msg == "SUCCEEDED":
                color = "green"
            elif msg == "FAILED":
                color = "red"
            elif msg == "INVALID":
                color = "yellow"
            elif msg == "SUBMITTED":  # pragma: no cover
                color = "cyan"
            elif msg == self._staged_message:  # pragma: no cover
                color = "magenta"
                blink = True

            if len(msg) > 20:  # pragma: no cover
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

    # ------------------
    # NGS SPECIFIC UTILS
    # ------------------

    def get_bedfile(self, experiment, bedfile_type="targets"):
        """Get targets or baits bedfile for experiment."""
        bedfile_key = f"{self.ASSEMBLY}_{bedfile_type}_bedfile"
        assert bedfile_type in {"targets", "baits"}, "Unsupported bedfile type."
        return experiment["technique"]["reference_data"][bedfile_key]["url"]

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
            ) as error:  # pragma: no cover
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

        # validate existing analyses
        with click.progressbar(
            existing_analyses,
            label=f"Validating tuples of the {len(existing_analyses)} existing analyses...\t\t",
        ) as bar:
            for i in bar:
                try:
                    self.validate_experiments(i["targets"], i["references"])
                    created_analyses.append(i)
                except (
                    exceptions.ValidationError,
                    AssertionError,
                ) as error:  # pragma: no cover
                    invalid_tuples.append((i, exceptions.ValidationError(*error.args)))

        # create new analyses and validate
        with click.progressbar(
            valid_tuples,
            label=f"Creating analyses for {len(valid_tuples)} tuples...\t\t",
        ) as bar:
            for i in bar:
                try:
                    targets, references, analyses, inputs, individual = i
                    # self.validate_species(targets + references)
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

        return created_analyses, invalid_tuples

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

        click.echo("Checking for existing analyses...")
        projects = {j["pk"] for i in tuples for j in i[0][0]["projects"]}
        cache = defaultdict(list)
        existing, missing = [], []
        targets_pks = ",".join(str(j.pk) for i in tuples for j in i[0])
        filters = dict(application=self.application["pk"])

        if projects:
            filters["projects__pk__in"] = projects

        if targets_pks:
            filters["targets__pk__in"] = targets_pks

        def get_cache_key(targets, references):
            return (
                tuple(sorted(set(i["pk"] for i in targets))),
                tuple(sorted(set(i["pk"] for i in references))),
            )

        for i in api.get_instances("analyses", limit=5000, **filters):
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

        if tuples_map:
            for i in api.get_analyses(
                application=self.application.pk,
                individual_level_analysis__pk__in=",".join(map(str, tuples_map)),
            ):
                individual = i.individual_level_analysis
                current_tuple = tuples_map[individual.pk]

                # make sure we only have one analysis per individual
                assert (
                    individual.pk not in existing
                ), f"Multiple analyses for {individual}"
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
        return f"{self.NAME} {self.VERSION}" + (
            f" ({self.ASSEMBLY})" if self.ASSEMBLY else ""
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

    def validate_same_platform(self, targets, references):
        """Validate targets and references are sequenced on the same platform."""
        tpla = {i["platform"]["slug"] for i in targets}
        rpla = {i["platform"]["slug"] for i in references}
        assert len(rpla) == 1, f"Expected one platform, got: {rpla}"
        assert len(tpla) == 1, f"Expected one platform, got: {tpla}"
        assert rpla == tpla, f"Same platforms required: {tpla}, {rpla}."

    def validate_species(self, experiments):
        """Validate experiments's species is same as application's setting."""
        msg = []

        for i in experiments:
            if self.ASSEMBLY and i["sample"]["individual"]["species"] != self.SPECIES:
                msg.append(f'{i["system_id"]} species not supported')

        assert not msg, "\n".join(msg)

    def validate_are_normals(self, experiments):
        """Raise error not all experiments come from NORMAL sample."""
        for i in experiments:
            msg = f"Experiment Sample {i.sample.system_id} is not NORMAL."
            assert i.sample.category == "NORMAL", msg

    def validate_individuals(self, targets, references):
        """
        Validate the individual of the targets and references.

        Validate pairs are of the same individual if the pipeline is matched.
        Validate pairs are of the different individuals if the pipeline is unmatched.
        """
        if references:
            targets_set = list(
                {v["sample"]["individual"]["pk"]: v for v in targets}.values()
            )
            references_set = list(
                {v["sample"]["individual"]["pk"]: v for v in references}.values()
            )

            assert len(targets_set) == 1, "One unique target individual is supported."
            assert (
                len(references_set) == 1
            ), "One unique reference individual is supported."

            tind = targets_set[0]["sample"]["individual"]
            rind = references_set[0]["sample"]["individual"]

            if (
                hasattr(self, "IS_UNMATCHED") and self.IS_UNMATCHED
            ):  # pylint: disable=no-member
                assert tind["pk"] != rind["pk"], (
                    "Different individuals required: "
                    f"{tind['system_id']} and {rind['system_id']} "
                    "are of the same individual."
                )
            else:  # application is designed for matched analyses
                assert tind["pk"] == rind["pk"], (
                    "Same individual required: "
                    f"{tind['system_id']} and {rind['system_id']} "
                    "are of different individuals."
                )

    def validate_source(self, experiments, source):
        """Validate experiments are from a specific source material."""
        for i in experiments:
            assert (
                i["sample"]["source"] == source
            ), f"Sample source for {i['sample']['system_id']} does not match {source}."

    # -------------------------
    # NOTIFICATION UTILS
    # -------------------------

    def notify_project_analyst(self, analysis, subject, message):
        """Send email notification to project's analysts of projects."""
        projects = []
        for target in analysis.targets:
            projects.extend(target.projects)
        analysts = set([project.analyst for project in projects if project.analyst])
        if not analysts:
            click.secho(
                "Skipping notification as projects have no registered analysts",
                err=True,
                fg="red",
            )
            return
        kwargs = {
            "data": {"recipients": analysts, "subject": subject, "content": message}
        }
        return api.api_request("post", url=f"/send_email", **kwargs)
