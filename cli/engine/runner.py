"""Engine runner."""

from itertools import chain
from os.path import join
import os
import random
import subprocess

from click import progressbar
import click

from cli import exceptions
from cli import system_settings
from cli import utils

from .creator import Creator


class Runner(Creator):

    """
    Analyses running logic.

    Implementors must define logic to build the command and get requirements.

    Attributes:

        _skip_status (set): set of analyses status to skip.
        _skip_exceptions (bool): expected exceptions that may occur.
    """
    _skip_status = {'FAILED', 'FINISHED', 'STARTED', 'SUBMITTED', 'SUCCEEDED'}
    _skip_exceptions = (
        click.UsageError,
        exceptions.MissingRequirementError,
        exceptions.ConfigurationError,
        exceptions.MissingOutputError)

    @staticmethod
    def get_patch_status_command(key, status):
        """Return a command to patch the `status` of a given analysis `key`."""
        return f'cli patch_status --key {key} --status {status}'

    @staticmethod
    def get_command(analysis):
        """Must return string with command to be run."""
        raise NotImplementedError

    @staticmethod
    def get_status(analysis):  # pylint: disable=W0613
        """Possible values are FINISHED and IN_PROGRESS."""
        return 'FINISHED'

    @staticmethod
    def get_job_name(analysis):
        """Get job name given an analysis dict."""
        targets = analysis['targets']
        references = analysis['references']
        methods = ' '.join({i['technique']['method'] for i in targets})
        projects = ' '.join({str(j) for i in targets for j in i['projects']})

        if len(targets) > 2 or not targets:
            targets = f'{len(targets)} samples.'
        else:
            targets = ' '.join([i['system_id'] for i in targets])

        if len(references) > 2 or not references:
            references = f'{len(references)} samples.'
        else:
            references = ' '.join([i['system_id'] for i in references])

        return (
            f'targets: {targets} | references: {references} | '
            f'methods: {methods} | analysis: {analysis["pk"]} | '
            f'projects: {projects} | rundir: {analysis["storage_url"]} | '
            f'pipeline: {analysis["pipeline"]["pk"]}')

    @staticmethod
    def submit_analyses(command_tuples):
        """
        Submit pipelines as arrays grouped by the target methods.

        Arguments:
            command_tuples (list): of (analysis, command) tuples.
        """
        label = f'Running {len(command_tuples)} analyses...'
        with progressbar(command_tuples, label=label) as bar:
            for _, i in bar:
                try:
                    subprocess.check_call(['bash', i])
                except subprocess.CalledProcessError:
                    pass

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

        if not tuples:
            return None

        analyses, invalid_tuples = self.get_or_create_analyses(tuples)
        command_tuples, skipped_tuples = [], []

        with progressbar(analyses, label='Building commands...\t\t') as bar:
            for i in bar:
                if force and i['status'] not in {'SUCCEEDED', 'FINISHED'}:
                    system_settings.TRASH_ANALYSIS_STORAGE_FUNCTION(i)
                    os.makedirs(i['storage_url'], exist_ok=True)
                elif i['status'] in self._skip_status:
                    skipped_tuples.append((i, i['status']))
                    continue

                try:
                    command = self._build_command_script(i)
                    command_tuples.append((i, command))
                except self._skip_exceptions as error:
                    skipped_tuples.append((i, error))

        if verbose:
            self._echo_summary(command_tuples, skipped_tuples, invalid_tuples)

        click.echo(
            f"RUNNING {len(command_tuples)} | "
            f"SKIPPED {len(skipped_tuples)} | "
            f"INVALID {len(invalid_tuples)}")

        if commit:
            self.submit_analyses(command_tuples)

        return command_tuples, skipped_tuples, invalid_tuples

    def _build_command_script(self, analysis):
        """
        Get analysis command as a path to a bash script.

        Arguments:
            analysis (leuktools.Analysis): the analysis object to be run.

        Returns:
            str: analysis command as a path to a bash script.
        """
        command = self.get_command(analysis)
        status = self.get_status(analysis)
        command_path = join(analysis['storage_url'], "head_job.sh")
        assert status in {'IN_PROGRESS', 'FINISHED'}

        # analyses not run but admin will be marked as succeeded after
        if not system_settings.is_admin_user and status == 'FINISHED':
            status = 'FINISHED'

        command = (
            f"sleep {random.uniform(0, 10):.3} && "  # avoid parallel API hits
            f"cd {analysis['storage_url']} && "
            f"{self.get_patch_status_command(analysis['pk'], 'STARTED')}"
            f"&& date && {command} && "
            f"{self.get_patch_status_command(analysis['pk'], status)}")

        with open(command_path, "w") as f:
            f.write(command)

        return command_path

    def _echo_summary(self, command_tuples, skipped_tuples, invalid_tuples):
        """
        Echo errors for error tuples such as `invalid`, `cant_run`.

        Arguments:
            invalid_tuples (list): of (tuple, error) tuples.
            command_tuples (list): of (analysis, command) tuples.
            skipped_tuples (list): of (analysis, skip reason) tuples.
        """
        summary_keys = ["PK", "PROJECTS", "TARGETS", "REFERENCES", "MESSAGE"]
        summary = ["\n"]
        command_tuples.sort(key=lambda i: i[1])
        invalid_tuples.sort(key=lambda i: i[1])
        skipped_tuples.sort(key=lambda i: i[1])
        tuples = chain(command_tuples, invalid_tuples, skipped_tuples)

        for i, msg in tuples:
            row = {k: "NA" for k in summary_keys}
            row["PK"] = getattr(i, "pk", "NA")
            row["MESSAGE"] = self._style_msg(msg)
            row["PROJECTS"] = self._style_projects(i['targets'])
            row["TARGETS"] = self._style_workflows(i['targets'])
            row["REFERENCES"] = self._style_workflows(i['references'])
            summary.append("\t".join(row[k] for k in summary_keys))

        summary = "\n".join(summary).expandtabs(12) + "\n"
        click.echo(f"{summary}\n")

    @staticmethod
    def _style_msg(msg):
        """Color message depending on its nature."""
        if isinstance(msg, exceptions.MissingRequirementError):
            color = "yellow"
        elif isinstance(msg, Exception):
            color = "red"
        elif msg == "SUCCEEDED":
            color = "green"
        else:
            color = "blue"
        return click.style(msg, fg=color)

    @staticmethod
    def _style_projects(targets):
        """Color message depending on its nature."""
        keys = set()

        for i in targets:
            for j in i['projects']:
                keys.add(f"{j} {i['technique']['method']}")

        return f"Projects {' '.join(keys)}"

    @staticmethod
    def _style_workflows(workflows):
        """Color message depending on its nature."""
        if len(workflows) > 2 or not workflows:
            message = f"{len(workflows)} workflows."
        else:
            message = " ".join([i['system_id'] for i in workflows])

        return message
