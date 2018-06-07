"""Engine runner."""

from os.path import join
import os
import subprocess

from click import progressbar
import click

from cli import exceptions
from cli import system_settings
from cli import utils
from cli import api

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
    def get_command(analysis, settings):
        """Must return string with command to be run."""
        raise NotImplementedError

    @staticmethod
    def get_status(analysis):  # pylint: disable=W0613
        """Possible values are FINISHED and IN_PROGRESS."""
        raise NotImplementedError

    @staticmethod
    def get_patch_status_command(key, status):
        """Return a command to patch the `status` of a given analysis `key`."""
        return f'cli patch_status --key {key} --status {status}'

    @staticmethod
    def submit_analyses(command_tuples):
        """
        Submit pipelines as arrays grouped by the target methods.

        Arguments:
            command_tuples (list): of (analysis, command) tuples.
        """
        label = f'Running {len(command_tuples)} analyses...'
        with progressbar(command_tuples, label=label) as bar:
            for i, j in bar:
                try:
                    log = join(i['storage_url'], 'head_job.log')
                    err = join(i['storage_url'], 'head_job.err')
                    with open(log, 'w') as stdout, open(err, 'w') as stderr:
                        subprocess.check_call(
                            args=['bash', j],
                            stdout=stdout,
                            stderr=stderr)

                except subprocess.CalledProcessError:
                    api.patch_instance(
                        endpoint='analyses',
                        identifier=i['pk'],
                        storage_usage=utils.get_tree_size(i['storage_url']),
                        status='FAILED')

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

        analyses, invalid_tuples = self.get_or_create_analyses(tuples)
        command_tuples, skipped_tuples = [], []

        with progressbar(analyses, label='Building commands...\t\t') as bar:
            for i in bar:
                if force and i['status'] not in {'SUCCEEDED', 'FINISHED'}:
                    system_settings.TRASH_ANALYSIS_STORAGE(i)
                    os.makedirs(i['storage_url'], exist_ok=True)
                elif i['status'] in self._skip_status:
                    skipped_tuples.append((i, i['status']))
                    continue

                try:
                    command = self._build_command_script(i)
                    command_tuples.append((i, command))
                except self._skip_exceptions as error:  # pragma: no cover
                    skipped_tuples.append((i, error))

        if verbose:
            self._echo_summary(command_tuples, skipped_tuples, invalid_tuples)

        click.echo(
            f"RUNNING {len(command_tuples)} | "
            f"SKIPPED {len(skipped_tuples)} | "
            f"INVALID {len(invalid_tuples)}\n")

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
        command = self.get_command(analysis, system_settings.PIPELINES_SETTINGS)
        status = self.get_status(analysis)
        command_path = join(analysis['storage_url'], 'head_job.sh')
        assert status in {'IN_PROGRESS', 'FINISHED'}

        # analyses not run but admin will be marked as succeeded after
        if system_settings.is_admin_user and status == 'FINISHED':
            status = 'SUCCEEDED'

        command = (
            f'cd {analysis["storage_url"]} && umask g+wrx && '
            f'{self.get_patch_status_command(analysis["pk"], "STARTED")} '
            f'&& date && {command} && '
            f'{self.get_patch_status_command(analysis["pk"], status)}')

        if system_settings.is_admin_user:
            command += f" && chmod -R a-w {analysis['storage_url']}"

        with open(command_path, "w") as f:
            f.write(command)

        return command_path

    @staticmethod
    def _echo_summary(command_tuples, skipped_tuples, invalid_tuples):
        """
        Echo errors for error tuples such as `invalid`, `cant_run`.

        Arguments:
            invalid_tuples (list): of (tuple, error) tuples.
            command_tuples (list): of (analysis, command) tuples.
            skipped_tuples (list): of (analysis, skip reason) tuples.
        """
        summary_keys = ['PK', 'PROJECTS', 'TARGETS', 'REFERENCES', 'MESSAGE']
        summary = ['\n']
        command_tuples.sort(key=lambda i: i[1])
        invalid_tuples.sort(key=lambda i: i[1])
        skipped_tuples.sort(key=lambda i: i[1])

        def _style_msg(msg):
            color = 'blue'
            if isinstance(msg, Exception):
                color = 'red'
            elif msg == 'SUCCEEDED':
                color = 'green'
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
