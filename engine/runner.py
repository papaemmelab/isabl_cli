"""Engine runner."""

from collections import defaultdict
from datetime import datetime
from itertools import chain
from os.path import abspath
from os.path import dirname
from os.path import join
from getpass import getpass
import os
import random
import re
import subprocess
import sys

from click import progressbar
import click

from cli import system_settings
from cli import utils
from cli import exceptions
from .creator import Creator


class Runner(Creator):

    """
    Analyses running logic.

    Implementors must define logic to build the command and get requirements.

    Attributes:
        _force (bool): if true, analyses are wiped before being submitted.
        _commit (bool): if true, analyses are started (`_force` overwrites).
        _quiet (bool): whether or not verbose output should be printed.
        _skip_status (set): set of analyses status to skip.
        _skip_exceptions (bool): expected exceptions that may occur.
    """

    _force = False
    _commit = False
    _quiet = True
    _skip_status = {'FAILED', 'FINISHED', 'STARTED', 'SUBMITTED', 'SUCCEEDED'}
    _skip_exceptions = (
        click.UsageError,
        exceptions.MissingRequirementError,
        exceptions.ConfigurationError,
        exceptions.MissingOutputError)

    @staticmethod
    def build_command(analysis):
        """Must return string, see pipeline.py module for documentation."""
        raise NotImplementedError

    @staticmethod
    def get_requirements(sequencing_methods):
        """Must return string, See pipeline.py module for documentation."""
        raise NotImplementedError

    def _run_tuples(self, tuples):
        """
        Run a list of tuples.

        Arguments:
            tuples (list): list of (targets, references, analyses) tuples.
                elements can be objects or identifiers.
        """
        tuples = list(tuples)
        utils.echo_title(f'Running {len(tuples)} tuples for {self.pipeline}')

        if not tuples:
            return None

        analyses, invalid_tuples = self.get_or_create_analyses(tuples)
        command_tuples, skipped_tuples = [], []

        with progressbar(analyses, label='Building commands...\t\t') as bar:
            for i in bar:
                if self._force and i['status'] not in {'SUCCEEDED', 'FINISHED'}:
                    system_settings.TRASH_ANALYSIS_STORAGE_FUNCTION(i)
                    os.makedirs(i['storage_url'])
                elif i['status'] in self._skip_status:
                    skipped_tuples.append((i, i['status']))
                    continue

                try:
                    command = self._build_command(i)
                    command_tuples.append((i, command))
                except self._skip_exceptions as error:
                    skipped_tuples.append((i, error))

        self._echo_summary(invalid_tuples, command_tuples, skipped_tuples)

        if self._force or self._commit:
            self._submit_analyses(command_tuples)
        else:
            click.secho('\nAdd --commit to submit.\n', fg='green', blink=True)

    def _build_command(self, analysis):
        """
        Get analysis command as a path to a bash script.

        Arguments:
            analysis (leuktools.Analysis): the analysis object to be run.

        Returns:
            str: analysis command as a path to a bash script.
        """
        command, status, env = self.build_command(analysis)  # pylint: disable=E1111
        job_path = join(analysis['storage_url'], "head_job.sh")
        env_path = join(analysis['storage_url'], "head_job.env")
        env['RUNDIR'] = analysis['storage_url']
        env['JOBNAME'] = self.get_jobname(analysis)
        status = 'FINISHED' if status in {'FINISHED', 'SUCCEEDED'} else status

        # add a random sleep so many parallel jobs don't hit API at same time
        command = (
            f"sleep {random.uniform(0, 10):.3} && "
            f"cd {analysis.rundir} && source {env_path} && "
            f"{analysis.get_patch_command('STARTED')} && date && "
            f"{command} && {analysis.get_patch_command(status)}")

        with open(job_path, "w") as f:
            f.write(command)

        with open(env_path, "w") as f:
            f.write("".join(f'export {i}={env[i]}\n' for i in sorted(env)))

        return job_path

    def _submit_analyses(self, command_tuples):
        """
        Submit pipelines as arrays grouped by the target methods.

        Arguments:
            command_tuples (list): of (analysis, command) tuples.
        """
        raise NotImplementedError

    def _echo_summary(self, invalid_tuples, command_tuples, skipped_tuples):
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

        if self._quiet:
            for i, msg in tuples:
                row = {k: "NA" for k in summary_keys}
                targets = getattr(i, "as_target_objects", i[0])
                references = getattr(i, "as_reference_objects", i[1])
                row["PK"] = getattr(i, "pk", "NA")
                row["MESSAGE"] = self._style_msg(msg)
                row["PROJECTS"] = self._style_projects(targets)
                row["TARGETS"] = self._style_workflows(targets)
                row["REFERENCES"] = self._style_workflows(references)
                summary.append("\t".join(row[k] for k in summary_keys))

        summary = "\n".join(summary).expandtabs(12) + "\n"
        summary += f"RUNNING {len(command_tuples)} | "
        summary += f"SKIPPED {len(skipped_tuples)} | "
        summary += f"INVALID {len(invalid_tuples)}"
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
            for j in i.projects:
                keys.add(f"{j} {i.technique.method}")

        return f"Projects {' '.join(keys)}"

    @staticmethod
    def _style_workflows(workflows):
        """Color message depending on its nature."""
        if len(workflows) > 2 or not workflows:
            message = f"{len(workflows)} workflows."
        else:
            message = " ".join([i.leukid for i in workflows])

        return message
