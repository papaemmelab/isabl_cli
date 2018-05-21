"""Engine runner."""

from collections import defaultdict
from datetime import datetime
from itertools import chain
from os.path import abspath
from os.path import dirname
from os.path import join
import getpass
import os
import random
import re
import subprocess
import sys

from click import progressbar
import click

from leuktools import Analysis
from leuktools import constants
from leuktools import exceptions
from leuktools import settings
from leuktools import utils
from leuktools.engine.creator import Creator


HEAD_COMAND_TEMPLATE = r"""#!/bin/bash
{{
    {command}
}} || {{
    exit 1
}}
"""


def print_title(title, color="cyan", blink=False):
    """Print a title."""
    title = "\n" + title.strip().upper() + "\n"
    title += "".join("-" for i in title.strip()) + "\n"
    click.secho(title, fg=color, file=sys.stderr, blink=blink)


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
    _skip_status = {"FAILED", "FINISHED", "STARTED", "SUBMITTED", "SUCCEEDED"}
    _skip_exceptions = (
        exceptions.UsageError,
        exceptions.MissingRequirementError,
        exceptions.ConfigurationError,
        exceptions.MissingOutputError,
        )

    @staticmethod
    def build_command(analysis):
        """Must return string, see pipeline.py module for documentation."""
        raise exceptions.NotImplementedError()

    @staticmethod
    def get_requirements(sequencing_methods):
        """Must return string, See pipeline.py module for documentation."""
        raise exceptions.NotImplementedError()

    def _run_tuples(self, tuples):
        """
        Run a list of tuples.

        Arguments:
            tuples (list): list of (targets, references, analyses) tuples.
                elements can be objects or identifiers.
        """
        tuples = list(tuples)
        print_title(f"Running {len(tuples)} tuples for {self.pipeline}")

        if not tuples:
            return None

        analyses, invalid_tuples = self._get_or_create_analyses(tuples)
        command_tuples, skipped_tuples = [], []

        with progressbar(analyses, label="Building commands...\t\t") as bar:
            for i in bar:
                if i.status in self._skip_status:
                    skipped_tuples.append((i, i.status))
                    continue

                try:
                    if self._force: i.wipe()
                    command = self._build_command(i)
                    command_tuples.append((i, command))
                except self._skip_exceptions as error:
                    skipped_tuples.append((i, error))

        self._echo_summary(invalid_tuples, command_tuples, skipped_tuples)

        if self._commit:
            self._submit_analyses(command_tuples)
        else:
            click.secho("\nAdd --commit to submit.\n", fg="green", blink=True)

    def _build_command(self, analysis):
        """
        Get analysis command as a path to a bash script.

        Arguments:
            analysis (leuktools.Analysis): the analysis object to be run.

        Returns:
            str: analysis command as a path to a bash script.
        """
        job_path = join(analysis.rundir, "head_job.sh")
        env_path = join(analysis.rundir, "head_job.env")
        move = "echo 'NO MOVE COMMAND'"

        # at this point, this variable is required by all pipelines
        # env["LOGS_LSF"] = join(env["RUNDIR"], "logs_lsf")
        # os.makedirs(env["LOGS_LSF"], exist_ok=True)

        # create logs dir for LSF and define an informative jobname
        command, status, env = self.build_command(analysis)  # pylint: disable=E1111
        env.update(dict(RUNDIR=analysis.rundir, JOBNAME=analysis.jobname))

        # final status is FINISHED unless being run by admin user
        if getpass.getuser() != settings.ADMIN_USER and status == "SUCCEEDED":
            status = "FINISHED"

            if analysis.rundir != analysis.outdir:
                move = utils.get_move_command(analysis.rundir, analysis.outdir)

        # add a random sleep so many parallel jobs don't hit api at same time
        command = (
            f"sleep {random.uniform(0, 5):.3} && "
            f"cd {analysis.rundir} && source {env_path} && "
            f"{analysis.get_patch_command('STARTED')} && date && "
            f"{command} && {move} && {analysis.get_patch_command(status)}"
            )

        with open(job_path, "w") as f:
            f.write(HEAD_COMAND_TEMPLATE.format(command=command))

        with open(env_path, "w") as f:
            f.write("".join(f'export {i}={env[i]}\n' for i in sorted(env)))

        return job_path

    def _submit_analyses(self, command_tuples):
        """
        Submit pipelines as arrays grouped by the target methods.

        Arguments:
            command_tuples (list): list of (analysis, command) tuples.
        """
        label = "Grouping analyses...\t\t"
        groups = defaultdict(list)

        # group analyses by the methods of its targets
        with progressbar(command_tuples, file=sys.stderr, label=label) as bar:
            for analysis, command in bar:
                analysis.command = command
                targets = analysis.as_target_objects
                key = tuple(sorted({i.technique.method for i in targets}))
                groups[key].append(analysis)

        # execute analyses on a methods basis
        for methods, analyses in groups.items():
            click.echo(f"Sumbitting {len(analyses)} {methods} jobs...")
            commands, projects = [], {}

            for i in analyses:
                commands.append((i.command, i.get_patch_command("FAILED")))
                projects.update(i.projects)

            jobname = (
                f"Array | pipeline: {self.pipeline} | "
                f"methods: {methods} | projects: {projects}"
                )

            try:
                requirements = self.get_requirements(methods)  # pylint: disable=E1111
                Analysis.update_analyses_status(analyses, "SUBMITTED")
                self._submit_array(commands, requirements, jobname)
            except Exception:  # pylint: disable=broad-except
                Analysis.update_analyses_status(analyses, "STAGED")
                raise Exception("Error during submission...")

    def _submit_array(self, commands, requirements, jobname):
        """
        Submit an array of bash scripts.

        Two other jobs will also be submitted:

            EXIT: run exit command if failure.
            CLEAN: clean temporary files and directories after completion.

        Arguments:
            commands (list): of (path to bash script, on exit command) tuples.
            requirements (str): string of LSF requirements.
            jobname (str): lsf array jobname.

        Returns:
            str: jobid of clean up job.
        """
        root = join(
            settings.PATHS.LEUKDC_DIR, ".runs", getpass.getuser(),
            datetime.now(constants.TIMEZONE).isoformat()
            )

        os.makedirs(root, exist_ok=True)
        jobname += " | rundir: {}".format(root)
        total = len(commands)
        index = 0

        for command, exit_command in commands:
            index += 1
            rundir = abspath(dirname(command))

            with open(join(root, "in.%s" % index), "w") as f:
                f.write(f"bash {command}")

            with open(join(root, "exit_cmd.%s" % index), "w") as f:
                f.write(exit_command)

            for j in "log", "err", "exit":
                src = join(rundir, "head_job.{}".format(j))
                dst = join(root, "{}.{}".format(j, index))
                open(src, "w").close()
                utils.force_symlink(src, dst)

        # submit array of commands
        cmd = (
            f'bsub -J "{jobname}[1-{total}]%50" {requirements} '
            f'-oo "{root}/log.%I" -eo "{root}/err.%I" -i "{root}/in.%I" bash'
            )

        jobid = subprocess.check_output(cmd, shell=True).decode("utf-8")
        jobid = re.findall("<(.*?)>", jobid)[0]

        # submit array of exit commands
        cmd = (
            f'bsub -J "EXIT: {jobname}[1-{total}]" -ti -o "{root}/exit.%I" '
            f'-w "exit({jobid}[*])" -i "{root}/exit_cmd.%I" bash '
            )

        jobid = subprocess.check_output(cmd, shell=True).decode("utf-8")
        jobid = re.findall("<(.*?)>", jobid)[0]

        # clean the execution directory
        cmd = f'bsub -J "CLEAN: {jobname}" -w "ended({jobid})" -ti rm -r {root}'
        jobid = subprocess.check_output(cmd, shell=True).decode("utf-8")
        jobid = re.findall("<(.*?)>", jobid)[0]

        return jobid

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
