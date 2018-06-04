"""An LSF pipeline engine."""

from .pipeline import AbstractPipeline

from collections import defaultdict
from datetime import datetime
from itertools import chain
from os.path import abspath
from os.path import dirname
from os.path import join
from getpass import getuser
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


class LsfPipeline(AbstractPipeline):

    @staticmethod
    def get_requirements(sequencing_methods):
        """
        Get submission requirements given a set of targets' methods.

        Arguments:
            sequencing_methods (set): targets sequencing methods.

        Returns:
            str: lsf requirements (e.g. -q test).
        """
        raise NotImplementedError()

    def submit_analyses(self, command_tuples):
        """
        Submit pipelines as arrays grouped by the target methods.

        Arguments:
            command_tuples (list): list of (analysis, command) tuples.
        """
        groups = defaultdict(list)

        # group analyses by the methods of its targets
        for analysis, command in command_tuples:
            targets = analysis['targets']
            key = tuple(sorted({i['technique']['method'] for i in targets}))
            groups[key].append((analysis, command))

        # execute analyses on a methods basis
        for methods, command_tuples in groups.items():
            click.echo(f"Sumbitting {len(command_tuples)} {methods} jobs...")
            commands, projects = [], set()

            for i, cmd in command_tuples:
                exit_cmd = f'cli patch_status --key {i["pk"]} --status FAILED'
                commands.append((cmd, exit_cmd))
                projects.update(j for k in i['targets'] for j in k['projects'])

            jobname = (
                f"Array | pipeline: {self.pipeline['pk']} | "
                f"methods: {methods} | projects: {projects}")

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
        root = join(settings.PATHS.LEUKDC_DIR, ".runs", getpass.getuser(),
                    datetime.now(constants.TIMEZONE).isoformat())

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
            f'-oo "{root}/log.%I" -eo "{root}/err.%I" -i "{root}/in.%I" bash')

        jobid = subprocess.check_output(cmd, shell=True).decode("utf-8")
        jobid = re.findall("<(.*?)>", jobid)[0]

        # submit array of exit commands
        cmd = (
            f'bsub -J "EXIT: {jobname}[1-{total}]" -ti -o "{root}/exit.%I" '
            f'-w "exit({jobid}[*])" -i "{root}/exit_cmd.%I" bash ')

        jobid = subprocess.check_output(cmd, shell=True).decode("utf-8")
        jobid = re.findall("<(.*?)>", jobid)[0]

        # clean the execution directory
        cmd = f'bsub -J "CLEAN: {jobname}" -w "ended({jobid})" -ti rm -r {root}'
        jobid = subprocess.check_output(cmd, shell=True).decode("utf-8")
        jobid = re.findall("<(.*?)>", jobid)[0]

        return jobid
