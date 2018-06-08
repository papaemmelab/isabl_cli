"""An LSF pipeline engine."""
# pylint: disable=W9008

from collections import defaultdict
from datetime import datetime
from getpass import getuser
from os.path import abspath
from os.path import dirname
from os.path import join
import os
import random
import re
import subprocess

import click

from cli import api
from cli import system_settings
from cli import utils

from .pipeline import AbstractPipeline


class LsfPipeline(AbstractPipeline):  # pragma: no cover

    """Submit jobs as job array throttled by `THROTTLE` argument."""

    THROTTLE = 50

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

    @staticmethod
    def get_job_name(analysis):
        """Get job name given an analysis dict."""
        targets = analysis['targets']
        references = analysis['references']
        methods = {i['technique']['method'] for i in targets}
        projects = {str(j['pk']) for i in targets for j in i['projects']}

        if len(targets) > 2 or not targets:
            targets = f'{len(targets)} samples.'
        else:
            targets = ' '.join([i['system_id'] for i in targets])

        if len(references) > 2 or not references:
            references = f'{len(references)} samples.'
        else:
            references = ' '.join([i['system_id'] for i in references])

        return ' | '.join([
            f'analysis: {analysis["pk"]}',
            f'targets: {targets}',
            f'references: {references} |',
            f'methods: {" ".join(methods)}',
            f'projects: {" ".join(projects)}',
            f'rundir: {analysis["storage_url"]}',
            f'pipeline: {analysis["pipeline"]["pk"]}'])

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
            commands, analyses, projects = [], [], set()

            for i, cmd in command_tuples:
                exit_cmd = f'cli patch_status --key {i["pk"]} --status FAILED'
                analyses.append(i)
                commands.append((cmd, exit_cmd))
                projects.update(
                    j['pk'] for k in i['targets'] for j in k['projects'])

            jobname = (
                f"Array | pipeline: {self.pipeline['pk']} | "
                f"methods: {methods} | projects: {projects}")

            try:
                requirements = self.get_requirements(methods)  # pylint: disable=E1111
                api.patch_analyses_status(analyses, 'SUBMITTED')
                self._submit_array(commands, requirements, jobname)
            except Exception:  # pylint: disable=broad-except
                api.patch_analyses_status(analyses, 'STAGED')
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
        assert system_settings.BASE_STORAGE_DIRECTORY

        root = join(
            system_settings.BASE_STORAGE_DIRECTORY, ".runs", getuser(),
            datetime.now(system_settings.TIME_ZONE).isoformat())

        os.makedirs(root, exist_ok=True)
        jobname += " | rundir: {}".format(root)
        total = len(commands)
        index = 0

        for command, exit_command in commands:
            index += 1
            rundir = abspath(dirname(command))

            with open(join(root, "in.%s" % index), "w") as f:
                # use random sleep to avoid parallel API hits
                f.write(f'sleep {random.uniform(0, 10):.3} && bash {command}')

            with open(join(root, "exit_cmd.%s" % index), "w") as f:
                f.write(exit_command)

            for j in "log", "err", "exit":
                src = join(rundir, "head_job.{}".format(j))
                dst = join(root, "{}.{}".format(j, index))
                open(src, "w").close()
                utils.force_symlink(src, dst)

        # submit array of commands
        cmd = (
            f'bsub -J "{jobname}[1-{total}]%{self.THROTTLE}" {requirements} '
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