# pylint: disable=W9008,R0201

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

from isabl_cli import api
from isabl_cli import utils
from isabl_cli.settings import system_settings
from isabl_cli.settings import perform_import

# settings from environ variables
LSF_GET_REQUIREMENTS = "ISABL_LSF_GET_REQUIREMENTS"
LSF_EXTRA_ARGS = "ISABL_LSF_EXTRA_ARGS"
LSF_THROTTLE_BY = "ISABL_LSF_THROTTLE_BY"


def submit_lsf(app, command_tuples):  # pragma: no cover
    """Submit applications as arrays grouped by the target methods."""
    groups = defaultdict(list)

    # group analyses by the methods of its targets
    for analysis, command in command_tuples:
        targets = analysis["targets"]
        key = tuple(sorted({i["technique"]["method"] for i in targets}))
        groups[key].append((analysis, command))

    # execute analyses on a methods basis
    for methods, cmd_tuples in groups.items():
        click.echo(f"Submitting {len(cmd_tuples)} {methods} jobs.")
        commands, analyses, projects = [], [], set()
        requirements = ""

        # ignore command in tuple and use script instead
        for i, _ in cmd_tuples:
            analyses.append(i)
            exit_cmd = app.get_patch_status_command(i["pk"], "FAILED")
            commands.append((app.get_command_script_path(i), exit_cmd))
            keys = [j["pk"] for k in i["targets"] for j in k["projects"]]
            projects.update(keys)

        if os.environ.get(LSF_GET_REQUIREMENTS):
            requirements = perform_import(
                val=os.environ.get(LSF_GET_REQUIREMENTS),
                setting_name=LSF_GET_REQUIREMENTS,
            )(app, methods)

        try:
            api.patch_analyses_status(analyses, "SUBMITTED")
            submit_lsf_array(
                commands=commands,
                requirements=requirements,
                extra_args=os.environ.get(LSF_EXTRA_ARGS, ""),
                throttle_by=os.environ.get(LSF_THROTTLE_BY, 50),
                jobname=f"Array | application: {app.primary_key} | "
                f"methods: {methods} | projects: {projects}",
            )

        except Exception:  # pylint: disable=broad-except
            api.patch_analyses_status(analyses, "STAGED")
            raise Exception("Error during submission...")

    return [(i, "SUBMITTED") for i, _ in command_tuples]


def submit_lsf_array(
    commands, requirements, jobname, extra_args="", throttle_by=50
):  # pragma: no cover
    """
    Submit an array of bash scripts.

    Two other jobs will also be submitted:

        EXIT: run exit command if failure.
        CLEAN: clean temporary files and directories after completion.

    Arguments:
        commands (list): of (path to bash script, on exit command) tuples.
        requirements (str): string of LSF requirements.
        jobname (str): lsf array jobname.
        extra_args (str): extra LSF args.
        throttle_by (int): max number of jobs running at same time.

    Returns:
        str: jobid of clean up job.
    """
    assert system_settings.BASE_STORAGE_DIRECTORY

    root = join(
        system_settings.BASE_STORAGE_DIRECTORY,
        ".runs",
        getuser(),
        datetime.now(system_settings.TIME_ZONE).isoformat(),
    )

    os.makedirs(root, exist_ok=True)
    jobname += " | rundir: {}".format(root)
    total = len(commands)
    index = 0

    for command, exit_command in commands:
        index += 1
        rundir = abspath(dirname(command))

        with open(join(root, "in.%s" % index), "w") as f:
            # use random sleep to avoid parallel API hits
            f.write(f"sleep {random.uniform(0, 10):.3} && bash {command}")

        with open(join(root, "exit_cmd.%s" % index), "w") as f:
            f.write(exit_command)

        for j in "log", "err", "exit":
            src = join(rundir, "head_job.{}".format(j))
            dst = join(root, "{}.{}".format(j, index))
            open(src, "w").close()
            utils.force_symlink(src, dst)

    # submit array of commands
    cmd = (
        f"bsub {requirements} {extra_args} "
        f'-J "{jobname}[1-{total}]%{throttle_by}" '
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
