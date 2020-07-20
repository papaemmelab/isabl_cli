# pylint: disable=W9008,R0201

from collections import defaultdict
from datetime import datetime
from getpass import getuser
from os.path import abspath
from os.path import dirname
from os.path import join
import os
import random
import subprocess

from slugify import slugify
import click

from isabl_cli import api
from isabl_cli import utils
from isabl_cli.settings import system_settings
from isabl_cli.settings import perform_import


def submit_slurm(app, command_tuples):  # pragma: no cover
    """Submit applications as arrays grouped by the target methods."""
    groups = defaultdict(list)

    # group analyses by the methods of its targets
    for analysis, command in command_tuples:
        targets = analysis["targets"]
        key = tuple(sorted({i["technique"]["method"] for i in targets}))
        groups[key].append((analysis, command))

    # execute analyses on a methods basis
    for methods, cmd_tuples in groups.items():
        click.echo(f"Submitting {len(cmd_tuples)} ({', '.join(methods)}) jobs.")
        commands, analyses, projects = [], [], set()
        requirements = ""

        # ignore command in tuple and use script instead
        for i, _ in cmd_tuples:
            analyses.append(i)
            exit_cmd = app.get_patch_status_command(i["pk"], "FAILED")
            commands.append((app.get_command_script_path(i), exit_cmd))
            keys = [j["pk"] for k in i["targets"] for j in k["projects"]]
            projects.update(keys)

        # get requirements as a function of the app and target methods
        get_requirements = system_settings.SUBMIT_CONFIGURATION.get("get_requirements")

        if get_requirements:
            name = "SUBMIT_CONFIGURATION.get_requirements"
            get_requirements = perform_import(val=get_requirements, setting_name=name)
            requirements = get_requirements(app, methods)

        try:
            api.patch_analyses_status(analyses, "SUBMITTED")
            submit_configuration = system_settings.SUBMIT_CONFIGURATION
            for i in api.chunks(commands, 10000):
                submit_slurm_array(
                    commands=i,
                    requirements=requirements or "",
                    extra_args=submit_configuration.get("extra_args", ""),
                    throttle_by=submit_configuration.get("throttle_by", 50),
                    jobname=(
                        f"application: {app} | "
                        f"methods: {', '.join(methods)} | "
                        f"projects: {', '.join(map(str, projects))}"
                    ),
                )

        except Exception:  # pylint: disable=broad-except
            api.patch_analyses_status(analyses, "STAGED")
            raise Exception("Error during submission...")

    return [(i, "SUBMITTED") for i, _ in command_tuples]


def submit_slurm_array(
    commands, requirements, jobname, extra_args=None, throttle_by=50, wait=False
):  # pragma: no cover
    """
    Submit an array of bash scripts.

    Two other jobs will also be submitted:

        EXIT: run exit command if failure.
        CLEAN: clean temporary files and directories after completion.

    Arguments:
        commands (list): of (path to bash script, on exit command) tuples.
        requirements (str): string of LSF requirements.
        jobname (str): slurm array jobname.
        extra_args (str): extra LSF args.
        throttle_by (int): max number of jobs running at same time.
        wait (bool): if true, wait until clean command finishes.

    Returns:
        str: jobid of clean up job.
    """
    extra_args = extra_args or ""
    assert system_settings.BASE_STORAGE_DIRECTORY

    root = join(
        system_settings.BASE_STORAGE_DIRECTORY,
        ".runs",
        getuser(),
        datetime.now(system_settings.TIME_ZONE).isoformat(),
    )

    wait = "-W" if wait else ""
    os.makedirs(root, exist_ok=True)
    jobname += "-rundir: {}".format(root)
    jobname = slugify(jobname)
    total = len(commands)
    index = 0

    with click.progressbar(commands, label="preparing submission...") as bar:
        for command, exit_command in bar:
            index += 1
            rundir = abspath(dirname(command))

            with open(join(root, "in.%s" % index), "w") as f:
                # submit a dependency job on failure
                # important when the scheduler kills the head job
                dependency = "${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
                afternotok = (
                    f"sbatch --dependency=afternotok:{dependency} --kill-on-invalid-dep yes "
                    f'-o {join(rundir, "head_job.exit")} -J "EXIT: {dependency}" '
                    f"<< EOF\n#!/bin/bash\n{exit_command}\nEOF\n"
                )

                # use random sleep to avoid parallel API hits
                f.write(
                    f"#!/bin/sh\nsleep {random.uniform(0, 10):.3} && "
                    f"({afternotok}) && bash {command}"
                )

            for j in "log", "err", "exit":
                src = join(rundir, f"head_job.{j}")
                dst = join(root, f"{j}.{index}")
                open(src, "w").close()
                utils.force_symlink(src, dst)

    with open(join(root, "in.sh"), "w") as f:
        f.write(f"#!/bin/sh\nbash {root}/in.$SLURM_ARRAY_TASK_ID")

    with open(join(root, "clean.sh"), "w") as f:
        f.write(f"#!/bin/sh\nrm -rf {root}")

    cmd = (
        f"sbatch {requirements} {extra_args} --array 1-{total}%{throttle_by} "
        f"-o '{root}/log.%a' -e '{root}/err.%a' "
        f'-J "ISABL: {jobname}" --parsable {root}/in.sh'
    )

    jobid = subprocess.check_output(cmd, shell=True).decode("utf-8").strip()

    cmd = (
        f"sbatch -J 'CLEAN: {jobname}' {wait} --kill-on-invalid-dep yes "
        f"-o /dev/null -e /dev/null --dependency=afterany:{jobid} --parsable {root}/clean.sh"
    )

    return subprocess.check_output(cmd, shell=True).decode("utf-8").strip()
