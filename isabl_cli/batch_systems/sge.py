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


# catch kill signal using -notify in qsub
COMMAND = """
#!/bin/sh

sigusr2()
{{
  echo kill signal catched
  ({exit_command}) &> {exit_log}
}}

# catch -USR2 signal
trap sigusr2 USR2
trap sigusr2 USR1

# run the command
{command}

# reset to be elegant
trap - USR2
trap - USR1
"""


def submit_sge(app, command_tuples):  # pragma: no cover
    """Submit applications as arrays grouped by the target methods."""
    groups = defaultdict(list)

    # group analyses by the methods of its targets
    for analysis, command in command_tuples:
        targets = analysis["targets"]
        key = tuple(sorted({i["technique"]["method"] for i in targets}))
        groups[key].append((analysis, command))

    # execute analyses on a methods basis
    for methods, cmd_tuples in groups.items():
        submit_configuration = system_settings.SUBMIT_CONFIGURATION
        click.echo(f"Submitting {len(cmd_tuples)} ({', '.join(methods)}) jobs.")
        commands, analyses, projects, extra_args_overrides = [], [], set(), []
        requirements = ""

        # ignore command in tuple and use script instead
        for i, _ in cmd_tuples:
            analyses.append(i)
            exit_cmd = app.get_patch_status_command(i["pk"], "FAILED")
            commands.append((app.get_command_script_path(i), exit_cmd))
            if hasattr(app, "EXTRA_ARGS_OVERRIDE"):
                extra_args_overrides.append(app.EXTRA_ARGS_OVERRIDE)
            else:
                extra_args_overrides.append(submit_configuration.get("extra_args", ""))
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
            for i, j in zip(api.chunks(commands, 10000),  api.chunks(extra_args_overrides, 10000)):
                submit_sge_array(
                    commands=i,
                    requirements=requirements or "",
                    extra_args=j[0],
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


def submit_sge_array(
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
        jobname (str): sge array jobname.
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

    wait = "-sync y" if wait else ""
    os.makedirs(root, exist_ok=True)
    jobname += "-rundir: {}".format(root)
    jobname = slugify(jobname)
    total = len(commands)
    index = 0
    base_args = f"-V -b n -terse -wd {root}"

    with click.progressbar(commands, label="preparing submission...") as bar:
        for command, exit_command in bar:
            index += 1
            rundir = abspath(dirname(command))

            with open(join(root, "in.%s" % index), "w") as f:
                # use random sleep to avoid parallel API hits
                f.write(
                    COMMAND.format(
                        exit_command=exit_command,
                        exit_log=join(rundir, "head_job.exit"),
                        command=f"sleep {random.uniform(0, 10):.3} && bash {command}",
                    )
                )

            for j in "log", "err", "exit":
                src = join(rundir, f"head_job.{j}")
                dst = join(root, f"{j}.{index}")
                open(src, "w").close()
                utils.force_symlink(src, dst)

    with open(join(root, "in.sh"), "w") as f:
        f.write("bash in.${SGE_TASK_ID}")

    with open(join(root, "clean.sh"), "w") as f:
        f.write(f"rm -rf {root}")

    cmd = (
        f"qsub {base_args} {requirements} {extra_args} -t 1-{total} "
        f'-N "ISABL-{jobname}" -notify -tc {throttle_by} '
        f"-o 'log.$TASK_ID' -e 'err.$TASK_ID' {root}/in.sh"
    )

    jobid = subprocess.check_output(cmd, shell=True).decode("utf-8")
    jobid = jobid.strip().split(".")[0]

    cmd = (
        f'qsub {base_args} -N "CLEAN-{jobname}" -hold_jid {jobid} {wait} '
        f"-o /dev/null -e /dev/null {root}/clean.sh"
    )

    return subprocess.check_output(cmd, shell=True).decode("utf-8").strip()
