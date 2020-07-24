"""Tests for pipelines package utils module."""

from os import path
from os.path import join
import copy
import os
import shutil
import subprocess
import tempfile
import time

import pytest

from isabl_cli.batch_systems import lsf
from isabl_cli.batch_systems import sge
from isabl_cli.batch_systems import slurm
import isabl_cli as ii

from .test_app import MockApplication
from .test_app import get_experiments_for_mock_app


@pytest.mark.skipif(not shutil.which("bsub"), reason="bsub not found")
def test_submit_array_lsf(tmpdir):
    _test_submit_commands(tmpdir, "lsf", lsf.submit_lsf_array)


@pytest.mark.skipif(not shutil.which("qsub"), reason="qsub not found")
def test_submit_array_sge(tmpdir):
    _test_submit_commands(tmpdir, "sge", sge.submit_sge_array)


@pytest.mark.skipif(not shutil.which("sbatch"), reason="sbatch not found")
def test_submit_array_slurm(tmpdir):
    _test_submit_commands(tmpdir, "slurm", slurm.submit_slurm_array)


def _test_submit_commands(tmpdir, scheduler, submit_array):
    test_dir = tmpdir.strpath
    commands = []
    total_jobs = 10
    jobname = "test_execute_headjob"

    # kill SGE with time limit instead of failing directly
    # because SGE's implementation uses SIGUSR2 signal to catch
    # the event when the scheduler kills the job.
    # The other schedulers allow for direct job dependencies.
    extra_args = "" if scheduler != "sge" else "-l h_rt=00:00:10"
    exit_cmd = "exit 1" if scheduler != "sge" else "sleep 20"

    for i in range(total_jobs):
        rundir = path.join(test_dir, str(i))
        os.makedirs(rundir, exist_ok=True)

        # the head job script echoes to stdout and stdout, both will be stored in
        # head_job.log and head_job.err respectively
        cmd_content = "echo 'FERMINA' && (>&2 echo 'VALENTINO') "
        command_path = join(rundir, "head_job.sh")

        if i % 2 == 0:
            cmd_content += f"&& {exit_cmd}"

        with open(command_path, "+w") as f:
            f.write(cmd_content)

        # the exit command echoes both stderr and stdout, however
        # both are meant to be stored in the same file, head_job.exit
        commands.append((command_path, "echo 'LORD' && (>&2 echo 'HENRY')"))

    submit_array(
        commands=commands,
        requirements="",
        extra_args=extra_args,
        throttle_by=total_jobs,
        jobname="Test Job Array",
        wait=True,
    )

    for i in range(total_jobs):
        rundir = path.join(test_dir, str(i))

        with open(join(rundir, "head_job.log")) as f:
            assert "FERMINA" in f.read()

        with open(join(rundir, "head_job.err")) as f:
            assert "VALENTINO" in f.read()

        # chech both stderr and stdout are redirected to head_job.exit
        with open(join(rundir, "head_job.exit")) as f:
            content = f.read()

            if i % 2 == 0:
                assert "LORD" in content
                assert "HENRY" in content

            if i % 2 != 0:
                assert not content


def test_slurm_restart_on_node_fail():
    [experiment], project = get_experiments_for_mock_app(1)
    application = MockApplication()
    application.application_settings[
        "submit_analyses"
    ] = "isabl_cli.batch_systems.submit_slurm"

    # the mock application is design to fail if the target identifier is `1`
    # but work if the restart argument is passed, which works great for this test
    experiment = api.patch_instance("experiments", experiment.pk, identifier="1")
    tuples = [([experiment], [])]

    # temporarily set restart on failed
    os.environ["ISABL_SLURM_RESTART_ON_REGEX_PATTERN"] = "FAILED"
    analyses, _, _ = application.run(tuples=tuples, commit=True)
    analysis = analyses[0][0]
    successful_test = False

    import ipdb

    ipdb.set_trace()

    for i in range(9):
        time.sleep(2 ** i)

        if ii.Analysis(analysis.pk).status == "SUCEEDED":
            # check if the analysis was restarted
            with open(join(ran_analyses[0][0].storage_url, "head_job.log")) as f:
                successful_test = "successfully restarted" in f.read()

    assert successful_test
