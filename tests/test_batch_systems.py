"""Tests for pipelines package utils module."""

from os import path
from os.path import join
import os
import shutil
import pytest

from isabl_cli.batch_systems import lsf
from isabl_cli.batch_systems import sge
from isabl_cli.batch_systems import slurm


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

    # kill SGE with time limit instead of failing directly
    # because SGE's implementation uses SIGUSR2 signal to catch
    # the event when the scheduler kills the job.
    # The other schedulers allow for direct job dependencies.
    extra_args = "" if scheduler != "sge" else "-l h_rt=00:00:10"
    exit_cmd = "exit 1" if scheduler != "sge" else "sleep 20"

    for i in range(total_jobs):
        rundir = path.join(test_dir, str(i))
        os.makedirs(rundir, exist_ok=True)

        # the head job script echoes to stdout and stderr, both will be stored in
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
