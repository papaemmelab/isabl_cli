from os.path import isfile
from os.path import join
import os

from click.testing import CliRunner
import click
import pytest

from isabl_cli import AbstractApplication
from isabl_cli import api
from isabl_cli import commands
from isabl_cli import exceptions
from isabl_cli import factories
from isabl_cli import options
from isabl_cli import utils
from isabl_cli.settings import _DEFAULTS
from isabl_cli.settings import get_application_settings
from isabl_cli.settings import system_settings


def besuhof_signal(instance):
    if "please fail" in instance.notes:
        raise Exception("I was told to fail...")
    elif "fail with different msg" in instance.notes:
        raise Exception("I was told to fail, but with a different msg...")


def test_failed_signal():
    analysis = api.create_instance("analyses", **factories.AnalysisFactory())
    get_kwargs = dict(
        target_endpoint="analyses", endpoint="signals", target_id=analysis.pk
    )

    # check signals work and nothing is created
    api._run_signals("analyses", analysis, [besuhof_signal])
    assert len(api.get_instances(**get_kwargs)) == 0

    # check signals failed
    analysis = api.patch_instance("analyses", analysis.pk, notes="please fail")
    api._run_signals("analyses", analysis, [besuhof_signal])
    instances = api.get_instances(**get_kwargs)
    assert len(instances) == 1
    assert "I was told to fail..." in instances[0].data["failure_traceback"]

    # assert that error traceback is updated
    runner = CliRunner()
    args = f"-fi target_endpoint analyses -fi target_id {analysis.pk}".split()
    api.patch_instance("analyses", analysis.pk, notes="fail with different msg")
    runner.invoke(commands.rerun_signals, args, catch_exceptions=False)
    instances = api.get_instances(**get_kwargs)
    assert len(instances) == 1
    assert "but with a different msg..." in instances[0].data["failure_traceback"]

    # assert that signal is deleted after no failure is detected
    api.patch_instance("analyses", analysis.pk, notes="")
    runner.invoke(commands.rerun_signals, args, catch_exceptions=False)
    assert len(api.get_instances(**get_kwargs)) == 0
