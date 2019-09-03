from click.testing import CliRunner

from isabl_cli import api
from isabl_cli import commands
from isabl_cli import data
from isabl_cli import factories
from .test_app import TestApplication
from isabl_cli.test import utils


def test_commands(tmpdir):
    analysis = api.create_instance(
        "analyses",
        project_level_analysis=factories.ProjectFactory(),
        storage_url=tmpdir.strpath,
        status="FINISHED",
        **factories.AnalysisFactory(ran_by=None),
    )

    path = tmpdir.join("test.path")
    path.write("not empty")

    runner = CliRunner()
    args = ["-fi", "pk", analysis["pk"]]
    runner.invoke(commands.process_finished, args, catch_exceptions=False)
    analysis = api.get_instance("analyses", analysis["pk"])

    assert analysis["status"] == "SUCCEEDED"
    assert analysis["storage_usage"]

    args = ["--key", analysis["pk"], "--status", "STAGED"]
    runner.invoke(commands.patch_status, args, catch_exceptions=False)
    analysis = api.get_instance("analyses", analysis["pk"])
    assert analysis["status"] == "STAGED"

    args = [
        "analyses",
        "-fi",
        "pk",
        analysis["pk"],
        "-f",
        "pk",
        "-f",
        "application.name",
        "-f",
        "application",
        "-f",
        "carlos",
        "-f",
        "invalid.nested_attr",
    ]

    result = runner.invoke(commands.get_metadata, args, catch_exceptions=False)
    assert analysis["application"]["name"] in result.output
    assert "application.name" in result.output
    assert "INVALID KEY (carlos)" in result.output
    assert "INVALID KEY (nested_attr)" in result.output
    result = runner.invoke(
        commands.get_metadata, args + ["--json"], catch_exceptions=False
    )

    args = ["analyses", "-fi", "pk", analysis["pk"], "--pattern", "*.path"]
    result = runner.invoke(commands.get_paths, args, catch_exceptions=False)
    assert tmpdir.strpath in result.output
    assert "test.path" in result.output

    args = ["analyses", "-fi", "pk", analysis["pk"]]
    result = runner.invoke(commands.get_paths, args, catch_exceptions=False)
    assert tmpdir.strpath in result.output

    args = ["analyses", "-fi", "pk", analysis["pk"]]
    result = runner.invoke(commands.get_count, args, catch_exceptions=False)
    assert "1" in result.output

    args = ["-fi", "pk", analysis["pk"]]
    result = runner.invoke(commands.get_outdirs, args, catch_exceptions=False)
    assert tmpdir.strpath in result.output
    result = runner.invoke(
        commands.get_outdirs, args + ["--pattern", "*.path"], catch_exceptions=False
    )
    assert "test.path" in result.output

    project = api.create_instance("projects", **factories.ProjectFactory())
    experiment = factories.ExperimentFactory(projects=[project])
    experiment["sample"]["individual"]["species"] = "HUMAN"
    analysis = utils.assert_run(
        application=TestApplication(),
        tuples=[([api.create_instance("experiments", **experiment)], [])],
        commit=True,
    )[0]

    args = ["--app-results", analysis.application.pk]
    result = runner.invoke(commands.get_results, args, catch_exceptions=False)
    assert "command_script" in result.output
    args = ["-fi", "pk", analysis.pk, "-r", "command_script"]
    result = runner.invoke(commands.get_results, args, catch_exceptions=False)
    assert "head_job.sh" in result.output

    args = ["-fi", "pk", analysis.pk, "--force"]
    result = runner.invoke(commands.patch_results, args, catch_exceptions=False)
    assert "Retrieving 1 from analyses API endpoint" in result.output
