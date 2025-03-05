from click.testing import CliRunner

from isabl_cli import api
from isabl_cli import commands
from isabl_cli import data
from isabl_cli import factories

from .test_app import MockApplication
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

    # use two experiments to increase coverage with project_results=
    project = api.create_instance("projects", **factories.ProjectFactory())
    experiment = factories.ExperimentFactory(projects=[project])
    experiment["sample"]["individual"]["species"] = "HUMAN"
    experiment_b = factories.ExperimentFactory(projects=[project])
    experiment_b["sample"] = experiment["sample"]
    analysis = utils.assert_run(
        application=MockApplication(),
        tuples=[
            ([api.create_instance("experiments", **experiment)], []),
            ([api.create_instance("experiments", **experiment_b)], []),
        ],
        commit=True,
        project_results=["project_result_file"],
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


def test_get_bed():
    runner = CliRunner()
    technique = api.create_instance("techniques", **factories.TechniqueFactory())
    args = [str(technique.pk)]
    result = runner.invoke(commands.get_bed, args, catch_exceptions=False)
    assert "No BED files" in result.output

    api.patch_instance(
        "techniques",
        technique.pk,
        reference_data={"test_targets_bedfile": {"url": "/hello/world"}},
    )

    result = runner.invoke(commands.get_bed, args, catch_exceptions=False)
    assert "/hello/world" in result.output

    api.patch_instance(
        "techniques",
        technique.pk,
        reference_data={
            "test_targets_bedfile": {"url": "/hello/world"},
            "another_targets_bedfile": {"url": "/hello/world"},
        },
    )

    result = runner.invoke(commands.get_bed, args, catch_exceptions=False)
    assert "Multiple BEDs" in result.output


def test_get_bams():
    runner = CliRunner()
    experiment = api.create_instance("experiments", **factories.ExperimentFactory())
    args = [str(experiment.pk)]
    result = runner.invoke(commands.get_bams, args, catch_exceptions=False)
    assert "No bams for" in result.output

    result = runner.invoke(
        commands.get_bams, args + ["--verbose"], catch_exceptions=False
    )
    assert experiment.system_id in result.output
    assert "None" in result.output

    api.patch_instance(
        "experiments",
        experiment.pk,
        bam_files={"grch": {"url": "/hello/world", "analysis": 1}},
    )

    result = runner.invoke(commands.get_bams, args, catch_exceptions=False)
    assert "/hello/world" in result.output

    api.patch_instance(
        "experiments",
        experiment.pk,
        bam_files={
            "a1": {"url": "/hello/world", "analysis": 1},
            "a2": {"url": "/hello/mars", "analysis": 2},
        },
    )

    result = runner.invoke(commands.get_bams, args, catch_exceptions=False)
    assert "Multiple bams" in result.output

    result = runner.invoke(
        commands.get_bams, args + ["--assembly", "a2"], catch_exceptions=False
    )
    assert "/hello/mars" in result.output


def test_get_data():
    runner = CliRunner()
    experiment = api.create_instance("experiments", **factories.ExperimentFactory())
    experiment = data.update_storage_url("experiments", experiment.pk)
    args = [str(experiment.pk)]
    result = runner.invoke(commands.get_data, args, catch_exceptions=False)
    assert "No data for" in result.output

    result = runner.invoke(
        commands.get_bams, args + ["--verbose"], catch_exceptions=False
    )
    assert experiment.system_id in result.output
    assert "None" in result.output

    api.patch_instance(
        "experiments",
        experiment.pk,
        raw_data=[
            {"file_url": "/hello/world", "file_type": "TXT"},
            {"file_url": "/hello/mars", "file_type": "PNG"},
        ],
    )

    result = runner.invoke(commands.get_data, args, catch_exceptions=False)
    assert "/hello/world" in result.output
    assert "/hello/mars" in result.output

    result = runner.invoke(
        commands.get_data, args + ["--dtypes", "TXT"], catch_exceptions=False
    )
    assert "/hello/mars" not in result.output


def test_login():
    runner = CliRunner()
    result = runner.invoke(commands.login, input="admin\nadmin", catch_exceptions=False)
    assert "Successful authorization! Token stored" in result.output


def test_run_signals(tmpdir):
    signal = "isabl_cli.data.symlink_analysis_to_targets"
    runner = CliRunner()
    experiment_dir = tmpdir.mkdir("experiment")
    analysis = api.create_instance(
        "analyses",
        **factories.AnalysisFactory(
            storage_url=str(tmpdir.mkdir("analysis")),
            targets=[factories.ExperimentFactory(storage_url=str(experiment_dir))],
        ),
    )

    api.patch_instance("analyses", analysis.pk, status="CREATED")
    args = ["analyses", "-fi", "pk", analysis.pk, "-s", signal]
    result = runner.invoke(commands.run_signals, args, catch_exceptions=False)
    assert "analyses" not in str(experiment_dir.listdir())

    api.patch_instance("analyses", analysis.pk, status="SUCCEEDED")
    args = ["analyses", "-fi", "pk", analysis.pk, "-s", signal]
    result = runner.invoke(commands.run_signals, args, catch_exceptions=False)
    assert analysis.application.name.lower() in str(
        experiment_dir.join("analyses").listdir()
    )


def test_run_web_signals():
    application = MockApplication().application
    analysis = api.create_instance(
        "analyses",
        **factories.AnalysisFactory(
            targets=[factories.ExperimentFactory()], application=application
        ),
    )

    api.create_instance(
        "signals",
        import_string="isabl_cli.signals.resume_analysis_signal",
        target_endpoint="analyses",
        target_id=analysis.pk,
    )

    runner = CliRunner()
    args = ["-fi", "target_id", analysis.pk, "-fi", "target_endpoint", "analyses"]
    result = runner.invoke(commands.run_web_signals, args, catch_exceptions=False)
    assert str(analysis.pk) in result.output
    assert "SUCCEEDED" in result.output

    api.create_instance(
        "signals",
        import_string="isabl_cli.signals.force_analysis_signal",
        target_endpoint="analyses",
        target_id=analysis.pk,
    )

    result = runner.invoke(commands.run_web_signals, args, catch_exceptions=False)
    assert str(analysis.pk) in result.output
    assert "SUCCEEDED" in result.output

    # increase coverage on get_result
    assert MockApplication().get_result(
        experiment=api.get_experiments([analysis.targets[0].pk])[0],
        application_key=application.pk,
        result_key="analysis_result_key",
        application_name=str(MockApplication),
    )


def test_process_finished_tags(tmpdir):
    """Check that a tagged analysis does NOT get processed to FINISHED."""
    analysis = api.create_instance(
        "analyses",
        project_level_analysis=factories.ProjectFactory(),
        storage_url=tmpdir.strpath,
        status="FINISHED",
        **factories.AnalysisFactory(ran_by=None, tags=[{"name":"PROCESSING FINISHED"}]),
    )
    runner = CliRunner()
    args = ["-fi", "pk", analysis["pk"]]
    runner.invoke(commands.process_finished, args, catch_exceptions=False)
    analysis = api.get_instance("analyses", analysis["pk"])
    assert analysis["status"] == "FINISHED"

    # strip the tags, check that is DOES get processed
    api.patch_instance("analyses", analysis.pk,
                       **factories.AnalysisFactory(ran_by=None, tags=[]))
    runner.invoke(commands.process_finished, args, catch_exceptions=False)
    analysis = api.get_instance("analyses", analysis["pk"])
    assert analysis["status"] == "SUCCEEDED"


def test_force_process_finished_tags(tmpdir):
    """Check that a tagged analysis does get processed to FINISHED when forced."""
    analysis = api.create_instance(
        "analyses",
        project_level_analysis=factories.ProjectFactory(),
        storage_url=tmpdir.strpath,
        status="FINISHED",
        **factories.AnalysisFactory(ran_by=None, tags=[{"name":"PROCESSING FINISHED"}]),
    )
    runner = CliRunner()
    args = ["-fi", "pk", analysis["pk"], "--force"]
    result = runner.invoke(commands.process_finished, args, catch_exceptions=True)
    print(result.output)
    analysis = api.get_instance("analyses", analysis["pk"])
    assert analysis["status"] == "SUCCEEDED"


def test_rejected_analysis(tmpdir):
    """Check an analysis can be rejected using the command cli."""
    analysis = api.create_instance(
        "analyses",
        status="FINISHED",
        storage_url=tmpdir.strpath,
        **factories.AnalysisFactory(ran_by=None),
    )
    runner = CliRunner()
    rejection_reason = "This analysis is FAKE."
    args = ["--key", analysis["pk"], "--reason", rejection_reason]
    runner.invoke(commands.reject_analysis, args, catch_exceptions=False)
    analysis = api.get_instance("analyses", analysis["pk"])
    assert analysis["status"] == "REJECTED"
    assert analysis["notes"] == rejection_reason
