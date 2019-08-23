import json
import uuid

from isabl_cli import api
from isabl_cli import factories


def create_experiment(
    bam,
    bedfile="/a/fake/bedfile",
    assembly="GRCh37",
    sequencing_data=None,
    method="TD",
    analyte="DNA",
    species="HUMAN",
    sample_class=None,
    technique_name=None,
    sample=None,
):
    """Easily create a sequencing experiment for testing purposes."""
    # create sequencing technique
    technique = factories.TechniqueFactory(
        method=method, analyte=analyte, name=technique_name or str(uuid.uuid4())
    )

    # create experiment
    technique["bed_files"][assembly] = dict(targets=bedfile, baits=bedfile)
    experiment = factories.ExperimentFactory(technique=technique)

    # force sample if provided
    if sample:
        experiment["sample"] = sample
    else:
        sample_class = sample_class or experiment["sample"]["sample_class"]
        experiment["sample"]["individual"]["species"] = species
        experiment["sample"]["sample_class"] = sample_class

    experiment["sequencing_data"] = sequencing_data or []
    experiment["bam_files"][assembly] = dict(url=bam, analysis=1)
    return api.create_instance("experiments", **experiment)


def create_pair(
    tumor_bam,
    normal_bam,
    bedfile="/a/fake/bedfile",
    assembly="GRCh37",
    method="TD",
    analyte="DNA",
    species="HUMAN",
):
    """Create pair."""
    pair = []
    technique = factories.TechniqueFactory(method=method, analyte=analyte)
    technique["bed_files"][assembly] = dict(targets=bedfile, baits=bedfile)

    for i in tumor_bam, normal_bam:
        experiment = factories.ExperimentFactory(technique=technique)
        experiment["sample"]["individual"]["species"] = species
        experiment["bam_files"][assembly] = dict(url=i, analysis=1)
        pair.append(api.create_instance("experiments", **experiment))

    return pair[0], pair[1]


def create_test_result(
    application=None, results=None, targets=None, references=None, analyses=None
):
    """Return an analysis object."""
    return api.create_instance(
        "analyses",
        **factories.AnalysisFactory(
            targets=targets or [],
            references=references or [],
            analyses=analyses or [],
            application=application or factories.ApplicationFactory(),
            status="SUCCEEDED",
            results=results or {},
        ),
    )


def assert_run(
    application,
    tuples,
    commit,
    results=None,
    project_results=None,
    assert_valid=True,
    assert_skipped=False,
    assert_invalid=False,
):
    """Run application, check results, and return analyses."""
    ret = []
    valid, skipped, invalid = application.run(tuples, commit=commit)

    if assert_valid:
        assert valid, "No valid RAN analyses."

        for i, status in valid:
            ret.append(i)
            results = results or []

            if not commit:
                assert status == application.STAGED_MSG
                continue

            for j in results:
                assert i["results"].get(j) is not None, (
                    f"Result {j} is missing in: "
                    + json.dumps(i["results"], sort_keys=True, indent=4)
                )

    if assert_skipped:
        assert skipped, "No SKIPPED analyses."

    if assert_invalid:
        assert invalid, "No INVALID analyses."

    if project_results:
        project = tuples[0][0][0].projects[0].pk
        analyses = api.get_instances(
            endpoint="analyses", project_level_analysis=project, limit=1
        )

        if analyses:
            assert (
                analyses[0]["status"] == "SUCCEEDED"
            ), f"Project level analysis {analyses[0].pk} has not SUCCEEDED status."

            for j in project_results:
                assert analyses[0]["results"].get(j) is not None, (
                    f"Result {j} is missing in: "
                    + json.dumps(analyses[0]["results"], sort_keys=True, indent=4)
                )

    return ret
