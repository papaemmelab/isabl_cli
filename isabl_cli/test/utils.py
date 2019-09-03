import json
import uuid

from isabl_cli import api
from isabl_cli import factories


def create_experiment(
    bam,
    bedfile="/a/fake/bedfile",
    assembly="GRCh37",
    raw_data=None,
    method="TD",
    species="HUMAN",
    category=None,
    technique_name=None,
    sample=None,
):
    """Easily create an experiment for testing purposes."""
    # create technique
    technique = factories.TechniqueFactory(
        method=method, name=technique_name or str(uuid.uuid4())
    )

    # create experiment
    bed_file_dict = dict(url=bedfile, description="")
    technique["reference_data"][f"{assembly}_targets_bedfile"] = bed_file_dict
    technique["reference_data"][f"{assembly}_baits_bedfile"] = bed_file_dict
    experiment = factories.ExperimentFactory(technique=technique)

    # force sample if provided
    if sample:
        experiment["sample"] = sample
    else:
        category = category or experiment["sample"]["category"]
        experiment["sample"]["individual"]["species"] = species
        experiment["sample"]["category"] = category

    experiment["raw_data"] = raw_data or []
    experiment["bam_files"][assembly] = dict(url=bam, analysis=1)
    return api.create_instance("experiments", **experiment)


def create_pair(
    tumor_bam,
    normal_bam,
    bedfile="/a/fake/bedfile",
    assembly="GRCh37",
    method="TD",
    species="HUMAN",
):
    """Create pair."""
    pair = []
    bed_file_dict = dict(url=bedfile, description="")
    technique = factories.TechniqueFactory(method=method)
    technique["reference_data"][f"{assembly}_targets_bedfile"] = bed_file_dict
    technique["reference_data"][f"{assembly}_baits_bedfile"] = bed_file_dict

    for (i, category) in [(tumor_bam, "TUMOR"), (normal_bam, "NORMAL")]:
        experiment = factories.ExperimentFactory(technique=technique)
        experiment["sample"]["individual"]["species"] = species
        experiment["sample"]["category"] = category
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
