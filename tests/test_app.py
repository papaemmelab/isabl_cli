from os.path import isfile
from os.path import join
import uuid

from cached_property import cached_property
from click.testing import CliRunner
import pytest

from isabl_cli import AbstractApplication
from isabl_cli import api
from isabl_cli import exceptions
from isabl_cli import factories
from isabl_cli import options
from isabl_cli.settings import _DEFAULTS
from isabl_cli.settings import get_application_settings


class NonSequencingApplication(AbstractApplication):

    NAME = str(uuid.uuid4())
    VERSION = "STILL_TESTING"


class ExperimentsFromDefaulCLIApplication(AbstractApplication):

    NAME = str(uuid.uuid4())
    VERSION = "STILL_TESTING"
    ASSEMBLY = "GRCh4000"
    SPECIES = "HUMAN"

    cli_options = [
        options.TARGETS,
        options.ANALYSES,
        options.PAIR,
        options.PAIRS,
        options.PAIRS_FROM_FILE,
        options.NULLABLE_REFERENCES,
    ]

    def validate_experiments(self, targets, references):
        assert not any("raise validation error" in target.notes for target in targets)

    def get_command(*_): # pylint: disable=no-method-argument
        return ""


class MockApplication(AbstractApplication):

    NAME = str(uuid.uuid4())
    VERSION = "STILL_TESTING"
    ASSEMBLY = "GRCh4000"
    SPECIES = "HUMAN"

    cli_help = "This is a test application"
    cli_options = [options.TARGETS]
    application_url = "http://www.fake-test-app.org"
    application_settings = {
        "foo": "bar",
        "test_reference": None,
        "from_system_settings": None,
    }
    application_inputs = {"bar": None}
    application_results = {
        "analysis_result_key": {
            "frontend_type": "number",
            "description": "A random description",
            "verbose_name": "The Test Result",
        }
    }
    application_project_level_results = {
        "project_result_key": {
            "frontend_type": "text-file",
            "description": "A random description",
            "verbose_name": "The Test Result",
        }
    }

    application_individual_level_auto_merge_results = {
        "individual_result_key": {
            "frontend_type": "text-file",
            "description": "A random description",
            "verbose_name": "The Test Result",
        }
    }

    # this is commented out to test get_experiments_from_default_cli_options
    # def get_experiments_from_cli_options(self, targets):
    #     return [([i], []) for i in targets]

    def validate_experiments(self, targets, references):
        self.validate_one_target_no_references(targets, references)

        if targets[0]["identifier"] == "0":
            raise AssertionError("Invalid Center ID")

        return True

    def get_dependencies(self, targets, references, settings):
        return [], {"bar": "foo"}

    def get_command(self, analysis, inputs, settings):
        if settings.restart:
            return "echo successfully restarted"

        if analysis["targets"][0]["identifier"] == "1":
            return "exit 1"

        assert inputs["bar"] == "foo"

        return f"echo {analysis['targets'][0]['system_id']}"

    def merge_project_analyses(self, analysis, analyses):
        assert len(analyses) == 2, f"Expected 2, got: {len(analyses)}"

        with open(join(analysis["storage_url"], "test.merge"), "w") as f:
            f.write(str(len(analyses)))

    def merge_individual_analyses(self, analysis, analyses):
        try:
            with open(join(analysis["storage_url"], "test.merge"), "w") as f:
                f.write(str(len(analyses)))
        except:
            pass

    def get_analysis_results(self, analysis):
        return {"analysis_result_key": 1}

    def get_project_analysis_results(self, analysis):
        # please note that ipdb wont work here as this function will
        # be submitted by a subprocess call
        return {"project_result_key": join(analysis["storage_url"], "test.merge")}

    def get_individual_analysis_results(self, analysis):
        # please note that ipdb wont work here as this function will
        # be submitted by a subprocess call
        return {"individual_result_key": join(analysis["storage_url"], "test.merge")}


class UniquePerIndividualApplication(AbstractApplication):

    NAME = "UNIQUE_PER_INDIVIDUAL"
    VERSION = "STILL_TESTING"
    ASSEMBLY = "GRCh4000"
    SPECIES = "HUMAN"

    cli_options = [options.TARGETS]
    unique_analysis_per_individual = True
    application_results = {
        "analysis_result_key": {
            "frontend_type": "string",
            "description": "A random description",
            "verbose_name": "The Test Result",
        }
    }

    def validate_experiments(self, targets, references):
        return True

    def get_experiments_from_cli_options(self, targets):
        return [([i], []) for i in targets]

    def get_command(self, analysis, inputs, settings):
        return f"echo {analysis.individual_level_analysis.system_id}"

    def get_analysis_results(self, analysis):
        return {"analysis_result_key": analysis.individual_level_analysis.system_id}


class UniquePerIndividualProtectResultsFalse(UniquePerIndividualApplication):

    NAME = "UNIQUE_PER_INDIVIDUAL_NO_PROTECT"
    application_protect_results = False


def submit_merge(instance, application, command):
    assert "merge-individual-analyses" in command
    assert isinstance(application, MockApplication)
    assert instance.pk == 1
    raise ValueError("stop testing here")


def test_custom_submit_merge():
    _DEFAULTS["SUBMIT_MERGE_ANALYSIS"] = "tests.test_app.submit_merge"

    with pytest.raises(ValueError) as error:
        MockApplication().submit_merge_analysis(api.IsablDict(pk=1, species="HUMAN"))

    assert "stop testing here" in str(error)
    _DEFAULTS["SUBMIT_MERGE_ANALYSIS"] = None


def test_assembly_is_not_required():
    app = NonSequencingApplication()
    assert str(app) == f"{app.NAME} {app.VERSION}"
    assert app.primary_key == NonSequencingApplication().primary_key


def test_get_application_settings():
    # test can use reference data id
    assert (
        "genome fasta path"
        == get_application_settings(
            defaults=dict(inner_dict=dict(reference="reference_data_id:genome_fasta")),
            settings=dict(),
            reference_data=dict(genome_fasta=dict(url="genome fasta path")),
            import_strings={},
        ).inner_dict.reference
    )

    # test can use import strings
    assert (
        join
        == get_application_settings(
            defaults=dict(inner_dict=dict(os_methods=[])),
            settings=dict(inner_dict=dict(os_methods=["os.path.join", "os.listdir"])),
            reference_data=dict(),
            import_strings={"os_methods"},
        ).inner_dict.os_methods[0]
    )

    # test inner keys are validated (unexpected settings)
    with pytest.raises(exceptions.ConfigurationError) as error:
        get_application_settings(
            defaults=dict(inner_dict=dict(correct_name=[])),
            settings=dict(inner_dict=dict(wrong_name=[])),
            reference_data=dict(),
            import_strings=set(),
        )

    assert "Got unexpected setting 'wrong_name' for 'inner_dict'" in str(error.value)

    # test can skip keys validation
    assert (
        get_application_settings(
            defaults=dict(inner_dict=dict(correct_name=[], skip_check=True)),
            settings=dict(inner_dict=dict(different_name_but_ok=None)),
            reference_data=dict(),
            import_strings=set(),
        ).inner_dict.different_name_but_ok
        is None
    )

    # test keys are not checked for list of dicts
    assert (
        get_application_settings(
            defaults=dict(dicts_list=[dict(any_setting_name_works=True)]),
            settings=dict(dicts_list=[dict(thats_right=True)]),
            reference_data=dict(),
            import_strings=set(),
        )
        .dicts_list[0]
        .thats_right
        is True
    )

    # test not implemented settings
    with pytest.raises(exceptions.ConfigurationError) as error:
        get_application_settings(
            defaults=dict(
                inner=dict(required_setting=NotImplemented),
                reference="reference_data_id:foo",
            ),
            settings=dict(inner=dict(required=None)),
            reference_data=dict(),
            import_strings=set(),
        )

    assert "Missing required setting: 'required_setting'" in str(error.value)
    assert "Missing required setting: 'reference'" in str(error.value)

    # expected dict, got other
    with pytest.raises(exceptions.ConfigurationError) as error:
        get_application_settings(
            defaults=dict(dict_setting=dict(foo="bar")),
            settings=dict(dict_setting="not a dict"),
            reference_data=dict(),
            import_strings=set(),
        )

    assert "Invalid setting expected dict, got: " in str(error.value)


def test_application_settings(tmpdir):
    application = MockApplication()

    # test __repr__
    assert (
        str(application)
        == f"{application.NAME} {application.VERSION} ({application.ASSEMBLY})"
    )

    # test get settings from reference data id
    application.application_settings["test_reference"] = "reference_data_id:test_id"
    application.assembly["reference_data"]["test_id"] = dict(url="FOO")
    assert application.settings.test_reference == "FOO"
    assert application.settings.from_system_settings is None
    application.application_settings["test_reference"] = None

    # assert can patch databased settings
    application.patch_application_settings(from_system_settings="set")
    assert application.settings.from_system_settings == "set"

    # run it again just to get coverage when evaluating that no API update is needed
    application.patch_application_settings(from_system_settings="set")

    # test bad settings
    application = MockApplication()
    application.application.settings = {}  # avoid default errors
    application.application_settings = {
        "test_reference": "reference_data_id:invalid_key",
        "needs_to_be_implemented": NotImplemented,
    }
    with pytest.raises(exceptions.ConfigurationError) as error:
        application.settings

    assert "Missing required setting: 'test_reference'" in str(error.value)
    assert "Missing required setting: 'needs_to_be_implemented'" in str(error.value)


def test_unique_analysis_per_individual_app(tmpdir):
    individual = factories.IndividualFactory(species="HUMAN")
    sample = factories.SampleFactory(individual=individual)
    project = api.create_instance("projects", **factories.ProjectFactory())
    experiments = [
        factories.ExperimentFactory(
            identifier=str(i), sample=sample, projects=[project]
        )
        for i in range(4)
    ]

    experiments = [api.create_instance("experiments", **i) for i in experiments]
    tuples = [(experiments, [])]
    application = UniquePerIndividualApplication()
    ran_analyses, _, __ = application.run(tuples, commit=True)

    assert len(ran_analyses) == 1
    assert "analysis_result_key" in ran_analyses[0][0]["results"]
    assert len(ran_analyses[0][0].targets) == 4
    assert (
        ran_analyses[0][0]["results"].analysis_result_key
        == experiments[0].sample.individual.system_id
    )

    application = UniquePerIndividualProtectResultsFalse()
    ran_analyses, _, __ = application.run(tuples, commit=True)
    assert len(ran_analyses) == 1
    assert "analysis_result_key" in ran_analyses[0][0]["results"]
    assert len(ran_analyses[0][0].targets) == 4

    # test application_protect_results false
    tuples = [(experiments[:2], [])]
    application = UniquePerIndividualProtectResultsFalse()
    ran_analyses, _, __ = application.run(tuples, commit=True)

    assert len(ran_analyses) == 1
    assert "analysis_result_key" in ran_analyses[0][0]["results"]
    assert len(ran_analyses[0][0].targets) == 2

    # test application_protect_results reduce add more samples - dont remove this test
    tuples = [(experiments, [])]
    application = UniquePerIndividualProtectResultsFalse()
    ran_analyses, _, __ = application.run(tuples, commit=True)

    assert len(ran_analyses) == 1
    assert "analysis_result_key" in ran_analyses[0][0]["results"]
    assert len(ran_analyses[0][0].targets) == 4

    runner = CliRunner()
    command = UniquePerIndividualApplication.as_cli_command()
    result = runner.invoke(
        command, ["-fi", "system_id", experiments[0].system_id], catch_exceptions=False
    )
    assert experiments[0].system_id in result.output


def test_engine(tmpdir):
    _DEFAULTS["DEFAULT_LINUX_GROUP"] = "not_a_group"

    individual = factories.IndividualFactory(species="HUMAN")
    sample = factories.SampleFactory(individual=individual)
    project = api.create_instance("projects", **factories.ProjectFactory())

    experiments = [
        factories.ExperimentFactory(
            identifier=str(i), sample=sample, projects=[project]
        )
        for i in range(4)
    ]

    experiments = [api.create_instance("experiments", **i) for i in experiments]
    tuples = [([i], []) for i in experiments]

    command = MockApplication.as_cli_command()
    application = MockApplication()
    ran_analyses, _, __ = application.run(tuples, commit=True)
    target = api.Experiment(ran_analyses[1][0].targets[0].pk)

    assert "analysis_result_key" in ran_analyses[1][0]["results"].keys()
    assert "analysis_result_key" in ran_analyses[2][0]["results"].keys()
    assert f'analysis: {ran_analyses[1][0]["pk"]}' in application.get_job_name(
        ran_analyses[1][0]
    )

    bam = join(tmpdir, "fake.bam")
    application.update_experiment_bam_file(experiments[0], bam, ran_analyses[0][0].pk)
    assert bam in application.get_bams([experiments[0]])

    # get coverage for when there is no need to update the bam again
    assert application.update_experiment_bam_file(
        experiments[0], bam, ran_analyses[0][0].pk
    )

    with pytest.raises(exceptions.ValidationError) as error:
        application.validate_bams(experiments)

    assert (
        f"{experiments[1].system_id} has no registered bam for "
        f"{application.assembly.name}" in str(error.value)
    )

    with pytest.raises(exceptions.ValidationError) as error:
        application.validate_bedfiles(experiments)

    assert f"{experiments[0].system_id} has no registered bedfile" in str(error.value)

    # test that get results work as expected
    assert (
        application.get_results(
            result_key="analysis_result_key",
            experiment=target,
            application_key=application.primary_key,
        )
        == [(1, ran_analyses[1][0].pk)]
    )

    # check assertion error is raised when an invalid result is searched for
    with pytest.raises(AssertionError) as error:
        application.get_results(
            result_key="invalid_result_key",
            experiment=target,
            application_key=application.primary_key,
        )

    assert "Result 'invalid_result_key' not found for analysis" in str(error.value)

    # test options
    runner = CliRunner()
    result = runner.invoke(command, ["--help"], catch_exceptions=False)
    assert "This is a test application" in result.output
    assert "--commit" in result.output
    assert "--force" in result.output
    assert "--quiet" in result.output
    assert "--restart" in result.output
    assert "--url" in result.output

    runner = CliRunner()
    result = runner.invoke(command, ["--url"], catch_exceptions=False)
    assert "http://www.fake-test-app.org" in result.output

    # test get experiments from default cli options
    runner = CliRunner()
    result = runner.invoke(
        command, ["-fi", "system_id", experiments[0].system_id], catch_exceptions=False
    )
    assert experiments[0].system_id in result.output

    # check project level results
    pks = ",".join(str(i["pk"]) for i in experiments)
    args = ["-fi", "pk__in", pks]
    result = runner.invoke(command, args, catch_exceptions=False)
    analysis = application.get_project_level_auto_merge_analysis(project)
    merged = join(analysis["storage_url"], "test.merge")

    assert analysis["status"] == "SUCCEEDED", f"Project Analysis failed {analysis}"
    assert "FAILED" in result.output
    assert "SUCCEEDED" in result.output
    assert "SKIPPED 3" in result.output
    assert "INVALID 1" in result.output
    assert isfile(merged)
    assert "project_result_key" in analysis["results"]

    with open(merged) as f:
        assert f.read().strip() == "2"

    # check individual level results
    analysis = application.get_individual_level_auto_merge_analysis(individual)
    merged = join(analysis["storage_url"], "test.merge")
    assert analysis["status"] == "SUCCEEDED", f"Individual Analysis failed {analysis}"
    assert isfile(merged)
    assert "individual_result_key" in analysis["results"]

    with open(merged) as f:
        assert f.read().strip() == "2"

    # check passing command line args
    args = ["-fi", "pk__in", pks, "--commit", "--force"]
    result = runner.invoke(command, args)
    assert "--commit not required when using --force" in result.output

    args = ["-fi", "pk__in", pks, "--restart", "--force"]
    result = runner.invoke(command, args)
    assert "cant use --force and --restart together" in result.output

    args = ["-fi", "pk__in", pks, "--force"]
    result = runner.invoke(command, args)
    assert "trashing:" in result.output

    args = ["-fi", "pk__in", pks, "--restart", "--quiet"]
    result = runner.invoke(command, args)
    assert "FAILED" not in result.output

    with open(join(ran_analyses[0][0].storage_url, "head_job.log")) as f:
        assert "successfully restarted" in f.read()

    MockApplication.cli_allow_force = False
    MockApplication.cli_allow_restart = False
    result = runner.invoke(
        MockApplication.as_cli_command(), ["--help"], catch_exceptions=False
    )
    assert "--force" not in result.output
    assert "--restart" not in result.output


def test_validate_is_pair():
    application = AbstractApplication()
    application.validate_is_pair([{"pk": 1}], [{"pk": 2}])

    with pytest.raises(AssertionError) as error:
        application.validate_is_pair([{"pk": 1}, {"pk": 2}], [{"pk": 3}])

    assert "Pairs only." in str(error.value)

    with pytest.raises(AssertionError) as error:
        application.validate_is_pair([{"pk": 1}], [{"pk": 1}])

    assert "Target can't be same as reference." in str(error.value)


def test_validate_reference_genome(tmpdir):
    reference = tmpdir.join("reference.fasta")
    required = ".fai", ".amb", ".ann", ".bwt", ".pac", ".sa"
    application = AbstractApplication()

    with pytest.raises(AssertionError) as error:
        application.validate_reference_genome(reference.strpath)

    assert "Missing indexes please run" in str(error.value)

    for i in required:
        tmp = tmpdir.join("reference.fasta" + i)
        tmp.write("foo")

    with pytest.raises(AssertionError) as error:
        application.validate_reference_genome(reference.strpath)

    assert "samtools dict -a" in str(error.value)


def test_validate_fastq_only():
    application = AbstractApplication()
    targets = [{"raw_data": [], "system_id": "FOO"}]

    with pytest.raises(AssertionError) as error:
        application.validate_has_raw_data(targets)

    assert "FOO" in str(error.value)

    targets = [
        {"raw_data": [{"file_type": "BAM"}], "system_id": "FOO"},
        {"raw_data": [{"file_type": "FASTQ_R1"}], "system_id": "BAR"},
    ]

    with pytest.raises(AssertionError) as error:
        application.validate_single_data_type(targets)

    assert "FOO" in str(error.value)

    targets = [{"raw_data": [{"file_type": "BAM"}], "system_id": "FOO"}]

    with pytest.raises(AssertionError) as error:
        application.validate_fastq_only(targets)

    assert "Only FASTQ supported" in str(error.value)


def test_validate_methods():
    application = AbstractApplication()
    targets = [{"technique": {"method": "FOO"}, "system_id": "FOO BAR"}]

    with pytest.raises(AssertionError) as error:
        application.validate_methods(targets, "BAR")

    assert "Only 'BAR' method(s) allowed" in str(error.value)


def test_validate_pdx_only():
    application = AbstractApplication()
    targets = [{"custom_fields": {"is_pdx": False}, "system_id": "FOO"}]

    with pytest.raises(AssertionError) as error:
        application.validate_pdx_only(targets)

    assert "is not PDX" in str(error.value)


def test_validate_are_normals():
    application = AbstractApplication()
    targets = [
        api.isablfy(
            {"sample": {"category": "TUMOR", "system_id": "FOO"}, "system_id": "FOO"}
        )
    ]

    with pytest.raises(AssertionError) as error:
        application.validate_are_normals(targets)

    assert "is not NORMAL" in str(error.value)


def test_validate_dna_rna_only():
    application = AbstractApplication()
    targets = [{"technique": {"category": "DNA"}, "system_id": "FOO"}]

    with pytest.raises(AssertionError) as error:
        application.validate_rna_only(targets)

    assert "is not RNA" in str(error.value)

    targets = [{"technique": {"category": "RNA"}, "system_id": "FOO"}]

    with pytest.raises(AssertionError) as error:
        application.validate_dna_only(targets)

    assert "is not DNA" in str(error.value)


def test_validate_species():
    application = AbstractApplication()
    targets = [{"sample": {"individual": {"species": "MOUSE"}}, "system_id": "FOO"}]

    with pytest.raises(AssertionError) as error:
        application.validate_species(targets)

    assert "species not supported" in str(error.value)


def test_validate_one_target_no_references():
    application = AbstractApplication()
    targets = [{}]
    references = []
    application.validate_one_target_no_references(targets, references)

    with pytest.raises(AssertionError) as error:
        references.append({})
        application.validate_one_target_no_references(targets, references)

    assert "No reference experiments" in str(error.value)


def test_validate_atleast_onetarget_onereference():
    application = AbstractApplication()
    targets = [{}]
    references = [{}]
    application.validate_at_least_one_target_one_reference(targets, references)

    with pytest.raises(AssertionError) as error:
        targets = []
        application.validate_at_least_one_target_one_reference(targets, references)

    assert "References and targets required" in str(error.value)


def test_validate_targets_not_in_references():
    application = AbstractApplication()
    targets = [{"pk": 1, "system_id": 1}]
    references = [{"pk": 2, "system_id": 2}]
    application.validate_targets_not_in_references(targets, references)

    with pytest.raises(AssertionError) as error:
        references = targets
        application.validate_targets_not_in_references(targets, references)

    assert "1 was also used as reference" in str(error.value)


def test_validate_dna_tuples():
    application = AbstractApplication()
    targets = [{"system_id": 1, "technique": {"category": "DNA"}}]
    references = [{"system_id": 2, "technique": {"category": "DNA"}}]
    application.validate_dna_only(targets + references)

    with pytest.raises(AssertionError) as error:
        targets[0]["technique"]["category"] = "RNA"
        application.validate_dna_only(targets + references)

    assert "category is not DNA" in str(error.value)


def test_validate_dna_pairs():
    application = AbstractApplication()
    targets = [{"pk": 1, "technique": {"category": "DNA"}}]
    references = [{"pk": 2, "technique": {"category": "DNA"}}]
    application.validate_dna_pairs(targets, references)


def test_validate_same_technique():
    application = AbstractApplication()
    targets = [{"system_id": 1, "technique": {"slug": "1"}}]
    references = [{"system_id": 2, "technique": {"slug": "1"}}]
    application.validate_same_technique(targets, references)

    with pytest.raises(AssertionError) as error:
        targets = [{"system_id": 1, "technique": {"slug": "2"}}]
        application.validate_same_technique(targets, references)

    assert "Same techniques required" in str(error.value)

    with pytest.raises(AssertionError) as error:
        references.append({"system_id": 3, "technique": {"slug": "2"}})
        application.validate_same_technique(targets, references)

    assert "Expected one technique, got:" in str(error.value)


def test_validate_same_platform():
    application = AbstractApplication()
    targets = [{"system_id": 1, "platform": {"slug": "1"}}]
    references = [{"system_id": 2, "platform": {"slug": "1"}}]
    application.validate_same_platform(targets, references)

    with pytest.raises(AssertionError) as error:
        targets = [{"system_id": 1, "platform": {"slug": "2"}}]
        application.validate_same_platform(targets, references)

    assert "Same platforms required" in str(error.value)

    with pytest.raises(AssertionError) as error:
        references.append({"system_id": 3, "platform": {"slug": "2"}})
        application.validate_same_platform(targets, references)

    assert "Expected one platform, got:" in str(error.value)


def test_get_experiments_from_default_cli_options(tmpdir):
    app = ExperimentsFromDefaulCLIApplication()

    experiments = []
    for i in range(4):
        experiment_factory = factories.ExperimentFactory()
        experiment_factory["sample"]["individual"]["species"] = "HUMAN"
        experiments.append(api.create_instance("experiments", **experiment_factory))

    analysis = api.create_instance(
        "analyses",
        **{
            **factories.AnalysisFactory(),
            "targets": experiments,
            "references": experiments,
        },
    )

    pairs_file = tmpdir.join("pairs.txt")
    pairs_file.write(experiments[1].system_id + "\t" + experiments[0].system_id + "\n")

    command = ExperimentsFromDefaulCLIApplication.as_cli_command()
    runner = CliRunner()

    args = [
        "--pair",
        experiments[0].system_id,
        experiments[1].system_id,
        "--pairs",
        experiments[2].system_id,
        experiments[3].system_id,
        "--targets-filters",
        "pk",
        experiments[3].pk,
        "--references-filters",
        "pk",
        experiments[2].pk,
        "--analyses-filters",
        "pk",
        analysis.pk,
        "--pairs-from-file",
        str(pairs_file),
    ]

    # get coverage for invalid experiments
    api.patch_instance(
        "experiments", experiments[0].system_id, notes="raise validation error"
    )
    result = runner.invoke(command, args, catch_exceptions=False)
    assert experiments[0].system_id in result.output
    assert "RAN 3 | SKIPPED 0 | INVALID 2" in result.output, result.output

    # revalidate experiment on existing analysis
    api.patch_instance("experiments", experiments[0].system_id, notes="")
    result = runner.invoke(command, args, catch_exceptions=False)
    assert experiments[0].system_id in result.output
    assert "RAN 5 | SKIPPED 0 | INVALID 0" in result.output

    # just get coverage for get_job_name
    assert ExperimentsFromDefaulCLIApplication.get_job_name(analysis)


def test_validate_individuals():
    # Test matched analyis
    matched_application = AbstractApplication()

    targets = [
        {"system_id": 1, "sample": {"individual": {"pk": 1, "system_id": "ind1"}}}
    ]
    references = [
        {"system_id": 2, "sample": {"individual": {"pk": 1, "system_id": "ind1"}}}
    ]
    matched_application.validate_individuals(targets, references)

    with pytest.raises(AssertionError) as error:
        targets = [
            {"system_id": 1, "sample": {"individual": {"pk": 2, "system_id": "ind2"}}}
        ]
        matched_application.validate_individuals(targets, references)
    assert "Same individual required:" in str(error.value)

    # Test unmatched analysis
    unmatched_application = AbstractApplication()
    unmatched_application.IS_UNMATCHED = True

    targets = [
        {"system_id": 1, "sample": {"individual": {"pk": 1, "system_id": "ind1"}}}
    ]
    references = [
        {"system_id": 2, "sample": {"individual": {"pk": 2, "system_id": "ind2"}}}
    ]
    unmatched_application.validate_individuals(targets, references)

    with pytest.raises(AssertionError) as error:
        targets = [
            {"system_id": 1, "sample": {"individual": {"pk": 2, "system_id": "ind2"}}}
        ]
        unmatched_application.validate_individuals(targets, references)
    assert "Different individuals required:" in str(error.value)


def test_notify_project_analyst():
    # Test notification for project analysts
    application = MockApplication()

    project = api.create_instance("projects", **factories.ProjectFactory())
    project["analyst"]
    experiments = [factories.ExperimentFactory(projects=[project]) for i in range(4)]
    experiments = [api.create_instance("experiments", **i) for i in experiments]
    analysis = api.create_instance(
        "analyses",
        **{
            **factories.AnalysisFactory(),
            "targets": experiments,
            "references": experiments,
        },
    )

    assert application.notify_project_analyst(analysis, "test", "test").ok


######## Test App Dependencies ###########


class MockFirstVersion(MockApplication):
    NAME = "PRE-PROCESSING"
    VERSION = "v1"

    def get_analysis_results(self, analysis):
        return {"analysis_result_key": 1}


class MockSecondVersion(MockApplication):
    NAME = "PRE-PROCESSING"
    VERSION = "v2"

    def get_analysis_results(self, analysis):
        return {"analysis_result_key": 2}


class PostMockApp(MockApplication):
    application_inputs = {"dependency_key": NotImplemented}

    def get_command(self, analysis, inputs, settings):
        assert inputs["dependency_key"] == 1
        return f"echo {analysis['targets'][0]['system_id']}"

    def get_analysis_results(self, analysis):
        return {"analysis_result_key": 3}


class PostMockAppWithFirstFixedDependency(PostMockApp):
    NAME = "POST-PROCESSING FIXED"
    VERSION = "DEPENDS ON PRE-PROCESSING v1"

    @cached_property
    def dependencies_results(self):
        return [
            {
                "app": MockFirstVersion(),
                "result": "analysis_result_key",
                "name": "dependency_key",
            }
        ]


class PostMockAppWithSecondFixedDependency(PostMockApp):
    NAME = "POST-PROCESSING FIXED"
    VERSION = "DEPENDS ON PRE-PROCESSING v2"

    @cached_property
    def dependencies_results(self):
        return [
            {
                "app": MockSecondVersion(),
                "result": "analysis_result_key",
                "name": "dependency_key",
            }
        ]


class PostMockAppWithFlexibleDependency(PostMockApp):
    NAME = "POST-PROCESSING FLEXIBLE"
    VERSION = "DEPENDS ON PRE-PROCESSING, ANY VERSION"

    @cached_property
    def dependencies_results(self):
        return [
            {
                "app_name": "PRE-PROCESSING",
                "result": "analysis_result_key",
                "name": "dependency_key",
            }
        ]


class PostMockAppWithFlexibleDependencyNoLinked(PostMockAppWithFlexibleDependency):
    VERSION = "DEPENDS ON PRE-PROCESSING, ANY VERSION. NO LINK"

    @cached_property
    def dependencies_results(self):
        return [
            {
                "app_name": "PRE-PROCESSING",
                "result": "analysis_result_key",
                "name": "dependency_key",
                "linked": False,
            }
        ]


def test_get_dependencies():
    individual = factories.IndividualFactory(species="HUMAN")
    sample = factories.SampleFactory(individual=individual)
    project = api.create_instance("projects", **factories.ProjectFactory())
    experiment = factories.ExperimentFactory(
        identifier="test-experiment", sample=sample, projects=[project]
    )
    experiment = api.create_instance("experiments", **experiment)
    tuples = [([experiment], [])]

    # Run dependency first
    application = MockFirstVersion()
    ran_analyses, _, __ = application.run(tuples, commit=True)
    preprocess_analysis, status = ran_analyses[0]
    assert status == "SUCCEEDED"

    experiment = api.get_instance("experiments", experiment.system_id)
    tuples = [([experiment], [])]

    # A) App with fixed dependency succeeds
    application = PostMockAppWithFirstFixedDependency()
    ran_analyses, _, __ = application.run(tuples, commit=True)
    _, status = ran_analyses[0]
    assert status == "SUCCEEDED"

    # B) App with newer version dependency invalid, because version is strict
    application = PostMockAppWithSecondFixedDependency()
    _, _, invalid_analyses = application.run(tuples, commit=True)
    _, validation_error = invalid_analyses[0]
    assert "No results found for application" in validation_error.message

    # C) App with newer version dependency succeeds, because version is flexible
    application = PostMockAppWithFlexibleDependency()
    ran_analyses, _, __ = application.run(tuples, commit=True)
    analysis, status = ran_analyses[0]
    assert status == "SUCCEEDED"
    assert analysis.analyses[0] == preprocess_analysis.pk

    # New analysis doesn't link the previous analyses to it's run
    application = PostMockAppWithFlexibleDependencyNoLinked()
    ran_analyses, _, __ = application.run(tuples, commit=True)
    analysis, status = ran_analyses[0]
    assert status == "SUCCEEDED"
    assert not analysis.analyses
