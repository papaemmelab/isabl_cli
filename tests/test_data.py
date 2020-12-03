import os
import re
import shutil
import uuid

from click.testing import CliRunner
import click
import pytest
import yaml

from isabl_cli import api
from isabl_cli import commands
from isabl_cli import data
from isabl_cli import factories


def test_trash_analysis_storage():
    # the rest is tested in test_engine
    with pytest.raises(click.UsageError) as error:
        data.trash_analysis_storage({"status": "SUCCEEDED"})
    assert "You can't wipe a succeeded analysis" in str(error.value)


def test__make_storage_directory(tmpdir):
    i = data._make_storage_directory(tmpdir.strpath, "test", 12345, use_hash=True)
    j = data._make_storage_directory(tmpdir.strpath, "test", 12345, use_hash=False)
    assert "/test/23/45/12345" in i
    assert "/test/12345" in j


def test_local_reference_genome_importer(tmpdir):
    runner = CliRunner()
    command = data.LocalReferenceGenomeImporter.as_cli_command()
    reference_test = tmpdir.join("test.fasta")
    reference_test.write("foo")
    identifier = str(uuid.uuid4())
    api.create_instance("assemblies", **factories.AssemblyFactory(name=identifier))
    args = [
        "--genome-path",
        reference_test.strpath,
        "--dont-index",
        "--assembly",
        identifier,
    ]
    result = runner.invoke(command, args, catch_exceptions=False)
    print(result.output)


def test_import_reference_data(tmpdir):
    for model, identifier, factory in [
        ("assemblies", str(uuid.uuid4()), factories.AssemblyFactory),
        ("techniques", str(uuid.uuid4()), factories.TechniqueFactory),
    ]:
        reference_test = tmpdir.join("test.fasta")
        reference_test.write("foo")
        api.create_instance(model, **factory(name=identifier))
        instance = data.LocalReferenceDataImporter.import_data(
            identifier=identifier,
            data_src=reference_test.strpath,
            description="test description",
            data_id="reference_link",
            symlink=True,
            model=model,
        )

        assert os.path.islink(instance["reference_data"]["reference_link"]["url"])
        assert (
            instance["reference_data"]["reference_link"]["description"]
            == "test description"
        )

        instance = data.LocalReferenceDataImporter.import_data(
            identifier=identifier,
            data_src=reference_test.strpath,
            description="test description",
            data_id="reference_move",
            model=model,
            symlink=False,
        )

        assert os.path.isfile(instance["reference_data"]["reference_move"]["url"])
        assert not os.path.islink(instance["reference_data"]["reference_move"]["url"])

        reference_test.write("foo")
        command = data.LocalReferenceDataImporter.as_cli_command()
        runner = CliRunner()
        args = [
            "--identifier",
            identifier,
            "--data-src",
            reference_test.strpath,
            "--data-id",
            "reference_move",
            "--description",
            "Test",
            "--model",
            model,
        ]

        result = runner.invoke(command, args, catch_exceptions=False)
        assert "has already reference data registered" in result.output

        args = [str(instance.pk), "--data-id", "reference_link", "--model", model]
        result = runner.invoke(commands.get_reference, args, catch_exceptions=False)
        assert instance["reference_data"]["reference_link"]["url"] in result.output

        args = [str(instance.pk), "--resources", "--model", model]
        result = runner.invoke(commands.get_reference, args, catch_exceptions=False)
        assert "test description" in result.output


@pytest.mark.skipif(
    not shutil.which("bgzip"), reason="bgzip required for this test, install samtools?"
)
def test_import_bedfiles(tmpdir):
    technique = api.create_instance("techniques", **factories.TechniqueFactory())
    targets = tmpdir.join("targets.bed")
    baits = tmpdir.join("baits.bed")
    targets.write("2\t1\t2\n1\t1\t2\n")
    baits.write("2\t1\t2\n1\t1\t2\n")
    species = "HUMAN"
    runner = CliRunner()
    technique = data.LocalBedImporter.import_bedfiles(
        species=species,
        technique=technique["pk"],
        targets_path=targets.strpath,
        baits_path=baits.strpath,
        assembly="AnAssembly",
        description="these are test BED files",
    )

    for i in "targets", "baits":
        bedfile = technique["reference_data"][f"AnAssembly_{i}_bedfile"]["url"]
        assert os.path.isfile(bedfile)
        assert os.path.isfile(bedfile + ".tbi")
        assert os.path.isfile(bedfile.replace(".gz", ""))
        args = [str(technique.pk), "--bed-type", i]
        result = runner.invoke(commands.get_bed, args, catch_exceptions=False)
        assert bedfile in result.output
        with open(bedfile.replace(".gz", ""), "r") as f:  # test bed is sorted
            assert next(f).startswith("1")

    command = data.LocalBedImporter.as_cli_command()
    args = [
        "--targets-path",
        targets.strpath,
        "--baits-path",
        baits.strpath,
        "--technique",
        technique["pk"],
        "--assembly",
        "AnAssembly",
        "--species",
        species,
        "--description",
        "Test",
    ]
    result = runner.invoke(command, args, catch_exceptions=False)
    assert "has registered BED files for" in result.output


def test_local_data_import(tmpdir):
    dirs = [tmpdir.strpath]
    projects = [api.create_instance("projects", **factories.ProjectFactory())]
    experiments = [factories.ExperimentFactory(projects=projects) for i in range(7)]
    experiments = [api.create_instance("experiments", **i) for i in experiments]
    keys = [i["pk"] for i in experiments]

    importer = data.LocalDataImporter()
    # test can't determine type of fastq
    with pytest.raises(click.UsageError) as error:
        path_1 = tmpdir.join(f'{experiments[0]["system_id"]}.fastq')
        path_1.write("foo")
        importer.import_data(directories=dirs, pk__in=keys)

    path_1.remove()
    assert "cant determine fastq type from" in str(error.value)

    # test imports fastq
    path_1 = tmpdir.join(f'{experiments[0]["system_id"]}_R1_foo.fastq')
    path_2 = tmpdir.join(f'{experiments[0]["system_id"]}_R2_foo.fastq')
    path_1.write("foo")
    path_2.write("foo")
    _, summary = importer.import_data(directories=dirs, pk__in=keys, commit=True)
    assert "samples matched: 1" in summary
    assert api.Experiment(experiments[0].pk).get_fastq()

    # test can exclude formats
    path_1 = tmpdir.join(f'{experiments[1]["system_id"]}_1.fastq')
    path_2 = tmpdir.join(f'{experiments[1]["system_id"]}.bam')
    path_1.write("foo")
    path_2.write("foo")
    _, summary = importer.import_data(directories=dirs, pk__in=keys, dtypes=["BAM"])
    assert "FASTQ_R1" not in str(summary)
    assert "BAM" in str(summary)

    # test can import multiple formats
    _, summary = importer.import_data(directories=dirs, pk__in=keys, commit=True)
    assert "FASTQ_R1" in str(summary)
    assert "BAM" in str(summary)

    # test raise error if duplicated ids
    with pytest.raises(click.UsageError) as error:
        api.patch_instance("experiments", experiments[2]["pk"], identifier="dup_id")
        api.patch_instance("experiments", experiments[3]["pk"], identifier="dup_id")
        importer.import_data(
            key=lambda x: x["identifier"], directories=dirs, pk__in=keys
        )

    assert "same identifier for" in str(error.value)

    # test summary
    path_1 = tmpdir.join(f'_{experiments[2]["system_id"]}_cram1_.cram')
    path_2 = tmpdir.join(f'_{experiments[2]["system_id"]}_cram2_.cram')
    path_3 = tmpdir.join(f'_{experiments[3]["system_id"]}_bam1_.bam')
    path_4 = tmpdir.join(f'_{experiments[3]["system_id"]}_bam2_.bam')
    path_1.write("foo")
    path_2.write("foo")
    path_3.write("foo")
    path_4.write("foo")
    imported, summary = importer.import_data(
        directories=dirs, commit=True, symlink=True, pk__in=keys
    )

    project = api.get_instance("projects", projects[0]["pk"])
    assert project["storage_url"]
    assert imported[0]["storage_usage"] > 0
    assert imported[0]["raw_data"]
    assert imported[1]["raw_data"]
    assert "experiments" in imported[1]["storage_url"]
    assert len(os.listdir(os.path.join(imported[1]["storage_url"], "data"))) == 2
    assert "samples matched: 2" in summary
    assert "samples skipped: 2" in summary

    # test import data from command line and files_data functionality
    path_1 = tmpdir.join(f'{experiments[1]["system_id"]}_1.fastq')
    path_2 = tmpdir.join(f'{experiments[1]["system_id"]}_2.fastq')
    path_1.write("foo")
    path_2.write("foo")
    api.patch_instance("experiments", experiments[1]["pk"], raw_data=None)
    file_data = tmpdir.join("file_data.yaml")

    with open(file_data.strpath, "w") as f:
        yaml.dump(
            {
                os.path.basename(path_1.strpath): {"PU": "TEST_PU"},
                os.path.basename(path_2.strpath): {"PU": "TEST_PU"},
            },
            f,
            default_flow_style=False,
        )

    command = data.LocalDataImporter.as_cli_command()
    runner = CliRunner()
    args = [
        "-di",
        tmpdir.strpath,
        "-id",
        "system_id",
        "-fi",
        "pk__in",
        keys,
        "--files-data",
        file_data.strpath,
        "--commit",
    ]

    result = runner.invoke(command, args, catch_exceptions=False)
    assert "samples matched: 1" in result.output
    experiments[1] = api.get_instance("experiments", experiments[1]["pk"])
    assert experiments[1]["raw_data"][0]["file_data"]["PU"] == "TEST_PU"
    assert experiments[1]["raw_data"][1]["file_data"]["PU"] == "TEST_PU"

    # test import using invalid identifier
    args = ["-di", tmpdir.strpath, "-id", "sample", "-fi", "pk__in", keys]
    result = runner.invoke(command, args)
    assert "invalid type for identifier" in result.output

    # test import data using files without read permissions
    path_1 = tmpdir.join(f'{experiments[4]["system_id"]}_1.fastq')
    path_2 = tmpdir.join(f'{experiments[4]["system_id"]}_2.fastq')
    path_1.write("foo")
    path_2.write("foo")
    os.chmod(path_1.strpath, 0o000)
    args = ["-di", tmpdir.strpath, "-id", "system_id", "-fi", "pk__in", keys]
    result = runner.invoke(command, args)
    assert "The following files are not readable by current user:" in result.output
    assert path_1.strpath in result.output
    assert path_2.strpath not in result.output

    # test proper importing from dir with symlinked dir
    fs1 = tmpdir.mkdir("fs1")
    dir_1 = fs1.mkdir("dir_1")
    path_1 = dir_1.join(f'{experiments[5]["system_id"]}_1.fastq')
    path_2 = dir_1.join(f'{experiments[5]["system_id"]}_2.fastq')
    path_1.write("foo")
    path_2.write("foo")

    fs2 = tmpdir.mkdir("fs2")
    dir_2 = fs2.join("dir_2")
    os.symlink(dir_1, dir_2, target_is_directory=True)
    _, summary = importer.import_data(directories=[fs2], pk__in=keys)
    assert "FASTQ_R1" in str(summary)
    assert "FASTQ_R2" in str(summary)

    # test copying data
    path_1 = tmpdir.join(f'{experiments[6]["system_id"]}_R1_foo.fastq')
    path_2 = tmpdir.join(f'{experiments[6]["system_id"]}_R2_foo.fastq')
    path_1.write("foo")
    path_2.write("foo")
    _, summary = importer.import_data(
        directories=dirs, pk=experiments[6]["pk"], commit=True, copy=True
    )
    assert "samples matched: 1" in summary
    assert api.Experiment(experiments[6].pk).get_fastq()
    assert path_1.exists()
    assert path_2.exists()


def test_get_dst():
    importer = data.LocalDataImporter()

    for i, j in [
        # sequencing
        ("sample.bam", "BAM"),
        ("sample.cram", "CRAM"),
        # imaging
        ("sample.png", "PNG"),
        ("sample.jpg", "JPEG"),
        ("sample.jpeg", "JPEG"),
        ("sample.tiff", "TIFF"),
        ("sample.dicom", "DICOM"),
        # text
        ("sample.tsv", "TSV"),
        ("sample.csv", "CSV"),
        ("sample.txt", "TXT"),
        ("sample.tsv.gz", "TSV"),
        ("sample.csv.gz", "CSV"),
        ("sample.txt.gz", "TXT"),
        # other
        ("sample.pdf", "PDF"),
        ("sample.html", "HTML"),
    ]:
        assert data.raw_data_inspector(i) == j
        assert not data.raw_data_inspector(i + "not raw data")

    for test, fq_type in [
        ("sample_R{}_moretext", "R"),
        ("sample_R{}_", "R"),
        ("sample_R{}", "R"),
        ("sample.R{}.more_text", "R"),
        ("sample.R{}.", "R"),
        ("sample.R{}", "R"),
        ("sample_{}", "R"),
        ("sample_S1_L001_I{}_001", "I"),
    ]:
        for index in [1, 2]:
            for fastq in [".fastq", ".fq"]:
                for gzipped in ["", ".gz"]:
                    assert (
                        data.raw_data_inspector(test.format(index) + fastq + gzipped)
                        == f"FASTQ_{fq_type}{index}"
                    )


def test_yaml_data_import(tmpdir):
    # create test data
    projects = [api.create_instance("projects", **factories.ProjectFactory())]
    experiment = factories.ExperimentFactory(projects=projects)
    experiment = api.create_instance("experiments", **experiment)
    data_files_yaml_path = os.path.join(tmpdir, "files_data.yaml")
    file1_path = os.path.join(tmpdir, "test_file_1.fastq.gz")
    file2_path = os.path.join(tmpdir, "test_file_2.fastq.gz")
    data_files_yaml = {
        file1_path: {"key1": "value1", "key2": "value2"},
        file2_path: {"key3": "value3", "key4": "value4"},
    }
    open(file1_path, "w").close()
    open(file2_path, "w").close()
    with open(data_files_yaml_path, "w") as outfile:
        yaml.dump(data_files_yaml, outfile, default_flow_style=False)

    # test import from command line
    command = data.LocalYamlDataImporter.as_cli_command()
    runner = CliRunner()
    args = [
        "-fi",
        "system_id",
        experiment.system_id,
        "--files-data",
        data_files_yaml_path,
        "--commit",
    ]
    runner.invoke(command, args, catch_exceptions=False)

    # reload experiment to get raw_data
    experiment = api.get_instance("experiments", experiment.identifier)

    # verify imported as expected
    assert experiment.raw_data is not None
    assert len(experiment.raw_data) == 2
    assert os.path.basename(file1_path) in experiment.raw_data[0].file_url
    assert experiment.raw_data[0].file_data["key1"] == "value1"
    assert experiment.raw_data[0].file_data["key2"] == "value2"
    assert experiment.raw_data[0].file_type == "FASTQ_R1"
    assert os.path.basename(file2_path) in experiment.raw_data[1].file_url
    assert experiment.raw_data[1].file_data["key3"] == "value3"
    assert experiment.raw_data[1].file_data["key4"] == "value4"
    assert experiment.raw_data[1].file_type == "FASTQ_R2"

    # clean up
    os.remove(data_files_yaml_path)
