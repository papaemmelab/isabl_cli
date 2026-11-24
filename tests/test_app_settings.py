"""Pure tests for the settings manifest loader in `isabl_cli.settings`."""

import importlib.util
import subprocess
import sys

import pytest

from isabl_cli.exceptions import ConfigurationError
from isabl_cli.settings import (
    check_paths,
    load_app_settings,
    setting_paths,
    settings_file,
    settings_record,
)

MSK = "https://isabl.mskcc.org/api/v1/"
SHAH = "https://isabl.shahlab.mskcc.org/api/v1/"

APPS = """
class Tool:
    application_settings = {
        "reference": str,
        "threads": int,
        "bindpaths": (str, type(None)),
        "nextflow_params": {"genes": str, "flag": bool},
    }


class ToolUnmatched(Tool):
    pass
"""

SECTION = """\
    head_job:
      WG: {partition: componc_cpu, cores: 2, memory_gb: 16, minutes: 10080}
    reference: /data1/ref.fasta
    threads: 4
    bindpaths: null
    nextflow_params: {genes: /data1/genes.bed, flag: true}
"""


def manifest(*blocks: str, schema: str = "isabl-app-settings/1") -> str:
    """Return a manifest text with one `Tool` section per instance URL."""
    text = f"schema: {schema}\n"
    for url in blocks:
        text += f"{url}:\n  Tool:\n{SECTION}"
    return text


def make_app(tmp_path, monkeypatch, text, name="Tool"):
    """Write a fixture app module and its manifest, then return the class."""
    (tmp_path / "apps.py").write_text(APPS)
    (tmp_path / "settings.yaml").write_text(text)
    module_name = f"fixture_apps_{tmp_path.name}"
    spec = importlib.util.spec_from_file_location(module_name, tmp_path / "apps.py")
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    monkeypatch.setitem(sys.modules, module_name, module)
    spec.loader.exec_module(module)
    return getattr(module, name)


def faults_of(app_class, api_url=MSK):
    """Return the fault lines of the error that the loader raises."""
    with pytest.raises(ConfigurationError) as info:
        load_app_settings(app_class, api_url)
    header, *lines = str(info.value).splitlines()
    assert header.startswith("Invalid settings manifest ")
    assert header.endswith(f"[{api_url} {app_class.__name__}]:")
    return [i.removeprefix("  - ") for i in lines]


def test_one_block(tmp_path, monkeypatch):
    app = make_app(tmp_path, monkeypatch, manifest(MSK))
    values = load_app_settings(app, MSK)
    assert values.reference == "/data1/ref.fasta"
    assert values.bindpaths is None
    assert values.nextflow_params.flag is True
    assert values.head_job.WG.memory_gb == 16


def test_two_blocks_and_trailing_slash(tmp_path, monkeypatch):
    text = manifest(MSK, SHAH).replace("threads: 4", "threads: 8", 1)
    app = make_app(tmp_path, monkeypatch, text)
    assert load_app_settings(app, MSK.rstrip("/")).threads == 8
    monkeypatch.setenv("ISABL_API_URL", SHAH)
    assert load_app_settings(app).threads == 4


def test_settings_file(tmp_path, monkeypatch):
    app = make_app(tmp_path, monkeypatch, manifest(MSK))
    assert settings_file(app) == str(tmp_path / "settings.yaml")
    (tmp_path / "settings.yaml").unlink()
    assert settings_file(app) is None
    assert "no settings.yaml next to" in faults_of(app)[0]


def test_missing_block(tmp_path, monkeypatch):
    app = make_app(tmp_path, monkeypatch, manifest(MSK, SHAH))
    other = "https://isabl.example.org/api/v1/"
    assert faults_of(app, other) == [f"no block for {other} (blocks: {MSK}, {SHAH})"]


def test_missing_section_and_unknown_class(tmp_path, monkeypatch):
    text = manifest(MSK) + "  Ghost: {}\n"
    app = make_app(tmp_path, monkeypatch, text, name="ToolUnmatched")
    assert faults_of(app) == [
        f"unknown class 'Ghost': {app.__module__} defines no such class",
        "no section for ToolUnmatched",
    ]


def test_subclass_inherits_declarations(tmp_path, monkeypatch):
    text = manifest(MSK) + "  ToolUnmatched:\n" + SECTION
    app = make_app(tmp_path, monkeypatch, text, name="ToolUnmatched")
    assert load_app_settings(app, MSK).threads == 4


def test_missing_unknown_and_wrong_type(tmp_path, monkeypatch):
    text = (
        manifest(MSK)
        .replace("    threads: 4\n", "    thread: 4\n")
        .replace("bindpaths: null", "bindpaths: 3")
        .replace("flag: true", "flag: 'yes', extra: 1")
    )
    app = make_app(tmp_path, monkeypatch, text)
    assert faults_of(app) == [
        "missing setting 'threads'",
        "unknown setting 'thread'",
        "setting 'bindpaths' must be str or NoneType, got int 3",
        "setting 'nextflow_params.flag' must be bool, got str 'yes'",
        "unknown setting 'nextflow_params.extra'",
    ]


def test_optional_value_and_nested_map(tmp_path, monkeypatch):
    text = (
        manifest(MSK)
        .replace("bindpaths: null", "bindpaths: /data1:/data1")
        .replace("{genes: /data1/genes.bed, flag: true}", "/data1/genes.bed")
    )
    app = make_app(tmp_path, monkeypatch, text)
    assert faults_of(app) == [
        "setting 'nextflow_params' must be dict, got str '/data1/genes.bed'"
    ]


def test_bool_is_not_int(tmp_path, monkeypatch):
    app = make_app(
        tmp_path, monkeypatch, manifest(MSK).replace("threads: 4", "threads: true")
    )
    assert faults_of(app) == ["setting 'threads' must be int, got bool True"]


def test_head_job_shape(tmp_path, monkeypatch):
    text = manifest(MSK).replace(
        "cores: 2, memory_gb: 16, minutes: 10080}",
        "cores: '2', memory_gb: 16, queue: x}\n      WE: 3",
    )
    app = make_app(tmp_path, monkeypatch, text)
    assert faults_of(app) == [
        "missing setting 'head_job.WG.minutes'",
        "setting 'head_job.WG.cores' must be int, got str '2'",
        "unknown setting 'head_job.WG.queue'",
        "setting 'head_job.WE' must be dict, got int 3",
    ]


def test_head_job_must_be_a_map(tmp_path, monkeypatch):
    text = manifest(MSK).replace(
        "head_job:\n      WG: {partition: componc_cpu, cores: 2, memory_gb: 16, "
        "minutes: 10080}",
        "head_job: WG",
    )
    app = make_app(tmp_path, monkeypatch, text)
    assert faults_of(app) == ["setting 'head_job' must be dict, got str 'WG'"]


def test_alias_and_anchor(tmp_path, monkeypatch):
    text = manifest(MSK).replace(
        "reference: /data1/ref.fasta", "reference: &ref /data1/ref.fasta"
    )
    text += "  ToolUnmatched:\n" + SECTION.replace(
        "reference: /data1/ref.fasta", "reference: *ref"
    )
    app = make_app(tmp_path, monkeypatch, text)
    assert faults_of(app) == [
        "YAML anchor at line 6: copy the value instead",
        "YAML alias at line 13: copy the value instead",
    ]


def test_duplicate_key(tmp_path, monkeypatch):
    app = make_app(tmp_path, monkeypatch, manifest(MSK) + "    threads: 5\n")
    assert faults_of(app) == ["duplicate key 'threads' at line 10"]


def test_wrong_schema(tmp_path, monkeypatch):
    app = make_app(tmp_path, monkeypatch, manifest(MSK, schema="isabl-app-settings/2"))
    assert faults_of(app) == [
        "schema must be 'isabl-app-settings/1', got 'isabl-app-settings/2'"
    ]


def test_yaml_syntax_error(tmp_path, monkeypatch):
    app = make_app(tmp_path, monkeypatch, manifest(MSK) + "  - [\n")
    assert faults_of(app)[0].startswith("YAML syntax error:")


def test_setting_paths():
    values = {
        "bind": "/data1:/data1,/scratch/x",
        "url": "https://isabl.mskcc.org",
        "threads": 4,
        "tools": ["/usr/bin/samtools", "bwa"],
        "nextflow_config": "process {\n  container = '/sif/a.sif'\n}\n"
        'params.genes="/data1/genes.bed"\n// a comment\n',
        "nested": {"dir": "--outdir=/data1/out"},
    }
    assert setting_paths(values) == [
        ("bind", "/data1"),
        ("bind", "/data1"),
        ("bind", "/scratch/x"),
        ("tools.0", "/usr/bin/samtools"),
        ("nextflow_config", "/sif/a.sif"),
        ("nextflow_config", "/data1/genes.bed"),
        ("nested.dir", "/data1/out"),
    ]


def test_check_paths(tmp_path):
    present = tmp_path / "present"
    present.touch()
    values = {"a": str(present), "b": {"c": f"{tmp_path}/absent"}}
    assert check_paths(values) == [f"path does not exist: b.c = {tmp_path}/absent"]
    assert check_paths(values, ["/data1/papaemme/isabl/"]) == [
        f"path outside production locations: a = {present}",
        f"path does not exist: b.c = {tmp_path}/absent",
        f"path outside production locations: b.c = {tmp_path}/absent",
    ]
    assert check_paths({"a": str(present)}, [f"{tmp_path}/"]) == []
    assert check_paths({"a": str(tmp_path)}, [f"{tmp_path}/"]) == []


def test_settings_record(tmp_path, monkeypatch):
    app = make_app(tmp_path, monkeypatch, manifest(MSK))
    values = {"tool": len, "items": ("/data1/a", 2), "plain": None}
    with pytest.raises(ConfigurationError, match="Run apps from a git checkout"):
        settings_record(app, values, MSK)

    def git(*args):
        subprocess.run(["git", "-C", str(tmp_path), *args], check=True)

    git("init", "-q")
    git("add", ".")
    git("-c", "user.name=t", "-c", "user.email=t@t", "commit", "-qm", "init")
    record = settings_record(app, values, MSK)
    assert len(record.pop("commit")) == 40
    assert record == {
        "file": "settings.yaml",
        "instance": MSK,
        "section": "Tool",
        "schema": "isabl-app-settings/1",
        "values": {"tool": "builtins.len", "items": ["/data1/a", 2], "plain": None},
        "dirty": False,
    }
    (tmp_path / "settings.yaml").write_text("changed")
    assert settings_record(app, values, MSK)["dirty"] is True
    (tmp_path / "settings.yaml").unlink()
    git("rm", "-q", "--cached", "settings.yaml")
    record = settings_record(app, values, MSK)
    assert record["file"] is None and record["schema"] is None
