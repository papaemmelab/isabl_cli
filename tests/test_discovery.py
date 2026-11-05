"""Tests for finding applications in installed packages."""

from textwrap import dedent
import importlib
import sys

import pytest

from isabl_cli import discovery

FAKE_PACKAGES = {
    "fakeapps/__init__.py": "",
    "fakeapps/apps/__init__.py": "",
    "fakeapps/apps/align.py": """
        from isabl_cli import AbstractApplication

        class Align(AbstractApplication):
            NAME = "ALIGN"
            VERSION = "1.0.0"
            cli_help = "Align reads."

        class AlignGRCh37(Align):
            ASSEMBLY = "GRCh37"
            SPECIES = "HUMAN"

        class AlignGRCh38(Align):
            ASSEMBLY = "GRCh38"
            SPECIES = "HUMAN"
            _PRIVATE = True

        class NeedsEnvGRCh37(AlignGRCh37):
            def __init__(self):
                raise AttributeError("ISABL_API_URL is not set")

        class FinaliseAlignGRCh37(AlignGRCh37):
            cli_help = "Finalise an alignment."

            def get_cli_command_name(self):
                return super().get_cli_command_name() + "_finalise"

        class NotAnApp:
            NAME = "NOT AN APP"
            VERSION = "1.0.0"
    """,
    "fakeapps/apps/call.py": """
        from isabl_cli import AbstractApplication
        from fakeapps.apps.align import AlignGRCh37
        from otherapps import Foreign

        class Call(AbstractApplication):
            NAME = "CALL"
            VERSION = "0.1"
            ASSEMBLY = "GRCh37"
            SPECIES = "HUMAN"
    """,
    "fakeapps/apps/broken.py": "raise RuntimeError('boom')",
    "fakeapps/apps/script.py": "import sys; sys.exit('usage: script.py <input>')",
    # imports fine only after fakeapps.apps.call has been imported
    "fakeapps/apps/a_needs_call.py": """
        import sys

        if "fakeapps.apps.call" not in sys.modules:
            raise ImportError("circular import")

        from isabl_cli import AbstractApplication

        class Retried(AbstractApplication):
            NAME = "RETRIED"
            VERSION = "1"
    """,
    # a folder without __init__.py can still be imported
    "fakeapps/apps/nested/apps.py": """
        from isabl_cli import AbstractApplication

        class Nested(AbstractApplication):
            NAME = "NESTED"
            VERSION = "1"
            ASSEMBLY = "GRCh37"
    """,
    "otherapps/__init__.py": """
        from isabl_cli import AbstractApplication

        class Foreign(AbstractApplication):
            NAME = "FOREIGN"
            VERSION = "1"
            ASSEMBLY = "GRCh37"
    """,
}


@pytest.fixture
def install_packages(tmp_path, monkeypatch):
    """Return a function that writes packages (and dist-info) into sys.path."""
    added = set()

    def install(files, entry_points=None, dist="fakeapps"):
        files = dict(files)

        if entry_points:  # make importlib.metadata find the entry points
            dist_info = f"{dist}-0.1.dist-info"
            files[f"{dist_info}/METADATA"] = (
                f"Metadata-Version: 2.1\nName: {dist}\nVersion: 0.1\n"
            )
            files[f"{dist_info}/entry_points.txt"] = "".join(
                f"[{group}]\n" + "".join(f"{k} = {v}\n" for k, v in items.items())
                for group, items in entry_points.items()
            )

        for path, content in files.items():
            path = tmp_path / path
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(dedent(content))
            added.add(path.relative_to(tmp_path).parts[0].split("-")[0])

        monkeypatch.syspath_prepend(str(tmp_path))
        importlib.invalidate_caches()
        return tmp_path

    yield install

    for name in list(sys.modules):
        if name.split(".")[0] in added:
            del sys.modules[name]


def test_find_applications(install_packages, capsys):
    install_packages(FAKE_PACKAGES)
    apps = discovery.find_applications("fakeapps.apps")

    # base classes are found unless more attributes are required, private and
    # foreign apps are skipped, and re-exported apps are listed once
    assert [i.__name__ for i in apps] == [
        "Retried",
        "Align",
        "AlignGRCh37",
        "FinaliseAlignGRCh37",
        "NeedsEnvGRCh37",
        "Call",
        "Nested",
    ]

    stderr = capsys.readouterr().err
    assert "Skipping fakeapps.apps.broken" in stderr
    assert "Skipping fakeapps.apps.script" in stderr
    assert "a_needs_call" not in stderr

    apps = discovery.find_applications(
        "fakeapps.apps", required=("NAME", "VERSION", "ASSEMBLY")
    )
    assert [i.__name__ for i in apps] == [
        "AlignGRCh37",
        "FinaliseAlignGRCh37",
        "NeedsEnvGRCh37",
        "Call",
        "Nested",
    ]

    apps = discovery.find_applications("fakeapps.apps", modules="fakeapps.apps.*.apps")
    assert [i.__name__ for i in apps] == ["Nested"]

    # only import some modules, re-exported apps are found where they are imported
    apps = discovery.find_applications("fakeapps.apps", modules="fakeapps.apps.call")
    assert [i.__name__ for i in apps] == ["AlignGRCh37", "Call"]
    assert apps[0].__module__ == "fakeapps.apps.align"

    # metadata points to where apps were found, the module meant to be imported
    specs = discovery.find_application_specs(
        "fakeapps.apps", modules="fakeapps.apps.call"
    )
    assert [i["import_string"] for i in specs] == [
        "fakeapps.apps.call.AlignGRCh37",
        "fakeapps.apps.call.Call",
    ]


def test_application_spec(install_packages):
    install_packages(FAKE_PACKAGES)
    from fakeapps.apps.align import AlignGRCh37  # pylint: disable=C0415
    from fakeapps.apps.align import FinaliseAlignGRCh37  # pylint: disable=C0415

    assert discovery.application_spec(AlignGRCh37) == {
        "import_string": "fakeapps.apps.align.AlignGRCh37",
        "name": "ALIGN",
        "version": "1.0.0",
        "assembly": "GRCh37",
        "command_name": "align-1.0.0",
        "help": "Align reads.",
    }

    # command names come from get_cli_command_name, which apps can override
    spec = discovery.application_spec(FinaliseAlignGRCh37)
    assert spec["command_name"] == "align-1.0.0_finalise"

    # apps that can't be instantiated get the default name
    from fakeapps.apps.align import NeedsEnvGRCh37  # pylint: disable=C0415

    assert discovery.get_command_name(NeedsEnvGRCh37) == "align-1.0.0"

    command_name = AlignGRCh37().get_cli_command_name()
    assert command_name == discovery.cli_command_name("ALIGN", "1.0.0") == "align-1.0.0"
    assert discovery.get_command_name(AlignGRCh37) == command_name
    assert discovery.get_command_name(dict) is None


def test_is_excluded():
    patterns = ["isabl_apps.apps.legacy.*", "bwa-mem-0.0.*"]
    assert discovery.is_excluded("isabl_apps.apps.legacy.apps.Old", "old-1", patterns)
    assert discovery.is_excluded(
        "isabl_apps.apps.bwa.apps.Bwa", "bwa-mem-0.0.1", patterns
    )
    assert not discovery.is_excluded(
        "isabl_apps.apps.bwa.apps.Bwa", "bwa-mem-1.0.0", patterns
    )
    assert not discovery.is_excluded("isabl_apps.apps.bwa.apps.Bwa", None, [])


def test_get_entry_points(install_packages):
    install_packages(
        {"fakeapps/__init__.py": ""},
        entry_points={"isabl_cli.test": {"fakeapps": "fakeapps:info"}},
    )

    entry_points = discovery.get_entry_points("isabl_cli.test")
    assert [(i.name, i.value) for i in entry_points] == [("fakeapps", "fakeapps:info")]
    assert discovery.get_entry_points("isabl_cli.missing") == []
