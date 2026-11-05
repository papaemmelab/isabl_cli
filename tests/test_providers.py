"""Tests for applications provided through entry points."""

import sys

from click.testing import CliRunner
import click
import pytest

from isabl_cli import cli
from isabl_cli import providers

from .test_discovery import FAKE_PACKAGES
from .test_discovery import install_packages  # pylint: disable=W0611

PROVIDER_PACKAGES = dict(
    FAKE_PACKAGES,
    **{
        "fakeapps/isabl.py": """
            from isabl_cli.discovery import find_applications

            def eager():
                return {
                    "applications": find_applications(
                        "fakeapps.apps", required=("NAME", "VERSION", "ASSEMBLY")
                    )
                }

            def lazy():
                return {
                    "applications": [
                        {
                            "import_string": "fakeapps.lazyapps.Slow",
                            "name": "SLOW",
                            "version": "2.0",
                            "assembly": "GRCh38",
                            "help": "A slow app.",
                            "added_in_a_future_version": True,
                        },
                        {
                            "import_string": "fakeapps.lazyapps.SlowFinalise",
                            "name": "SLOW",
                            "version": "2.0",
                            "assembly": "GRCh38",
                            "command_name": "slow-2.0_finalise",
                            "help": "Finalise a slow app.",
                        },
                    ]
                }

            def broken():
                raise RuntimeError("provider boom")
        """,
        "fakeapps/lazyapps.py": """
            from isabl_cli import AbstractApplication

            class Slow(AbstractApplication):
                NAME = "SLOW"
                VERSION = "2.0"
                ASSEMBLY = "GRCh38"
                cli_help = "A slow app."

            class SlowFinalise(Slow):
                cli_help = "Finalise a slow app."

                def get_cli_command_name(self):
                    return super().get_cli_command_name() + "_finalise"
        """,
    },
)

ENTRY_POINTS = {
    providers.ENTRY_POINT_GROUP: {
        "eager": "fakeapps.isabl:eager",
        "lazy": "fakeapps.isabl:lazy",
        "broken": "fakeapps.isabl:broken",
    }
}


@pytest.fixture
def fake_providers(install_packages):  # pylint: disable=W0621
    """Install a package with three providers."""
    return install_packages(PROVIDER_PACKAGES, entry_points=ENTRY_POINTS)


def names(apps):
    """Return class or lazy application names."""
    return [getattr(i, "__name__", None) or i.import_string for i in apps]


def test_collect_applications(fake_providers, capsys):  # pylint: disable=W0621
    from fakeapps.apps.align import AlignGRCh37  # pylint: disable=C0415

    # providers are ignored unless discovery is enabled
    eager, lazy = providers.collect_applications([AlignGRCh37], discover=False)
    assert names(eager) == ["AlignGRCh37"] and lazy == []

    # installed apps come first and are not duplicated
    eager, lazy = providers.collect_applications([AlignGRCh37])
    assert names(eager) == [
        "AlignGRCh37",
        "FinaliseAlignGRCh37",
        "NeedsEnvGRCh37",
        "Call",
        "Nested",
    ]
    assert names(lazy) == ["fakeapps.lazyapps.Slow", "fakeapps.lazyapps.SlowFinalise"]
    assert lazy[0].provider == "lazy" and lazy[0].command_name == "slow-2.0"
    assert "Failed to load applications from 'broken'" in capsys.readouterr().err

    # lazy apps are not imported while collecting
    assert "fakeapps.lazyapps" not in sys.modules

    # exclude by command name or import string
    eager, lazy = providers.collect_applications(
        [AlignGRCh37], excluded=["call-*", "fakeapps.lazyapps.*"]
    )
    assert names(eager) == [
        "AlignGRCh37",
        "FinaliseAlignGRCh37",
        "NeedsEnvGRCh37",
        "Nested",
    ]
    assert lazy == []

    # invalid items are reported and skipped
    eager, _ = providers.collect_applications(
        ["fakeapps.apps.align.Missing"], discover=False
    )
    assert eager == []
    assert "Invalid application" in capsys.readouterr().err


def test_lazy_application_group(fake_providers):  # pylint: disable=W0621
    _, lazy = providers.collect_applications()
    group = providers.ApplicationGroup(
        name="apps-grch38", lazy_applications={i.command_name: i for i in lazy}
    )

    runner = CliRunner()
    result = runner.invoke(group, ["--help"])
    assert result.exit_code == 0
    assert "slow-2.0" in result.output and "A slow app." in result.output
    assert "fakeapps.lazyapps" not in sys.modules

    result = runner.invoke(group, ["slow-2.0", "--help"])
    assert result.exit_code == 0, result.output
    assert "A slow app." in result.output
    assert "fakeapps.lazyapps" in sys.modules
    assert "slow-2.0" in group.commands

    # command names can come from metadata, matching get_cli_command_name
    result = runner.invoke(group, ["slow-2.0_finalise", "--help"])
    assert result.exit_code == 0, result.output
    assert "Finalise a slow app." in result.output
    assert "out of date" not in result.output


def test_lazy_application_errors(fake_providers, capsys):  # pylint: disable=W0621
    missing = providers.LazyApplication("fakeapps.lazyapps.Missing", "MISSING", "1")

    with pytest.raises(click.ClickException) as error:
        missing.load_command()

    assert "Failed to load fakeapps.lazyapps.Missing" in str(error.value.message)

    # metadata that doesn't match the class is reported
    stale = providers.LazyApplication(
        "fakeapps.lazyapps.Slow", "SLOW", "2.0", assembly="GRCh37"
    )
    assert stale.load_command().name == "slow-2.0"
    assert "is out of date" in capsys.readouterr().err


def test_add_apps_groups_with_lazy_applications(
    fake_providers,
):  # pylint: disable=W0621
    eager, lazy = providers.collect_applications()
    cli.add_apps_groups(eager, lazy)

    result = CliRunner().invoke(cli.main, ["--help"])
    assert "apps-grch37" in result.output and "apps-grch38" in result.output

    result = CliRunner().invoke(cli.main, ["apps-grch37", "--help"])
    assert "align-1.0.0" in result.output and "call-0.1" in result.output

    # a lazy app is skipped if the same command was already built
    cli.add_apps_groups(
        eager,
        lazy + [providers.LazyApplication("x.Y", "CALL", "0.1", assembly="GRCh37")],
    )
    assert cli.main.commands["apps-grch37"].lazy_applications == {}

    for i in ("apps-grch37", "apps-grch38"):
        cli.main.commands.pop(i)
