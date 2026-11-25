from click.testing import CliRunner
import pytest

from isabl_cli import cli
from isabl_cli import settings

from .test_app import MockApplication


def test_main(capsys):
    settings._DEFAULTS["INSTALLED_APPLICATIONS"] = [
        # valid app
        "tests.test_app.MockApplication",
        # invalid app
        "tests.test_app.ImaginaryApplication",
        # used to increase coverage
        "isabl_cli.api.IsablDict",
    ]
    # Prints error trying to import invalid app
    cli.add_apps_groups(settings.system_settings.INSTALLED_APPLICATIONS)
    captured = capsys.readouterr()
    assert "Could not import 'tests.test_app.ImaginaryApplication'" in captured.err

    runner = CliRunner()
    result = runner.invoke(cli.main)
    assert "process-finished" in result.output
    assert "import-data" in result.output
    assert "import-bedfiles" in result.output
    assert "import-reference-data" in result.output

    # Still imports valid app
    assert f"{MockApplication.ASSEMBLY} applications." in result.output

    settings._DEFAULTS["INSTALLED_APPLICATIONS"] = []
