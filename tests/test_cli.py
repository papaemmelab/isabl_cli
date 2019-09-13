from click.testing import CliRunner

from isabl_cli import cli
from isabl_cli import settings

from .test_app import MockApplication


def test_main():
    settings._DEFAULTS["INSTALLED_APPLICATIONS"] = [
        "tests.test_app.MockApplication",
        # used to increase coverage
        "isabl_cli.api.IsablDict",
    ]
    cli.add_apps_groups(settings.system_settings.INSTALLED_APPLICATIONS)
    runner = CliRunner()
    result = runner.invoke(cli.main)
    assert "process-finished" in result.output
    assert "import-data" in result.output
    assert "import-bedfiles" in result.output
    assert "import-reference-data" in result.output
    assert f"{MockApplication.ASSEMBLY} applications." in result.output
    settings._DEFAULTS["INSTALLED_APPLICATIONS"] = []
