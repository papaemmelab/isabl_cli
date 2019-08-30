from click.testing import CliRunner

from isabl_cli import cli
from .test_app import TestApplication


def test_main():
    application = TestApplication()
    cli.add_apps_groups([application])
    runner = CliRunner()
    result = runner.invoke(cli.main, [])
    assert "process-finished" in result.output
    assert "import-data" in result.output
    assert "import-bedfiles" in result.output
    assert "import-reference-data" in result.output
    assert f"{application.ASSEMBLY} applications." in result.output
