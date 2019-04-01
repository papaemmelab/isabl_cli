from click.testing import CliRunner

from isabl_cli import cli


def test_main():
    runner = CliRunner()
    result = runner.invoke(cli.main, [])
    assert "process-finished" in result.output
    assert "import-data" in result.output
    assert "import-bedfiles" in result.output
    assert "import-reference-data" in result.output
