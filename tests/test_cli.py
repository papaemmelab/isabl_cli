from click.testing import CliRunner

from isabl import cli


def test_main():
    runner = CliRunner()
    result = runner.invoke(cli.main, [])
    assert "patch_status" in result.output
    assert "processed_finished" in result.output
    assert "import_data" in result.output
    assert "import_bedfiles" in result.output
    assert "import_reference_data" in result.output
