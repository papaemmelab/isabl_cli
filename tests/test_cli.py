"""cli cli tests."""

from click.testing import CliRunner

from cli import cli


def test_main():
    runner = CliRunner()
    result = runner.invoke(cli.main, [])
    assert 'patch_status' in result.output
    assert 'processed_finished' in result.output
