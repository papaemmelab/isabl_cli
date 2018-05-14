"""cli cli tests."""

from click.testing import CliRunner
import pytest

from cli import cli

def test_main():
    """Sample test for main command."""
    message = "This is a test message for the Universe."
    runner = CliRunner()
    params = ["message", message]
    result = runner.invoke(cli.main, params)
    assert message in result.output

