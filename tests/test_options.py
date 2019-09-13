from click.testing import CliRunner
import click

from isabl_cli import api
from isabl_cli import factories
from isabl_cli import options

from .test_app import MockApplication


def test_get_analyses_filters_option():
    runner = CliRunner()
    application = MockApplication().application
    analysis = api.create_instance(
        "analyses", **factories.AnalysisFactory(application=application)
    )

    @click.command()
    @options.get_analyses_filters_option([MockApplication])
    def command(analyses_filters):
        for i in analyses_filters:
            print(i.pk)

    result = runner.invoke(command, ["-fi", "pk", analysis.pk])
    assert str(analysis.pk) in result.output

    result = runner.invoke(command, ["--help"])
    assert str(MockApplication) in result.output
