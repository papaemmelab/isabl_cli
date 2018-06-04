from click.testing import CliRunner
import click

from cli import _DEFAULTS
from cli import api
from cli.engine import AbstractPipeline

from . import factories

class TestPipeline(AbstractPipeline):

    NAME = 'A Test Pipeline'
    VERSION = 'A Test Version'

    @staticmethod
    def get_cli_help():
        return "This is a test pipeline"

    @staticmethod
    def get_cli_options():
        return [click.option('--filters', type=(str, str), multiple=True)]

    @staticmethod
    def get_tuples(filters, **cli_options):
        tuples = []
        filters = dict(filters)

        for i in api.get_instances('workflows', **filters):
            tuples.append(([i], [], []))

        return tuples

    def validate_tuple(self, targets, references, analyses):
        self.validate_onetarget_noreferences(targets, references, analyses)
        return True

    @staticmethod
    def get_command(analysis):
        return (
            f"echo {analysis['targets'][0]['system_id']} > "
            f"{analysis['storage_url']}/head_job.log")


def test_engine(tmpdir):
    data_storage_directory = tmpdir.mkdir('data_storage_directory')
    _DEFAULTS['BASE_STORAGE_DIRECTORY'] = data_storage_directory.strpath

    # keys = [1, 2, 3]
    # workflows = api.get_instances('workflows', keys)

    # if not workflows:
    workflows = [factories.WorkflowFactory() for i in range(3)]
    workflows = [api.create_instance('workflows', **i) for i in workflows]
    tuples = [([i], [], []) for i in workflows]
    pipeline = TestPipeline()

    # runner = CliRunner()
    # command = pipeline.get_cli_command()
    # result = runner.invoke(command, ["--help"])

    # assert "This is a test pipeline" in result.output
    # assert "--commit" in result.output
    # assert "--force" in result.output
    # assert "--verbose" in result.output

    pipeline.run_tuples(tuples, commit=True)

    # result = runner.invoke(
    #     command, ["--filters", "pk__in", "1,2,3", "--verbose"],
    #     catch_exceptions=False)

    # print('CARLOS', result.output)
    # print('CARLOS', result.exit_code)
    # print('CARLOS', result.exc_info)
    # print('CARLOS', result.exception)
