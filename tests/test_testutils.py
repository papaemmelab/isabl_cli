from isabl_cli import api
from isabl_cli import factories
from isabl_cli.test import utils

from .test_app import TestApplication


def test_create_utils():
    # for now just get some coverage
    assert utils.create_experiment("/a/bam")
    assert utils.create_pair("/tumor/bam", "/tumor/bam")
    assert utils.create_test_result()


def test_assert_run():
    # test assert_run utility
    data = api.create_instance("projects", **factories.ProjectFactory())
    data = factories.ExperimentFactory(projects=[data])
    data["sample"]["individual"]["species"] = "HUMAN"
    assert utils.assert_run(
        application=TestApplication(),
        tuples=[([api.create_instance("experiments", **data)], [])],
        commit=True,
    )
