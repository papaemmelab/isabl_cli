import pytest

from isabl_cli.settings import _DEFAULTS


@pytest.fixture(scope="session", autouse=True)
def set_data_dir(tmpdir_factory):
    """Set data dir for all tests."""
    _DEFAULTS["BASE_STORAGE_DIRECTORY"] = tmpdir_factory.mktemp("data").strpath
