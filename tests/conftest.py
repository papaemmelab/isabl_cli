import pytest

from isabl_cli.api import create_instance, patch_instance
from isabl_cli.settings import _DEFAULTS


@pytest.fixture(scope="session", autouse=True)
def set_data_dir(tmpdir_factory):
    """Set data dir for all tests."""
    _DEFAULTS["BASE_STORAGE_DIRECTORY"] = tmpdir_factory.mktemp("data").strpath


@pytest.fixture
def use_test_client(monkeypatch):
    """Create client settings for testing."""
    client_id = "test-client"

    # Get or create the test client and patch its settings
    client = create_instance("clients", slug=client_id)
    client = patch_instance(
        "clients",
        client.pk,
        settings={"EXTRA_RAW_DATA_FORMATS": [(r"\.fkf(\.gz)?$", "FAKEFORMAT")]},
    )
    # Set the test environment variable
    monkeypatch.setenv("ISABL_CLIENT_ID", client_id)
    return client_id
