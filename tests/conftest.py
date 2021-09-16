import pytest

from isabl_cli.api import create_instance, patch_instance, get_instance
from isabl_cli.settings import _DEFAULTS


@pytest.fixture(scope="session", autouse=True)
def set_data_dir(tmpdir_factory):
    """Set data dir for all tests."""
    _DEFAULTS["BASE_STORAGE_DIRECTORY"] = tmpdir_factory.mktemp("data").strpath


@pytest.fixture(scope="session")
def use_test_client():
    """Configure cli and backend clients with extra settings."""
    client_cli_id = "test-cli-client"
    client_api_id = "default-backend-settings"
    new_settings_cli = {"EXTRA_RAW_DATA_FORMATS": [(r"\.fkf(\.gz)?$", "FAKEFORMAT")]}
    new_settings_api = {"EXTRA_RAW_DATA_FORMATS": [("FAKEFORMAT", "FAKEFORMAT")]}

    # Get or create the test cli client and patch its settings
    client_cli = create_instance("clients", slug=client_cli_id)
    client_cli = patch_instance("clients", client_cli.pk, settings=new_settings_cli)

    # Get or create the test backend client and patch its settings
    client_api = get_instance("clients", client_api_id)
    client_api = patch_instance("clients", client_api.pk, settings=new_settings_api)

    return client_cli_id
