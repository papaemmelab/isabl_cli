"""
Test configuration.

The tests create, change and delete records and replace the settings of the
`default-backend-settings` client, which configures the whole API. They must only run
against a disposable local API, see tests/test-with-compose.sh. pytest.ini sets
ISABL_API_URL to it, and the check below stops the tests before isabl_cli is imported
if the URL is not local anyway (e.g. pytest-env is not installed).
"""

from urllib.parse import urlparse
import os

import pytest

LOCAL_API_HOSTS = {"localhost", "127.0.0.1", "::1"}


def check_test_api_url():
    """Raise a UsageError unless the tests would use a local API."""
    url = os.environ.get("ISABL_API_URL", "http://localhost:8000/api/v1/")

    if urlparse(url).hostname not in LOCAL_API_HOSTS:
        raise pytest.UsageError(
            f"Refusing to run the tests against {url}, they create, change and "
            "delete records and replace client settings. Run them against a local "
            "API, see tests/test-with-compose.sh."
        )


check_test_api_url()

# pylint: disable=wrong-import-position
from isabl_cli.api import create_instance, patch_instance, get_instance
from isabl_cli.settings import _DEFAULTS


@pytest.hookimpl(tryfirst=True)
def pytest_runtest_setup(item):
    """Check the API URL again before each test and its fixtures are set up."""
    check_test_api_url()


def pytest_sessionstart():
    """Run before any test starts."""
    check_test_api_url()  # again, right before the clients are changed
    client_cli_id = "test-cli-client"
    client_api_id = "default-backend-settings"
    new_settings_cli = {"EXTRA_RAW_DATA_FORMATS": [(r"\.tf(\.gz)?$", "TEST_FORMAT")]}
    new_settings_api = {"EXTRA_RAW_DATA_FORMATS": [("TEST_FORMAT", "TEST_FORMAT")]}

    # Get or create the test cli client and patch its settings
    client_cli = create_instance("clients", slug=client_cli_id)
    client_cli = patch_instance("clients", client_cli.pk, settings=new_settings_cli)

    # Get or create the test backend client and patch its settings
    client_api = get_instance("clients", client_api_id)
    client_api = patch_instance("clients", client_api.pk, settings=new_settings_api)


@pytest.fixture(scope="session", autouse=True)
def set_data_dir(tmpdir_factory):
    """Set data dir for all tests."""
    _DEFAULTS["BASE_STORAGE_DIRECTORY"] = tmpdir_factory.mktemp("data").strpath

