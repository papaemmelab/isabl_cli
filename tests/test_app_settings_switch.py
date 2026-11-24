"""Pure tests for how `AbstractApplication` reads a settings manifest."""

import pytest

from isabl_cli import api
from isabl_cli.exceptions import ConfigurationError
from tests.settings_app.apps import SettingsFileApp
from tests.settings_app.apps import SettingsFileAppMissingPath

API_URL = "http://localhost:8000/api/v1/"


@pytest.fixture(autouse=True)
def no_api(monkeypatch: pytest.MonkeyPatch) -> None:
    """Point at the fixture block and fail on any API request."""

    def _fail(*args: object, **kwargs: object) -> None:
        raise AssertionError(f"unexpected API request: {args} {kwargs}")

    monkeypatch.setenv("ISABL_API_URL", API_URL)
    # Every request to the API, the token request too, goes through it.
    monkeypatch.setattr(api, "retry_request", _fail)


def test_settings_reads_the_file() -> None:
    settings = SettingsFileApp().settings

    assert settings.greeting == "hello-from-settings-yaml"
    assert settings.reference == "/tmp"
    assert settings.head_job.WG.cores == 1
    assert settings.get("submit_analyses") is None


def test_run_raises_on_a_missing_path() -> None:
    with pytest.raises(ConfigurationError) as error:
        SettingsFileAppMissingPath().run([([], [])], commit=False)

    assert "settings_app/settings.yaml [SettingsFileAppMissingPath]:" in str(
        error.value
    )
    assert "path does not exist: reference = /no/such/settings/reference" in str(
        error.value
    )
