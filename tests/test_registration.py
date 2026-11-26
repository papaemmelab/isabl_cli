"""Tests for registering applications and comparing them with the database."""

from types import SimpleNamespace
import uuid

from click.testing import CliRunner
import pytest

from isabl_cli import AbstractApplication
from isabl_cli.app import strip_results_whitespace
from isabl_cli.app import strip_whitespace
from isabl_cli import api
from isabl_cli import commands
from isabl_cli import exceptions
from isabl_cli import registration

pytestmark = pytest.mark.no_auto_register

UNIQUE = str(uuid.uuid4())


class Primary(AbstractApplication):
    NAME = f"REGISTRATION {UNIQUE}"
    VERSION = "1.0.0"
    ASSEMBLY = "GRCh4000"
    SPECIES = "HUMAN"
    application_description = "Primary application."
    application_results = {
        "vcf": {
            "frontend_type": "text-file",
            "description": "A vcf.",
            "verbose_name": "VCF",
            "pattern": "*.vcf",
        }
    }

    def merge_project_analyses(self, analysis, analyses):  # pragma: no cover
        pass


class Step(Primary):
    application_primary = False

    def get_cli_command_name(self):
        return super().get_cli_command_name() + "_step"


class OtherPrimary(Primary):
    pass


class Unregistered(AbstractApplication):
    NAME = f"UNREGISTERED {UNIQUE}"
    VERSION = "1.0.0"


def installed(*apps):
    """Return Installed tuples for application classes."""
    return [
        registration.Installed(
            i, i.get_application_class_path(), i().get_cli_command_name(), None
        )
        for i in apps
    ]


def record(app, **overrides):
    """Return a fake database object matching an application class."""
    instance = app()
    data = dict(
        name=app.NAME,
        version=app.VERSION,
        assembly=SimpleNamespace(name=app.ASSEMBLY) if app.ASSEMBLY else None,
        application_class=app.get_application_class_path(),
        description=app.application_description,
        url=None,
        # the API drops unknown fields such as pattern and fills defaults with None
        results={
            key: {
                "frontend_type": value.get("frontend_type"),
                "description": value.get("description"),
                "verbose_name": value.get("verbose_name"),
                "external_link": None,
                "optional": None,
            }
            for key, value in instance._application_results.items()
        },
    )
    data.update(overrides)
    return SimpleNamespace(**data)


def records_of(*items):
    return {registration.record_key(i): i for i in items}


def merge_record(app):
    return record(
        app,
        name=app().project_level_auto_merge_name,
        results={i: {} for i in app()._application_project_level_results},
    )


def by_class(statuses):
    return {i.import_string.rsplit(".", 1)[-1]: (i.status, i.details) for i in statuses}


def test_statuses_registered_and_shared():
    records = records_of(record(Primary), merge_record(Primary))
    statuses = registration.get_statuses(installed(Primary, Step), records)
    assert by_class(statuses) == {
        "Primary": ("registered", []),
        "Step": ("shared", [f"uses {Primary.get_application_class_path()}"]),
    }


def test_statuses_not_registered():
    statuses = registration.get_statuses(installed(Primary, Step, Unregistered), {})
    assert by_class(statuses) == {
        "Primary": ("not registered", []),
        "Step": ("not registered", ["registered by its primary class"]),
        "Unregistered": ("not registered", []),
    }


def test_statuses_changed():
    stored = record(
        Primary,
        application_class=Step.get_application_class_path(),
        description="old",
        results={
            "vcf": {"description": "An old vcf.", "frontend_type": "text-file"},
            "old": {},
        },
    )
    statuses = registration.get_statuses(installed(Primary), records_of(stored))
    status, details = by_class(statuses)["Primary"]
    assert status == "changed"
    assert details == [
        f"application class is {Step.get_application_class_path()}",
        "results not in the database: command_err, command_log, command_script",
        "results only in the database: old",
        "results that differ: vcf",
        "description differs",
        f"'{Primary().project_level_auto_merge_name}' not registered",
    ]


class Padded(Primary):
    application_description = """
    Primary application, written as a triple-quoted string.
    """
    application_url = ""
    application_results = {
        "vcf": {"frontend_type": "text-file", "description": "\nA vcf.\n"}
    }


def test_strip_results_whitespace_only_strips_text():
    results = {
        "vcf": {"description": " A vcf.\n", "pattern": " *.vcf ", "optional": True}
    }
    assert strip_results_whitespace(results) == {
        "vcf": {"description": "A vcf.", "pattern": " *.vcf ", "optional": True}
    }
    assert strip_whitespace("\n Text \n") == "Text"
    assert strip_whitespace(None) is None


def test_statuses_compare_text_as_stored():
    # the API trims text, and register sends it trimmed too
    stored = record(
        Padded,
        description="Primary application, written as a triple-quoted string.",
        url=None,
    )
    stored.results["vcf"]["description"] = "A vcf."
    records = records_of(stored, merge_record(Padded))
    assert by_class(registration.get_statuses(installed(Padded), records)) == {
        "Padded": ("registered", [])
    }

    stored.description = "Another description."
    status, details = by_class(registration.get_statuses(installed(Padded), records))[
        "Padded"
    ]
    assert (status, details) == ("changed", ["description differs"])


def test_statuses_class_path_compared_by_import():
    # a different import string of the same class is not a change
    stored = record(Primary, application_class="tests.test_registration.Primary")
    records = records_of(stored, merge_record(Primary))
    assert by_class(registration.get_statuses(installed(Primary), records)) == {
        "Primary": ("registered", [])
    }

    stored = record(Primary, application_class="tests.nowhere.Primary")
    records = records_of(stored, merge_record(Primary))
    status, details = by_class(registration.get_statuses(installed(Primary), records))[
        "Primary"
    ]
    assert status == "changed"
    assert details == ["application class tests.nowhere.Primary can't be imported"]


def test_statuses_conflict_invalid_and_not_installed():
    broken = registration.Installed(None, "nowhere.Broken", "broken-1.0", "boom")
    other = record(Unregistered)
    statuses = registration.get_statuses(
        installed(Primary, OtherPrimary) + [broken],
        records_of(record(Primary), merge_record(Primary), other),
        show_all=True,
    )
    by_status = {(i.status, i.command_name) for i in statuses}
    assert ("conflict", Primary().get_cli_command_name()) in by_status
    assert ("invalid", "broken-1.0") in by_status
    assert ("not installed", Unregistered.NAME) in by_status
    # the merge application of an installed app is not listed as not installed
    assert all(i.status != "not installed" or i.key[0] == other.name for i in statuses)


def test_find_statuses():
    statuses = registration.get_statuses(installed(Primary, Step), {})
    command = Step().get_cli_command_name()
    found, missing = registration.find_statuses(statuses, [command, "Primary", "nope"])
    assert [i.app for i in found] == [Step, Primary]
    assert missing == ["nope"]


def test_application_is_not_created_when_run():
    app = Unregistered()

    with pytest.raises(exceptions.ApplicationNotRegisteredError) as error:
        app.application

    assert "isabl register-apps" in str(error.value)
    assert app.get_registered_application() is None


def test_register_and_running_doesnt_change_it():
    app = Primary()
    application = app.register()
    assert application.application_class == Primary.get_application_class_path()
    assert application.description == Primary.application_description
    assert set(application.results) == {
        "vcf",
        "command_err",
        "command_log",
        "command_script",
    }
    assert application.settings == {}  # register never writes settings

    merge = app.project_level_auto_merge_application
    assert merge.application_class == Primary.get_application_class_path()
    assert merge.name == f"{Primary.NAME} Project Application"

    # changes in the database are kept when an app is used
    api.patch_instance("applications", application.pk, description="edited")
    fresh = Primary()
    assert fresh.application.description == "edited"
    assert fresh.settings is not None
    assert api.get_instance("applications", application.pk).description == "edited"

    # steps use the application of the primary class
    assert Step().application.pk == application.pk

    # registering again updates it, and keeps the settings
    settings = {"7": {"reference": "/a/path"}}
    api.patch_instance("applications", application.pk, settings=settings)
    application = Primary().register()
    assert application.description == Primary.application_description
    assert application.settings == settings


def test_register_step_fails():
    with pytest.raises(exceptions.ConfigurationError, match="primary class"):
        Step().register()


class Fresh(AbstractApplication):
    NAME = f"FRESH {UNIQUE}"
    VERSION = "1.0.0"


def test_register_apps_command(monkeypatch):
    monkeypatch.setattr(
        registration, "get_installed_applications", lambda: installed(Fresh, Step)
    )
    runner = CliRunner()
    command = Fresh().get_cli_command_name()

    result = runner.invoke(commands.register_apps, [command])
    assert result.exit_code == 0, result.output
    assert "Would register" in result.output
    assert Fresh().get_registered_application() is None

    result = runner.invoke(commands.register_apps, [command, "--commit"])
    assert result.exit_code == 0, result.output
    assert Fresh().get_registered_application() is not None

    result = runner.invoke(commands.list_apps, ["--status", "registered"])
    assert result.exit_code == 0, result.output
    assert command in result.output

    result = runner.invoke(commands.register_apps, [Step().get_cli_command_name()])
    assert "not the primary class" in result.output

    result = runner.invoke(commands.register_apps, ["nope"])
    assert result.exit_code == 2
    assert "No installed application matches: nope" in result.output


def test_class_mismatch_rules():
    primary = Primary.get_application_class_path()
    step = Step.get_application_class_path()
    other = OtherPrimary.get_application_class_path()
    foreign = type(
        "Foreign",
        (Primary,),
        {"application_primary": False, "__module__": "other.apps"},
    )

    # the primary class only uses an application registered to itself
    assert Primary().get_class_mismatch(primary) is None
    assert "registered to" in OtherPrimary().get_class_mismatch(primary)
    assert "can't be imported" in Primary().get_class_mismatch("nowhere.Primary")
    assert "no registered class" in Primary().get_class_mismatch(None)

    # records that point to a step or alias are fixed by registering the primary
    assert "register-apps <its command>" in Primary().get_class_mismatch(step)

    # steps and aliases use the application of a primary class of their package
    assert Step().get_class_mismatch(primary) is None
    assert Step().get_class_mismatch(other) is None
    assert Step().get_class_mismatch(step)
    assert foreign().get_class_mismatch(primary)


def test_registration_conflict_rules():
    def record(path):
        return SimpleNamespace(application_class=path)

    primary = Primary()
    own = Primary.get_application_class_path()
    step = Step.get_application_class_path()
    assert primary.get_registration_conflict(None) is None
    assert primary.get_registration_conflict(record(None)) is None
    assert primary.get_registration_conflict(record(own)) is None

    # a record written by a step or alias of this package is corrected
    assert primary.get_registration_conflict(record(step)) is None

    # another class owns it
    other = OtherPrimary.get_application_class_path()
    assert other in primary.get_registration_conflict(record(other))
    assert "can't be imported" in primary.get_registration_conflict(record("nowhere.X"))


class Owner(AbstractApplication):
    NAME = f"OWNERSHIP {UNIQUE}"
    VERSION = "1.0.0"


class Intruder(Owner):
    pass


def test_register_and_run_only_as_the_owner(monkeypatch):
    Owner().register()

    # another class can't register or use the application
    with pytest.raises(exceptions.ApplicationClassMismatchError, match="--takeover"):
        Intruder().register()
    with pytest.raises(exceptions.ApplicationClassMismatchError, match="registered to"):
        Intruder().application

    # the command skips it unless --takeover is passed
    monkeypatch.setattr(
        registration, "get_installed_applications", lambda: installed(Intruder)
    )
    command = Intruder().get_cli_command_name()
    runner = CliRunner()
    result = runner.invoke(commands.register_apps, [command, "--commit"])
    assert result.exit_code == 0, result.output
    assert "Use --takeover to replace it" in result.output
    assert Owner().application.application_class == Owner.get_application_class_path()

    result = runner.invoke(commands.register_apps, [command, "--takeover", "--commit"])
    assert result.exit_code == 0, result.output
    assert "Registered 1 applications" in result.output
    path = Intruder.get_application_class_path()
    assert Intruder().application.application_class == path

    # now the previous owner can't use it, until it takes it back
    with pytest.raises(exceptions.ApplicationClassMismatchError):
        Owner().application
    application = Owner().register(takeover=True)
    assert application.application_class == Owner.get_application_class_path()


class PaddedText(AbstractApplication):
    NAME = f"PADDED TEXT {UNIQUE}"
    VERSION = "1.0.0"
    application_description = """
    Written as a triple-quoted string.
    """
    application_results = {
        "vcf": {"frontend_type": "text-file", "description": "\nA vcf.\n"}
    }


def test_register_sends_text_as_stored():
    application = PaddedText().register()
    assert application.description == "Written as a triple-quoted string."
    assert application.results["vcf"]["description"] == "A vcf."

    statuses = registration.get_statuses(
        installed(PaddedText), registration.get_application_records()
    )
    assert by_class(statuses) == {"PaddedText": ("registered", [])}
