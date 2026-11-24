"""
Compare the applications installed in this environment with their database objects.

Running an application never creates or changes its database object, it must be
registered first with `isabl register-apps`. `isabl list-apps` shows which
installed applications are registered and which differ from their database object.
"""

from collections import defaultdict
from collections import namedtuple

from isabl_cli import api
from isabl_cli.discovery import get_command_name
from isabl_cli.discovery import import_string as get_import_string
from isabl_cli.settings import import_from_string
from isabl_cli.settings import system_settings

REGISTERED = "registered"
CHANGED = "changed"
NOT_REGISTERED = "not registered"
SHARED = "shared"
CONFLICT = "conflict"
INVALID = "invalid"
NOT_INSTALLED = "not installed"

# statuses `isabl register-apps` acts on by default
NEEDS_REGISTRATION = {CHANGED, NOT_REGISTERED}

Installed = namedtuple("Installed", ["app", "import_string", "command_name", "error"])


class RegistrationStatus:

    """The registration status of an installed application or a database object."""

    def __init__(self, status, installed=None, record=None, key=None, details=()):
        """Store the status of an installed application and its database object."""
        self.status = status
        self.installed = installed
        self.record = record
        self.key = key
        self.details = list(details)

    @property
    def app(self):
        """Return the application class, None if not installed or not imported."""
        return self.installed.app if self.installed else None

    @property
    def command_name(self):
        """Return the command name, or the name of the database object."""
        if self.installed:
            return self.installed.command_name
        return self.key[0] if self.key else ""

    @property
    def import_string(self):
        """Return the import string of the installed class or the database object."""
        if self.installed:
            return self.installed.import_string
        return getattr(self.record, "application_class", None) or ""

    @property
    def assembly(self):
        """Return the assembly name."""
        return (self.key[2] if self.key else None) or ""

    @property
    def primary(self):
        """Return True if the installed class is the primary class."""
        return bool(self.app and self.app.application_primary)


def get_installed_applications():
    """
    Return the applications of this environment, the same ones the CLI lists.

    Applications provided as metadata by packages are imported here, failures are
    returned with the error instead of the class.

    Returns:
        list: of Installed tuples.
    """
    from isabl_cli import providers  # pylint: disable=C0415

    eager, lazy = providers.collect_applications(
        installed=system_settings.INSTALLED_APPLICATIONS,
        excluded=system_settings.EXCLUDED_APPLICATIONS,
        discover=system_settings.DISCOVER_APPLICATIONS,
    )

    installed = [
        Installed(i, get_import_string(i), get_command_name(i), None) for i in eager
    ]

    for i in lazy:
        try:
            app = import_from_string(i.import_string)
            installed.append(Installed(app, i.import_string, i.command_name, None))
        except Exception as error:  # pylint: disable=W0703
            installed.append(Installed(None, i.import_string, i.command_name, error))

    # a class can be found under more than one import string, e.g. re-exported
    unique, seen = [], set()
    for i in installed:
        if i.app is None or i.app not in seen:
            seen.add(i.app)
            unique.append(i)

    return unique


def get_application_records():
    """Return all application database objects by (name, version, assembly)."""
    return {
        record_key(i): i
        for i in api.get_instances("applications", verbose=False)
    }


def record_key(record):
    """Return the (name, version, assembly name) of an application database object."""
    return (record.name, record.version, record.assembly.name if record.assembly else None)


def app_key(app, name=None):
    """Return the (name, version, assembly name) of an application class."""
    return (name or app.NAME, app.VERSION, app.ASSEMBLY or None)


def get_statuses(installed, records, show_all=False):
    """
    Return the registration status of each installed application.

    Arguments:
        installed (list): of Installed tuples, see `get_installed_applications`.
        records (dict): application database objects, see `get_application_records`.
        show_all (bool): also return database objects with no installed class.

    Returns:
        list: of RegistrationStatus.
    """
    statuses, primaries, merge_keys = [], defaultdict(list), set()

    for i in installed:
        if i.app is not None and i.app.application_primary:
            primaries[app_key(i.app)].append(i)

    for i in installed:
        if i.app is None:
            details = [f"failed to import: {i.error}"]
            statuses.append(RegistrationStatus(INVALID, i, details=details))
            continue

        key = app_key(i.app)

        try:
            instance = i.app()
            merges = _get_merges(instance)
        except Exception as error:  # pylint: disable=W0703
            details = [f"invalid: {error!r}"]
            statuses.append(RegistrationStatus(INVALID, i, key=key, details=details))
            continue

        merge_keys.update(app_key(i.app, name) for name, _ in merges)
        record = records.get(key)

        if i.app.application_primary and len(primaries[key]) > 1:
            others = [j.import_string for j in primaries[key] if j is not i]
            details = [f"also the primary class: {', '.join(others)}"]
            statuses.append(RegistrationStatus(CONFLICT, i, record, key, details))
        elif record is None:
            details = (
                []
                if i.app.application_primary
                else ["registered by its primary class"]
            )
            statuses.append(RegistrationStatus(NOT_REGISTERED, i, None, key, details))
        elif i.app.application_primary:
            details = compare(instance, record, records)
            status = CHANGED if details else REGISTERED
            statuses.append(RegistrationStatus(status, i, record, key, details))
        else:
            details = _diff_results(instance._application_results, record.results)
            status = CHANGED if details else SHARED
            details.insert(0, f"uses {record.application_class}")
            statuses.append(RegistrationStatus(status, i, record, key, details))

    if show_all:
        installed_keys = {i.key for i in statuses if i.key}

        for key, record in sorted(records.items(), key=lambda i: str(i[0])):
            if key not in installed_keys and key not in merge_keys:
                statuses.append(
                    RegistrationStatus(NOT_INSTALLED, record=record, key=key)
                )

    return statuses


def compare(app, record, records):
    """
    Return how an application database object differs from the installed class.

    Arguments:
        app (AbstractApplication): an instance of the primary class.
        record (types.SimpleNamespace): its database object.
        records (dict): all database objects, to check auto merge applications.

    Returns:
        list: of differences, empty if they match.
    """
    details = []
    class_diff = _diff_class(app, record.application_class)

    if class_diff:
        details.append(class_diff)

    details += _diff_results(app._application_results, record.results)

    if (app.application_description or "") != (record.description or ""):
        details.append("description differs")

    if (app.application_url or None) != (record.url or None):
        details.append("url differs")

    for name, results in _get_merges(app):
        merge = records.get(app_key(app, name))

        if merge is None:
            details.append(f"'{name}' not registered")
        else:
            details += [f"{name}: {i}" for i in _diff_results(results, merge.results)]

    return details


def _get_merges(app):
    merges = []

    if app.has_individual_auto_merge:
        merges.append(
            (
                app.individual_level_auto_merge_name,
                app._application_individual_level_results,
            )
        )

    if app.has_project_auto_merge:
        merges.append(
            (app.project_level_auto_merge_name, app._application_project_level_results)
        )

    return merges


def _diff_class(app, stored):
    """Compare the stored application class with the class, by what it imports."""
    if not stored:
        return "no application class in the database"

    if stored == app.get_application_class_path():
        return None

    try:
        if import_from_string(stored) is type(app):
            return None
    except ImportError:
        return f"application class {stored} can't be imported"

    return f"application class is {stored}"


def _diff_results(expected, stored):
    """
    Compare the results specification of a class with the stored one.

    The API keeps only the fields it knows of each result (e.g. it drops `pattern`)
    and fills the others with None, so only the stored fields are compared, and a
    missing field equals None.
    """
    stored = stored or {}
    details = []
    added = sorted(set(expected) - set(stored))
    removed = sorted(set(stored) - set(expected))
    changed = sorted(
        key
        for key in set(expected) & set(stored)
        if any(
            (expected[key] or {}).get(field) != value
            for field, value in (stored[key] or {}).items()
        )
    )

    if added:
        details.append(f"results not in the database: {', '.join(added)}")

    if removed:
        details.append(f"results only in the database: {', '.join(removed)}")

    if changed:
        details.append(f"results that differ: {', '.join(changed)}")

    return details


def find_statuses(statuses, identifiers):
    """
    Return the statuses of installed applications matching any identifier.

    Identifiers are command names, import strings or class names.

    Returns:
        tuple: matching statuses and identifiers that matched nothing.
    """
    matched, missing = [], []

    for identifier in identifiers:
        found = [
            i
            for i in statuses
            if i.installed
            and identifier
            in {
                i.command_name,
                i.import_string,
                i.import_string.rsplit(".", 1)[-1],
            }
        ]

        if not found:
            missing.append(identifier)

        matched += [i for i in found if i not in matched]

    return matched, missing
