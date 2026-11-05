"""
Find applications in installed packages.

These helpers let an application package (e.g. isabl_apps) list its applications
without a hand-maintained INSTALLED_APPLICATIONS setting.
"""

from fnmatch import fnmatch
from pathlib import Path
import importlib
import inspect

from slugify import slugify
import click

try:
    from importlib import metadata as importlib_metadata
except ImportError:  # pragma: no cover, python < 3.8
    import importlib_metadata


def cli_command_name(name, version):
    """Return the command name of an application given its NAME and VERSION."""
    return f"{slugify(name)}-{slugify(version, separator='.')}"


def import_string(obj):
    """Return the dotted path of a class or function."""
    return f"{obj.__module__}.{obj.__qualname__}"


def get_entry_points(group):
    """Return the installed entry points of a group, on any python version."""
    entry_points = importlib_metadata.entry_points()

    if hasattr(entry_points, "select"):
        selected = entry_points.select(group=group)
    else:  # pragma: no cover, python < 3.10
        selected = entry_points.get(group, [])

    # the same distribution can be found twice in sys.path
    unique = {}
    for i in selected:
        unique.setdefault((i.name, i.value), i)

    return list(unique.values())


def is_excluded(path, command_name, patterns):
    """
    Check if an application matches any of the EXCLUDED_APPLICATIONS patterns.

    Patterns are shell-style globs matched against the application import string
    (e.g. `isabl_apps.apps.legacy.*`) and its command name (e.g. `bwa-mem-0.1.0`).
    """
    names = [i for i in (path, command_name) if i]
    return any(fnmatch(name, pattern) for pattern in patterns or [] for name in names)


def get_command_name(app):
    """
    Return the command name of an application class, or None if not an app.

    Applications can override `get_cli_command_name`, e.g. to give a finalize
    step its own command, so the name is taken from an instance. If the app
    can't be instantiated here (e.g. it needs an environment variable), the
    default name is used.
    """
    from isabl_cli.app import AbstractApplication  # pylint: disable=C0415

    if not (
        inspect.isclass(app)
        and issubclass(app, AbstractApplication)
        and app.NAME
        and app.VERSION
    ):
        return None

    try:
        return app().get_cli_command_name()
    except Exception as error:  # pylint: disable=W0703
        _warn(f"Using the default command name of {import_string(app)}: {error!r}")
        return cli_command_name(app.NAME, app.VERSION)


def find_applications(package, modules=None, required=("NAME", "VERSION")):
    """
    Import the modules of a package and return the applications they hold.

    A class is returned if it subclasses AbstractApplication, is defined inside
    `package` (applications imported from other packages are skipped), sets
    every attribute in `required` and doesn't set `_PRIVATE = True`. Modules
    that fail to import are reported and skipped.

    Arguments:
        package (str): package to search, e.g. `"isabl_apps.apps"`.
        modules (str): glob matched against module names to only import some of
            them, e.g. `"isabl_apps.apps.*.apps"`. Defaults to all modules, so
            set it if the package holds scripts that shouldn't be imported.
        required (tuple): attributes that must be set, add `"ASSEMBLY"` to skip
            base classes that are only specialized by assembly in subclasses.

    Returns:
        list: application classes, sorted by the module they were found in.
    """
    return list(_find_applications(package, modules, required).values())


def find_application_specs(package, modules=None, required=("NAME", "VERSION")):
    """
    Find applications like `find_applications` and return their metadata.

    Import strings point to the module where each application was found (e.g.
    an `apps.py` that imports it) instead of the module that defines it. That
    is the module the package meant to be imported, importing a defining module
    first can fail on circular imports.

    Returns:
        list: dicts made with `application_spec`.
    """
    found = _find_applications(package, modules, required)
    return [application_spec(app, path) for path, app in found.items()]


def _find_applications(package, modules, required):
    """Return a dict of import string (where found) to application class."""
    root = importlib.import_module(package)
    names = [root.__name__] + _get_module_names(root)

    if modules:
        names = [i for i in names if fnmatch(i, modules)]

    found, seen = {}, set()
    for module in _import_modules(names):
        for name, obj in vars(module).items():
            if _is_application(obj, root.__name__, required) and obj not in seen:
                seen.add(obj)
                found[f"{module.__name__}.{name}"] = obj

    return {i: found[i] for i in sorted(found)}


def application_spec(app, path=None):
    """
    Return the metadata needed to list an application without importing it.

    Application packages can save these dictionaries (e.g. in a JSON file
    generated at release time) and serve them to Isabl CLI.

    Arguments:
        app (class): an application.
        path (str): import string to use, defaults to where `app` is defined.
    """
    help_text = getattr(app, "cli_help", "")

    return {
        "import_string": path or import_string(app),
        "name": app.NAME,
        "version": app.VERSION,
        "assembly": app.ASSEMBLY,
        "command_name": get_command_name(app),
        "help": help_text if isinstance(help_text, str) else "",
    }


def _import_modules(names):
    """
    Import modules, retrying failures while others succeed.

    Circular imports can fail or not depending on what was imported first, so a
    module that fails is tried again after the rest. Scripts that call exit at
    import are skipped too.
    """
    modules, pending, errors = [], list(names), {}

    while pending:
        failed = []

        for name in pending:
            try:
                modules.append(importlib.import_module(name))
            except (Exception, SystemExit) as error:  # pylint: disable=W0703
                errors[name] = error
                failed.append(name)

        if len(failed) == len(pending):
            break

        pending = failed

    for name in pending:
        _warn(f"Skipping {name}, failed to import it: {errors[name]!r}")

    return modules


def _get_module_names(package):
    """List the modules inside a package without importing them."""
    names = set()

    for path in getattr(package, "__path__", []):
        for i in Path(path).rglob("*.py"):
            parts = i.relative_to(path).with_suffix("").parts

            if parts[-1] == "__init__":
                parts = parts[:-1]

            if parts and all(j.isidentifier() for j in parts):
                names.add(".".join((package.__name__,) + parts))

    return sorted(names)


def _is_application(obj, package, required):
    from isabl_cli.app import AbstractApplication  # pylint: disable=C0415

    return (
        inspect.isclass(obj)
        and issubclass(obj, AbstractApplication)
        and _is_defined_in(obj, package)
        and all(getattr(obj, i, None) for i in required)
        and not getattr(obj, "_PRIVATE", False)
    )


def _is_defined_in(obj, package):
    return obj.__module__ == package or obj.__module__.startswith(f"{package}.")


def _warn(msg):
    click.secho(msg, err=True, fg="yellow")
