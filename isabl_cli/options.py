"""Click options decorators."""

import click

from isabl_cli import api
from isabl_cli import validators


def _get_experiments(filter_tuples):
    return api.get_instances("experiments", **dict(filter_tuples))


NULLABLE_IDENTIFIERS = click.argument("identifiers", nargs=-1, required=False)

IDENTIFIER = click.option(
    "--identifier",
    "-id",
    show_default=True,
    help="identifier to be used (traverse with dot, e.g. `sample.system_id`)",
    callback=lambda _, __, i: i.split("."),
    required=True,
)

FIELDS = click.option(
    "--field",
    "-f",
    show_default=True,
    multiple=True,
    help="fields to be retrieved (traverse with dot, e.g. `sample.disease`)",
    callback=lambda _, __, i: [j.split(".") for j in i],
    required=False,
)

ANALYSIS_STATUS = click.option(
    "--status",
    show_default=True,
    help="analysis status",
    required=True,
    type=click.Choice(
        [
            "CREATED",
            "FAILED",
            "FINISHED",
            "IN_PROGRESS",
            "STAGED",
            "STARTED",
            "SUBMITTED",
            "SUCCEEDED",
        ]
    ),
)

ENDPOINT = click.argument(
    "endpoint",
    required=True,
    type=click.Choice(
        [
            "experiments",
            "analyses",
            "projects",
            "techniques",
            "samples",
            "individuals",
            "assemblies",
            "diseases",
            "centers",
            "platforms",
            "submissions",
            "applications",
            "groups",
        ]
    ),
)

ANALYSIS_PRIMARY_KEY = click.option(
    "--key",
    help="analysis primary key",
    show_default=True,
    type=click.INT,
    required=True,
)

TECHNIQUE_IDENTIFIER = click.option(
    "--technique",
    help="technique identifier",
    show_default=True,
    type=click.STRING,
    required=True,
)

BED_TYPE = click.option(
    "--bed-type",
    help="technique BED file type.",
    show_default=True,
    default="targets",
    type=click.Choice(["targets", "baits"]),
)

COMMIT = click.option(
    "--commit", help="commit results", show_default=True, is_flag=True
)

NO_HEADERS = click.option(
    "--no-headers", help="do not output headers", show_default=True, is_flag=True
)

FILE_PATTERN = click.option(
    "--pattern",
    help="paths pattern to be matched",
    default=None,
    show_default=True,
    type=click.STRING,
    required=False,
)

FILTERS = click.option(
    "--filters",
    "-fi",
    multiple=True,
    help="API filters",
    show_default=True,
    type=(str, str),
    callback=lambda _, __, i: dict(i),
    required=True,
)

NULLABLE_FILTERS = click.option(
    "--filters",
    "-fi",
    multiple=True,
    help="API filters",
    show_default=True,
    type=(str, str),
    callback=lambda _, __, i: dict(i),
)

ANALYSES = click.option(
    "--analyses-filters",
    "-fi",
    multiple=True,
    help="API filters analyses instances",
    show_default=True,
    type=(str, str),
    callback=lambda _, __, i: api.get_instances("analyses", **dict(i)),
    required=True,
)

NORMAL_TARGETS = click.option(
    "--normal-targets-filters",
    "-fi",
    "targets",
    multiple=True,
    help="API filters for NORMAL target experiments",
    show_default=True,
    type=(str, str),
    callback=lambda _, __, i: _get_experiments(i + ("sample__sample_class", "NORMAL")),
    required=True,
)

TARGETS = click.option(
    "--targets-filters",
    "-fi",
    "targets",
    multiple=True,
    help="API filters for target experiments",
    show_default=True,
    type=(str, str),
    callback=lambda _, __, i: _get_experiments(i),
    required=True,
)

REFERENCES = click.option(
    "--references-filters",
    "-rfi",
    "references",
    multiple=True,
    help="API filters for references experiments",
    show_default=True,
    type=(str, str),
    callback=lambda _, __, i: _get_experiments(i),
    required=True,
)

NULLABLE_REFERENCES = click.option(
    "--references-filters",
    "-rfi",
    "references",
    multiple=True,
    help="API filters for references experiments",
    show_default=True,
    type=(str, str),
    callback=lambda _, __, i: _get_experiments(i) if i else [],
    required=False,
)

PAIR = click.option(
    "--pair",
    "-p",
    show_default=True,
    type=(str, str),
    callback=lambda _, __, i: validators.validate_pairs([i])[0],
    help="Pass one tumor normal pair identifiers (e.g. -tn 1 2).",
)

PAIRS = click.option(
    "--pairs",
    "-p",
    show_default=True,
    type=(str, str),
    multiple=True,
    callback=lambda _, __, i: validators.validate_pairs(i),
    help="Pass one or more tumor normal pairs (e.g. -tn 1 2 -tn 3 4).",
)

PAIRS_FROM_FILE = click.option(
    "--pairs-from-file",
    "-pf",
    show_default=True,
    callback=validators.validate_pairs_from_file,
    help="Tab separated file with 2 columns: tumor, normal system ids.",
    type=click.Path(
        exists=True, file_okay=True, dir_okay=False, writable=False, readable=True
    ),
)

DIRECTORIES = click.option(
    "--directories",
    "-di",
    show_default=True,
    multiple=True,
    help="directories to be processed",
    type=click.Path(
        exists=True,
        file_okay=False,
        dir_okay=True,
        writable=False,
        readable=True,
        resolve_path=True,
    ),
)

SYMLINK = click.option("--symlink", help="use symlink", show_default=True, is_flag=True)

REFERENCE_DATA_SOURCE = click.option(
    "--data-src",
    show_default=True,
    type=click.Path(resolve_path=True, exists=True, readable=True),
    help="path to reference data",
)

TARGETS_PATH = click.option(
    "--targets-path",
    show_default=True,
    help="path to targets bedfile",
    type=click.Path(
        resolve_path=True,
        exists=True,
        file_okay=True,
        dir_okay=False,
        writable=False,
        readable=True,
    ),
)

BAITS_PATH = click.option(
    "--baits-path",
    show_default=True,
    help="path to baits bedfile",
    type=click.Path(
        resolve_path=True,
        exists=True,
        file_okay=True,
        dir_okay=False,
        writable=False,
        readable=True,
    ),
)

FILES_DATA = click.option(
    "--files-data",
    help="a yaml file with extra annotation for imported files.",
    show_default=True,
    default=None,
    type=click.Path(
        resolve_path=True,
        exists=True,
        file_okay=True,
        dir_okay=False,
        writable=False,
        readable=True,
    ),
)

QUIET = click.option(
    "--quiet", help="Don't display verbose output", show_default=False, is_flag=True
)

VERBOSE = click.option(
    "--verbose", help="Display verbose output", show_default=True, is_flag=True
)


def get_analyses_filters_option(application_classes=None, **defaults):
    """Get analyses filters with `defaults`."""
    msg = ""
    application_classes = application_classes or []

    if defaults:
        msg += " - the following filters will also be included: "
        msg += ", ".join(f"{i}={j}" for i, j in defaults.items())

    if application_classes:
        msg += " - analyses will be limited to the following applications: "
        msg += ", ".join(map(str, application_classes))

    def callback(tuples):
        if application_classes:
            analyses = []
            for i in application_classes:
                defaults["application__pk"] = str(i.primary_key)
                click.secho(
                    f"Getting completed analyses for {i.NAME} {i.VERSION}:", fg="cyan"
                )
                analyses += api.get_instances(
                    "analyses", **{**dict(tuples), **defaults}
                )
            return analyses
        return api.get_instances("analyses", **{**dict(tuples), **defaults})

    return click.option(
        "--analyses-filters",
        "-fi",
        multiple=True,
        help="API filters for analyses instances" + msg,
        show_default=True,
        type=(str, str),
        callback=lambda _, __, i: callback(i),
    )
