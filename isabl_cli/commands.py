"""commands logic."""

from collections import OrderedDict
from glob import glob
from os.path import join
from requests.exceptions import HTTPError
from subprocess import CalledProcessError
import json
import os
import shutil
import subprocess
import tempfile

import click

from isabl_cli import api
from isabl_cli import exceptions
from isabl_cli import options
from isabl_cli import utils
from isabl_cli.settings import import_from_string
from isabl_cli.settings import user_settings


def cb_app_results_keys(ctx, param, value):
    """Print applications results keys."""
    if not value or ctx.resilient_parsing:  # pragma: no cover
        return

    click.echo(
        "\n".join(
            f"{click.style(i, fg='green')}\t{j.description}"
            for i, j in sorted(api.get_instance("applications", value).results.items())
        ).expandtabs(30)
    )
    ctx.exit()


def _filters_or_identifiers(endpoint, identifiers, filters, fields=None):
    if filters and identifiers:  # pragma: no cover
        raise click.UsageError("Can't combine filters and identifiers.")

    if fields:
        filters["fields"] = fields
        filters["limit"] = 100_000

    return (
        [api.get_instance(endpoint, i, fields=fields) for i in identifiers]
        if identifiers
        else api.get_instances(endpoint, verbose=True, **filters)
    )


@click.command()
def login():  # pragma: no cover
    """Login with isabl credentials."""
    user_settings.api_token = None
    api.get_token_headers.cache_clear()
    api.get_token_headers()


@click.command()
@click.option("--project", help="primary key of project to merge by", type=int)
@click.option("--application", help="analyses application primary key", type=int)
def merge_project_analyses(project, application):  # pragma: no cover
    """Merge analyses by project."""
    project = api.get_instance("projects", project)
    application = api.get_instance("applications", application)
    application = import_from_string(application["application_class"])()
    application.run_project_merge(project)


@click.command()
@click.option("--individual", help="primary key of individual to merge by", type=int)
@click.option("--application", help="analyses application primary key", type=int)
def merge_individual_analyses(individual, application):  # pragma: no cover
    """Merge analyses by individual."""
    individual = api.get_instance("individuals", individual)
    application = api.get_instance("applications", application)
    application = import_from_string(application["application_class"])()
    application.run_individual_merge(individual)


@click.command()
@click.option("--force", help="Update previously patched results.", is_flag=True)
@options.NULLABLE_FILTERS
def process_finished(filters, force):
    """Process and update finished analyses."""
    utils.check_admin()
    filters.update(status="FINISHED", fields="pk")
    tag = "PROCESSING FINISHED"

    # refetch analysis to avoid race conditions
    n_not_patched = 0
    n_patched = 0
    for i in api.get_instances("analyses", verbose=True, **filters):
        i = api.Analysis(i.pk)

        if i.status == "FINISHED":
            if not force and tag in {j.name for j in i.tags}:
                n_not_patched += 1
            else:
                api.patch_instance("analyses", i.pk, tags=i.tags + [{"name": tag}])
                try:
                    api.patch_analysis_status(i, "SUCCEEDED")
                    n_patched += 1
                    patch_error = None
                except (PermissionError, AssertionError, CalledProcessError) as error:
                    patch_error = error

                api.patch_instance("analyses", i.pk, tags=i.tags)

                if patch_error:
                    raise Exception(patch_error)
    if n_not_patched:
        click.echo(
            f"Successfully processed {n_patched} analyses, but skipped "
            + f"{n_not_patched} analyses in use by other processes"
        )


@click.command()
@options.FILTERS
@click.option("--force", help="Update previously patched results.", is_flag=True)
def patch_results(filters, force):
    """Update the results field of many analyses."""
    utils.check_admin()
    skipped = []

    with click.progressbar(
        api.get_instances("analyses", verbose=True, **filters),
        label="Patching analyses...",
    ) as bar:
        for i in bar:
            if force or not i.results:
                results = api._get_analysis_results(i, raise_error=False)
                api.patch_instance("analyses", i.pk, results=results)
            else:  # pragma: no cover
                skipped.append(i)

    if skipped:  # pragma: no cover
        click.echo(f"{len(skipped)} analyses had results, use --force to update...")


@click.command(hidden=True)
@options.ANALYSIS_PRIMARY_KEY
@options.ANALYSIS_STATUS
def patch_status(key, status):
    """Patch status of a given analysis."""
    analysis = api.get_instance("analyses", key)
    api.patch_analysis_status(analysis, status)


@click.command(
    epilog="Learn more about fx: "
    "https://github.com/antonmedv/fx/blob/master/doc/doc.md#interactive-mode"
)
@options.ENDPOINT
@options.FIELDS
@options.NULLABLE_FILTERS
@options.NO_HEADERS
@options.NULLABLE_IDENTIFIERS
@click.option("--json", "json_", help="Print as JSON", is_flag=True)
@click.option("--fx", "use_fx", help="Visualize json with fx", is_flag=True)
@click.option("--pretty", "pretty", help="prettified output", is_flag=True)
@click.option("--all", "output_all", help="include all fields in the tabular output", is_flag=True)
def get_metadata(identifiers, endpoint, field, filters, no_headers, json_, use_fx, output_all, pretty):
    """Retrieve metadata for multiple instances."""
    if not field and not (json_ or use_fx or output_all):  # pragma: no cover
        raise click.UsageError("Pass --field or use --json/--fx/--all")

    if use_fx and not shutil.which("fx"):  # pragma: no cover
        raise click.UsageError("fx is not installed")

    instances = _filters_or_identifiers(
        endpoint=endpoint,
        identifiers=identifiers,
        filters=filters,
        fields=",".join([i[0] for i in field]),
    )

    if field:  # if fields were passed, update the instances list
        instances = [
            OrderedDict([(".".join(j), utils.traverse_dict(i, j)) for j in field])
            for i in instances
        ]
    elif output_all:
        field = [[i] for i in instances[0].__dict__.keys()]

    if json_:
        click.echo(json.dumps(instances, sort_keys=True, indent=4))
    elif use_fx:  # pragma: no cover
        fp = tempfile.NamedTemporaryFile("w+", delete=False)
        json.dump(instances, fp)
        fp.close()
        subprocess.check_call(["fx", fp.name])
        os.unlink(fp.name)
    else:
        result = [] if no_headers else ["\t".join(".".join(i) for i in field)]
        result += ["\t".join(map(str, i.values())) for i in instances]
        if pretty:
            click.echo("\n".join(result).expandtabs(30))
        else:
            click.echo("\n".join(result))


@click.command()
@options.ENDPOINT
@options.NULLABLE_FILTERS
def get_count(endpoint, filters):
    """Get count of database instances."""
    click.echo(api.get_instances_count(endpoint, **filters))


@click.command()
@options.ENDPOINT
@options.FILE_PATTERN
@options.NULLABLE_FILTERS
@options.NULLABLE_IDENTIFIERS
def get_paths(endpoint, pattern, filters, identifiers):
    """Get storage directories, use `pattern` to match files inside the directories."""
    for i in _filters_or_identifiers(
        endpoint=endpoint,
        identifiers=identifiers,
        filters=filters,
        fields="storage_url",
    ):
        if i["storage_url"]:
            if pattern:
                click.echo("\n".join(glob(join(i["storage_url"], pattern))))
            else:
                click.echo(i["storage_url"])


@click.command()
@options.FILE_PATTERN
@options.NULLABLE_FILTERS
@options.NULLABLE_IDENTIFIERS
def get_outdirs(pattern, filters, identifiers):
    """Get analyses outdirs, use `pattern` to match files inside the directories."""
    for i in _filters_or_identifiers(
        endpoint="analyses",
        identifiers=identifiers,
        filters=filters,
        fields="storage_url",
    ):
        if i["storage_url"]:
            if pattern:
                click.echo("\n".join(glob(join(i["storage_url"], pattern))))
            else:
                click.echo(i["storage_url"])


@click.command()
@options.NULLABLE_FILTERS
@options.NULLABLE_IDENTIFIERS
@options.VERBOSE
@click.option("--dtypes", help="Limit data types to be printed.", multiple=True)
def get_data(filters, identifiers, verbose, dtypes):
    """Get file paths for experiments raw data."""
    for i in _filters_or_identifiers(
        endpoint="experiments",
        identifiers=identifiers,
        filters=filters,
        fields="raw_data,system_id",
    ):
        system_id = i["system_id"]

        if not i["raw_data"] and not verbose:
            raise click.UsageError(f"No data for {system_id}, ignore with --verbose")

        for j in i["raw_data"] or ["None"]:
            if not dtypes or j["file_type"] in dtypes:
                click.echo(j["file_url"] if not verbose else f"{system_id}\t{j}")


@click.command()
@options.BED_TYPE
@click.option("--assembly", help="required if multiple options for assembly")
@click.argument("technique", required=True)
def get_bed(technique, bed_type, assembly):
    """Get a BED file for a given Sequencing Tehcnique."""
    instance = api.get_instance("techniques", technique)
    data_id = f"{assembly}_{bed_type}_bedfile"
    paths = {}

    for i, j in instance.reference_data.items():
        if i.endswith(f"{bed_type}_bedfile"):
            paths[i] = j["url"]

    if not paths:
        raise click.UsageError("No BED files registered yet...")
    elif len(paths) > 1 and not assembly:
        raise click.UsageError(f"Multiple BEDs for {technique}, pass --assembly")

    click.echo(paths[data_id] if assembly else list(paths.values())[0])


@click.command()
@click.argument("identifier", required=True)
@click.option(
    "--data-id",
    help="data identifier of the reference resource",
    default="genome_fasta",
    show_default=True,
)
@click.option(
    "--resources",
    help="Print the list of available reference files for this assembly.",
    is_flag=True,
)
@click.option(
    "--model",
    help="default model is assemblies",
    type=click.Choice(["assemblies", "techniques"]),
    default="assemblies",
)
def get_reference(identifier, data_id, resources, model):
    """Retrieve reference data from assemblies (default) or techniques."""
    instance = api.get_instance(model, identifier)

    if resources:
        click.echo(
            "\n".join(
                f"{click.style(i, fg='green')}\t{j.description}"
                for i, j in sorted(instance["reference_data"].items())
            ).expandtabs(30)
        )
    else:
        try:
            click.echo(instance["reference_data"][data_id]["url"])
        except KeyError:  # pragma: no cover
            raise click.UsageError(f"No {data_id} reference for {instance['name']}.")


@click.command()
@options.NULLABLE_IDENTIFIERS
@options.NULLABLE_FILTERS
@options.VERBOSE
@click.option("-r", "--result-key", help="result identifier", required=True)
@click.option(
    "--app-results",
    help="Provide an application primary key to get a list of available results.",
    callback=cb_app_results_keys,
    expose_value=False,
    is_eager=True,
    type=click.INT,
)
def get_results(filters, identifiers, result_key, verbose):
    """Get analyses results."""
    for i in _filters_or_identifiers(
        endpoint="analyses",
        identifiers=identifiers,
        filters=filters,
        fields="results,pk",
    ):
        if result_key in i["results"]:
            result = i["results"][result_key]
            click.echo(result if not verbose else f"{i['pk']}\t{result}")
        elif verbose:  # pragma: no cover
            click.echo(f"{i['pk']}\t- No result available")
        else:  # pragma: no cover
            raise click.UsageError(
                f"No '{result_key}' for {i['pk']}, ignore with --verbose"
            )


@click.command()
@options.NULLABLE_IDENTIFIERS
@options.NULLABLE_FILTERS
@options.VERBOSE
@click.option("--assembly", help="required if multiple options for assembly")
def get_bams(filters, assembly, verbose, identifiers):
    """Get storage directories, use `pattern` to match files inside dirs."""
    for i in _filters_or_identifiers(
        endpoint="experiments",
        identifiers=identifiers,
        filters=filters,
        fields="bam_files,system_id",
    ):
        bam_path = None
        system_id = i["system_id"]

        if assembly:
            bam_path = i["bam_files"][assembly]["url"]
        elif len(i["bam_files"]) == 1:
            bam_path = list(i["bam_files"].values())[0]["url"]

        if bam_path or verbose:
            click.echo(bam_path if not verbose else f"{system_id}\t{bam_path}")
        elif not i["bam_files"]:
            raise click.UsageError(f"No bams for {system_id}, ignore with --verbose")
        else:
            raise click.UsageError(f"Multiple bams for {system_id}, pass --assembly")


@click.command()
@options.NULLABLE_FILTERS
def rerun_signals(filters):
    """Rerun failed signals."""
    for i in api.get_instances(
        "signals", pk__gt=0, data__failure_traceback__isnull=False, **filters
    ):
        click.secho(f"Rerunning signal: {i.slug}", fg="yellow")

        try:
            instance = api.get_instance(i.target_endpoint, i.target_id)
            api._run_signals(
                endpoint=i.target_endpoint,
                instance=instance,
                signals=[import_from_string(i.import_string)],
                raise_error=True,
            )

            api.delete_instance("signals", i.pk)
        except HTTPError as error:
            # Delete the signal if the object doesn't exist anymore
            if "Object not found try a different ID" in error.response.text:
                api.delete_instance("signals", i.pk)
        except exceptions.AutomationError:
            pass


@click.command()
@options.NULLABLE_FILTERS
def run_web_signals(filters):
    """Run signals triggered from the frontend."""
    for i in api.get_instances(
        "signals",
        import_string__in=[
            "isabl_cli.signals.resume_analysis_signal",
            "isabl_cli.signals.force_analysis_signal",
        ],
        **filters,
    ):
        click.secho(f"Running web signal: {i.slug}", fg="yellow")
        instance = api.get_instance(i.target_endpoint, i.target_id)

        try:
            api._run_signals(
                endpoint=i.target_endpoint,
                instance=instance,
                signals=[import_from_string(i.import_string)],
                raise_error=True,
            )

            api.delete_instance("signals", i.pk)
        except exceptions.AutomationError:
            pass


@click.command()
@options.NULLABLE_FILTERS
@click.option(
    "--signals",
    "-s",
    help="Signal import string (e.g. isabl_cli.data.trigger_analyses_merge)",
    required=True,
    multiple=True,
)
@click.argument(
    "endpoint", required=True, type=click.Choice(["experiments", "analyses"])
)
def run_signals(endpoint, filters, signals):
    """Run any arbitrary signal on analyses or experiments using import strings."""
    count = api.get_instances_count(endpoint, **filters)
    click.secho(f"Running {', '.join(signals)} for {count} {endpoint}...", fg="blue")
    signals = [import_from_string(i) for i in signals]

    for i in api.get_instances(endpoint, **filters):
        api._run_signals(endpoint, i, signals, raise_error=True, create_record=False)


@click.command()
@options.FAILED_ANALYSES
@options.FORCE
@options.RESTART
def run_failed_analyses(failed_analyses_filters, force, restart):
    """Command to run failed analyses in batch."""
    utils.check_admin()
    # group analyses per application
    analyses = {}
    for i in failed_analyses_filters:
        app_class = i.application.application_class
        if app_class in analyses:
            analyses[app_class].append(i)
        else:
            analyses[app_class] = [i]

    for app_class, app_analyses in analyses.items():
        tuples = [(i.targets, i.references) for i in app_analyses]
        app = import_from_string(app_class)
        app().run(
            tuples=tuples,
            commit=False,
            restart=restart,
            force=force,
            run_args=i.data.get("run_args", {}),
        )


@click.command(hidden=True)
@options.ANALYSIS_PRIMARY_KEY
@click.option("--reason", help="Rejection reason. (Will be stored in Analysis.notes)")
def reject_analysis(key, reason):
    """Patch an analysis status as REJECTED, providing the rejection reason."""
    analysis = api.get_instance("analyses", key)
    api.patch_analysis_status(analysis, "REJECTED")
    api.patch_instance("analyses", key, notes=reason)
