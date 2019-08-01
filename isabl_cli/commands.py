"""commands logic."""

from collections import OrderedDict
from glob import glob
from os.path import join
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
    if not value or ctx.resilient_parsing:
        return

    click.echo(
        "\n".join(
            f"{click.style(i, fg='green')}\t{j.description}"
            for i, j in sorted(api.get_instance("applications", value).results.items())
        ).expandtabs(30)
    )
    ctx.exit()


def _filters_or_identifiers(endpoint, identifiers, filters, fields=None):
    if filters and identifiers:
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
@options.NULLABLE_FILTERS
def process_finished(filters):
    """Process and update finished analyses."""
    utils.check_admin()
    filters.update(status="FINISHED")

    for i in api.get_instances("analyses", verbose=True, **filters):
        if i["status"] == "FINISHED":
            api.patch_analysis_status(i, "SUCCEEDED")


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
            else:
                skipped.append(i)

    if skipped:
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
    "https://github.com/antonmedv/fx/blob/master/docs.md#interactive-mode"
)
@options.ENDPOINT
@options.FIELDS
@options.NULLABLE_FILTERS
@options.NO_HEADERS
@options.NULLABLE_IDENTIFIERS
@click.option("--json", "json_", help="Print as JSON", is_flag=True)
@click.option("--fx", "use_fx", help="Visualize json with fx", is_flag=True)
def get_metadata(identifiers, endpoint, field, filters, no_headers, json_, use_fx):
    """Retrieve metadata for multiple instances."""
    if not field and not (json_ or use_fx):
        raise click.UsageError("Pass --field or use --json/--fx")

    if use_fx and not shutil.which("fx"):
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

    if json_:
        click.echo(json.dumps(instances, sort_keys=True, indent=4))
    elif use_fx:
        fp = tempfile.NamedTemporaryFile("w+", delete=False)
        json.dump(instances, fp)
        fp.close()
        subprocess.check_call(["fx", fp.name])
        os.unlink(fp.name)
    else:
        result = [] if no_headers else ["\t".join(".".join(i) for i in field)]
        result += ["\t".join(map(str, i.values())) for i in instances]
        click.echo("\n".join(result).expandtabs(30))


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
def get_data(filters, identifiers, verbose):
    """Get file paths for experiments sequencing data."""
    for i in _filters_or_identifiers(
        endpoint="experiments",
        identifiers=identifiers,
        filters=filters,
        fields="sequencing_data,system_id",
    ):
        system_id = i["system_id"]

        if not i["sequencing_data"] and not verbose:
            raise click.UsageError(f"No data for {system_id}, ignore with --verbose")

        for j in i["sequencing_data"] or ["None"]:
            click.echo(j["file_url"] if not verbose else f"{system_id}\t{j}")


@click.command()
@options.BED_TYPE
@click.option("--assembly", help="required if multiple options for assembly")
@click.argument("technique", required=True)
def get_bed(technique, bed_type, assembly):
    """Get a BED file for a given Sequencing Tehcnique."""
    instance = api.get_instance("techniques", technique)

    if not instance["bed_files"]:
        raise click.UsageError("No BED files registered yet...")
    elif len(instance["bed_files"]) > 1 and not assembly:
        raise click.UsageError(f"Multiple BEDs for {technique}, pass --assembly")
    else:
        assembly = list(instance["bed_files"].keys())[0]

    click.echo(instance["bed_files"][assembly][bed_type])


@click.command()
@click.argument("assembly", required=True)
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
def get_reference(assembly, data_id, resources):
    """Get reference resource for an Assembly."""
    try:
        assembly = api.get_instance("assemblies", assembly)
    except KeyError:
        click.UsageError(f"No {data_id} reference for {assembly['name']}.")

    if resources:
        click.echo(
            "\n".join(
                f"{click.style(i, fg='green')}\t{j.description}"
                for i, j in sorted(assembly["reference_data"].items())
            ).expandtabs(30)
        )
    else:
        click.echo(assembly["reference_data"][data_id]["url"])


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
        elif verbose:
            click.echo(f"{i['pk']}\t- No result available")
        else:
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
