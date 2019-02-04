"""commands logic."""

from glob import glob
from os.path import join
from collections import OrderedDict
import os
import json
import shutil
import subprocess
import tempfile

import click

from isabl_cli import api
from isabl_cli import options
from isabl_cli import utils
from isabl_cli.settings import import_from_string
from isabl_cli.settings import user_settings


@click.command()
def login():  # pragma: no cover
    """Login with isabl credentials."""
    user_settings.api_token = None
    api.get_token_headers()


@click.command()
@click.option("--project", help="primary key of project to merge by", type=int)
@click.option("--application", help="analyses application primary key", type=int)
def merge_project_analyses(project, application):  # pragma: no cover
    """Merge analyses by project primary key."""
    project = api.get_instance("projects", project)
    application = api.get_instance("applications", application)
    application = import_from_string(application["application_class"])()
    application.run_project_merge(project)


@click.command()
@options.FILTERS
def processed_finished(filters):
    """Process and update finished analyses."""
    utils.check_admin()
    filters.update(status="FINISHED")

    for i in api.get_instances("analyses", **filters):
        api.patch_analysis_status(i, "SUCCEEDED")


@click.command()
@options.FILTERS
def patch_results(filters):
    """Update the results field of many analyses."""
    utils.check_admin()

    for i in api.get_instances("analyses", **filters):
        application = import_from_string(i["application"]["application_class"])()
        results = application._get_analysis_results(i)
        api.patch_instance("analyses", i["pk"], results=results)


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
@click.option("--fx", help="Visualize json with fx", is_flag=True)
def get_metadata(
    identifiers, endpoint, field, filters, no_headers, json_, fx
):  # pylint: disable=invalid-name
    """Retrieve metadata for multiple instances."""
    if not field and not (json_ or fx):
        raise click.UsageError("Pass --field or use --json/--fx")

    if fx and not shutil.which("fx"):
        raise click.UsageError("fx is not installed")

    if filters and identifiers:
        raise click.UsageError("can't combine filters and identifiers")

    fields = [i[0] for i in field]  # first level fields
    identifiers = identifiers or None

    if identifiers:
        instances = [api.get_instance(endpoint, i, fields=fields) for i in identifiers]
    else:
        filters.update({"fields": ",".join(fields)} if field else {})
        instances = api.get_instances(endpoint, identifiers, verbose=True, **filters)

    results = instances

    if field:  # if fields were passed, update the results list
        results = [
            OrderedDict([(".".join(j), utils.traverse_dict(i, j)) for j in field])
            for i in instances
        ]

    if json_:
        click.echo(json.dumps(results, sort_keys=True, indent=4))
    elif fx:
        fp = tempfile.NamedTemporaryFile("w+", delete=False)
        json.dump(results, fp)
        fp.close()
        subprocess.check_call(["fx", fp.name])
        os.unlink(fp.name)
    else:
        result = [] if no_headers else ["\t".join(".".join(i) for i in field)]
        result += ["\t".join(map(str, i.values())) for i in results]
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
@options.FILTERS
def get_paths(endpoint, pattern, filters):
    """Get storage directories, use `pattern` to match files inside dirs."""
    filters.update(fields="storage_url", limit=100_000)
    for i in api.get_instances(endpoint, verbose=True, **filters):
        if i["storage_url"]:
            if pattern:
                click.echo("\n".join(glob(join(i["storage_url"], pattern))))
            else:
                click.echo(i["storage_url"])


@click.command()
@options.FILTERS
@options.VERBOSE
def get_data(filters, verbose):
    """Get file paths for experiments sequencing data."""
    filters.update(fields="sequencing_data,system_id", limit=100_000)
    for i in api.get_instances("experiments", verbose=True, **filters):
        system_id = i["system_id"]

        if not i["sequencing_data"] and not verbose:
            raise click.UsageError(f"No data for {system_id}, ignore with --verbose")

        for j in i["sequencing_data"] or ["None"]:
            click.echo(j["file_url"] if not verbose else f"{system_id} {j}")


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
@click.option("--data-id", help="data identifier", default="genome_fasta")
def get_reference(assembly, data_id):
    """Get reference resource for an Assembly."""
    try:
        assembly = api.get_instance("assemblies", assembly)
        click.echo(assembly["reference_data"][data_id]["url"])
    except KeyError:
        click.UsageError(f"No {data_id} reference for {assembly['name']}.")


@click.command()
@options.FILTERS
@options.VERBOSE
@click.option("--assembly", help="required if multiple options for assembly")
def get_bams(filters, assembly, verbose):
    """Get storage directories, use `pattern` to match files inside dirs."""
    filters.update(fields="bam_files,system_id", limit=100_000)
    for i in api.get_instances("experiments", verbose=True, **filters):
        bam_path = None
        system_id = i["system_id"]

        if assembly:
            bam_path = i["bam_files"][assembly]["url"]
        elif len(i["bam_files"]) == 1:
            bam_path = list(i["bam_files"].values())[0]["url"]

        if bam_path:
            click.echo(bam_path if not verbose else f"{system_id} {bam_path}")
        elif not i["bam_files"]:
            raise click.UsageError(f"No bams for {system_id}, ignore with --verbose")
        else:
            raise click.UsageError(f"Multiple bams for {system_id}, pass --assembly")
