"""commands logic."""

from glob import glob
from os.path import join

import click

from isabl_cli import api
from isabl_cli import options
from isabl_cli import utils
from isabl_cli.settings import import_from_string


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
        results = application.get_analysis_results(i)
        api.patch_instance("analyses", i["pk"], results=results)


@click.command(hidden=True)
@options.ANALYSIS_PRIMARY_KEY
@options.ANALYSIS_STATUS
def patch_status(key, status):
    """Patch status of a given analysis."""
    analysis = api.get_instance("analyses", key)
    api.patch_analysis_status(analysis, status)


@click.command()
@options.ENDPOINT
@options.FIELDS
@options.NULLABLE_FILTERS
@options.NO_HEADERS
def get_attributes(endpoint, field, filters, no_headers):
    """Get database attributes from API."""
    result = [] if no_headers else ["\t".join(".".join(i) for i in field)]
    filters["fields"] = ",".join(i[0] for i in field)

    for i in api.get_instances(endpoint, verbose=True, **filters):
        values = [utils.traverse_dict(i, j, serialize=True) for j in field]
        result.append("\t".join(values))

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
def get_sequencing_data(filters, verbose):
    """Get storage directories, use `pattern` to match files inside dirs."""
    filters.update(fields="sequencing_data,system_id", limit=100_000)
    for i in api.get_instances("experiments", verbose=True, **filters):
        system_id = i["system_id"]

        if not i["sequencing_data"] and not verbose:
            raise click.UsageError(f"No data for {system_id}, ignore with --verbose")

        for j in i["sequencing_data"] or ["None"]:
            click.echo(j if not verbose else f"{system_id} {j}")


@click.command()
@options.FILTERS
@options.VERBOSE
@click.option("--assembly", help="when bams available for multiple assemblies")
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
