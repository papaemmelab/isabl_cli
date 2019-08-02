
import click

from isabl_cli.settings import import_from_string


def resume_analysis_signal(analysis):
    """Signal to resume analysis execution."""
    click.secho(f"Resuming analysis {analysis.slug}", fg="yellow")
    run_web_signals(analysis, force=False)


def force_analysis_signal(analysis):
    """Signal to wipe analysis and restart."""
    click.secho(f"Forcing analysis {analysis.slug}", fg="yellow")
    run_web_signals(analysis, force=True)


def run_web_signals(analysis, force=False):
    """Signal to trigger analyses executions."""
    tuples = [(analysis.targets, analysis.references)]
    app = import_from_string(analysis.application.application_class)
    app().run(tuples=tuples, commit=True, force=force)
