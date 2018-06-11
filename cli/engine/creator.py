"""Engine analysis creator."""

from collections import defaultdict
import os
import sys

from cached_property import cached_property
import click

from cli import system_settings
from cli import api

from .validator import Validator


class Creator(Validator):

    """Analysis creation logic."""

    NAME = None
    VERSION = None
    ASSEMBLY = None

    @cached_property
    def pipeline(self):
        """Get pipeline database object."""
        created_by = system_settings.api_username

        if not self.NAME or not self.VERSION or not self.ASSEMBLY:
            raise NotImplementedError("NAME, VERSION and ASSEMBLY must be set.")

        return api.create_instance(
            endpoint='pipelines',
            name=self.NAME,
            version=self.VERSION,
            created_by=created_by,
            assembly={'name': self.ASSEMBLY, 'created_by': created_by})

    @staticmethod
    def validate_tuple(targets, references, analyses):
        """Raise error tuple isn't valid, see pipeline.py for documentation."""
        raise NotImplementedError('`validate_tuples` not implemented.')

    def get_or_create_analyses(self, tuples):
        """
        Get or create analysis for a pipeline and a list of tuples.

        Arguments:
            tuples (list): list of (targets, references, analyses) tuples.

        Returns:
            tuple: list of analyses, invalid tuples [(tuple, error), ...].
        """
        label = f'Getting analyses for {len(tuples)} tuples...\t\t'
        existing_analyses, tuples = self.get_existing_analyses(tuples)
        created_analyses = []
        invalid_tuples = []

        with click.progressbar(tuples, file=sys.stderr, label=label) as bar:
            for i in bar:
                try:
                    assert self.validate_tuple(*i)
                    analysis = api.create_instance(
                        endpoint='analyses',
                        pipeline=self.pipeline,
                        targets=i[0],
                        references=i[1],
                        analyses=i[2],
                        created_by=system_settings.api_username)

                    storage_url = system_settings.GET_STORAGE_DIRECTORY(
                        endpoint='analyses',
                        primary_key=analysis['pk'])

                    analysis = api.patch_instance(
                        endpoint='analyses',
                        identifier=analysis['pk'],
                        storage_url=storage_url,
                        storage_usage=0)

                    os.makedirs(storage_url, exist_ok=True)
                    created_analyses.append(analysis)
                except click.UsageError as error:
                    invalid_tuples.append((i, error))

        return existing_analyses + created_analyses, invalid_tuples

    def get_existing_analyses(self, tuples):
        """
        Get list of existing analyses for a list of `tuples`.

        Check targets and references for an exact match. Check whether all the
        requested analyses are in `analysis.analyses`.

        Arguments:
            tuples (list): list of (targets, references, analyses) tuples.

        Returns:
            tuple: list of existing analyses, list of tuples without analysis
        """
        click.echo("Checking for existing analyses...", file=sys.stderr)
        projects = {j['pk'] for i, _, __ in tuples for j in i[0]['projects']}
        filters = dict(pipeline=self.pipeline['pk'], projects__pk__in=projects)
        cache = defaultdict(list)
        existing, missing = [], []

        def get_cache_key(targets, references):
            return (
                tuple(sorted(set(i['pk'] for i in targets))),
                tuple(sorted(set(i['pk'] for i in references))))

        for i in api.get_instances('analyses', **filters):
            cache[get_cache_key(i['targets'], i['references'])].append(i)

        for targets, references, analyses in tuples:
            expected = set(analyses)
            found_existing = False

            for i in cache[get_cache_key(targets, references)]:
                if expected.issubset(set(i['analyses'])):
                    found_existing = True
                    existing.append(i)
                    break

            if not found_existing:
                missing.append((targets, references, analyses))

        return existing, missing
