"""Engine analysis creator."""

from collections import defaultdict
import sys

from cached_property import cached_property
import click

from cli import api
from cli import data
from cli import exceptions

from .validator import Validator


class Creator(Validator):

    """Analysis creation logic."""

    @cached_property
    def pipeline(self):
        """Get pipeline database object."""
        pipeline_class = f'{self.__module__}.{self.__class__.__name__}'

        if not all([self.NAME, self.VERSION, self.ASSEMBLY, self.SPECIES]):
            raise NotImplementedError(
                f"NAME must be set: {self.NAME}\n"
                f"VERSION must be set: {self.VERSION}\n"
                f"ASSEMBLY must be set: {self.ASSEMBLY}\n"
                f"SPECIES must be set: {self.SPECIES}\n")

        pipeline = api.create_instance(
            endpoint='pipelines',
            name=self.NAME,
            version=self.VERSION,
            assembly={'name': self.ASSEMBLY, 'species': self.SPECIES},
            pipeline_class=pipeline_class)

        if pipeline['pipeline_class'] != pipeline_class:
            api.patch_instance(
                endpoint='pipelines',
                identifier=pipeline['pk'],
                pipeline_class=pipeline_class)

        return pipeline

    @cached_property
    def assembly(self):
        """Get assembly database object."""
        return self.pipeline['assembly']

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
                    self._validate_tuple(*i)
                    analysis = api.create_instance(
                        endpoint='analyses',
                        pipeline=self.pipeline,
                        targets=i[0],
                        references=i[1],
                        analyses=i[2])

                    created_analyses.append(
                        data.update_storage_url(
                            endpoint='analyses',
                            identifier=analysis['pk'],
                            use_hash=True))
                except exceptions.ValidationError as error:
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

    def _validate_tuple(self, targets, references, analyses):
        self.validate_species(targets, references, analyses)
        self.validate_tuple(targets, references, analyses)
