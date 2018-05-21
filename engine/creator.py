"""Engine analysis creator."""

from collections import defaultdict
from itertools import chain
import sys

from cached_property import cached_property
import click

from leuktools import Analysis
from leuktools import Pipeline
from leuktools import Workflow
from leuktools import api
from leuktools import constants
from leuktools import exceptions
from leuktools import get_analyses
from leuktools import get_workflows

from .validator import Validator


class Creator(Validator):

    """Analysis creation logic."""

    @cached_property
    def pk(self):
        """Get pipeline primary key."""
        raise exceptions.NotImplementedError()

    @cached_property
    def pipeline(self):
        """Get pipeline database object."""
        return Pipeline(self.pk)

    @staticmethod
    def validate_tuple(targets, references, analyses):
        """Raise error tuple isn't valid, see pipeline.py for documentation."""
        raise exceptions.NotImplementedError()

    def _get_or_create_analyses(self, tuples):
        """
        Get or create analysis for a pipeline and a list of tuples.

        Arguments:
            tuples (list): list of (targets, references, analyses) tuples.

        Returns:
            tuple: list of analyses, invalid tuples [(tuple, error), ...].
        """
        click.echo("Checking for existing analyses...", file=sys.stderr)
        tuples = self._coerce_to_objects(tuples)
        analyses, tuples = self._get_existing_analyses(tuples)
        label = "Creating analyses...\t\t"
        invalid = []

        with click.progressbar(tuples, file=sys.stderr, label=label) as bar:
            for i in bar:
                try:
                    analysis = self._create_analysis(*i)
                    analysis.as_target_objects = i[0]
                    analysis.as_reference_objects = i[1]
                    analysis.analyses_objects = i[2]
                    analyses.append(analysis)
                except exceptions.UsageError as error:
                    invalid.append(i, error)
                    continue

        return analyses, invalid

    def _get_existing_analyses(self, tuples):
        """
        Get list of existing analyses for a list of `tuples`.

        Check targets and references for an exact match. Check whether all the
        requested analyses are in `analysis.analyses`.

        Arguments:
            tuples (list): list of (targets, references, analyses) tuples.

        Returns:
            tuple: list of existing analyses, list of tuples without analysis
        """
        existing, missing = [], []
        cache = self._build_cache(tuples)

        for targets, references, analyses in tuples:
            expected = set(analyses)
            analysis = None

            for i in cache.get(self._get_cache_key(targets, references), []):
                if expected.issubset(set(i.analyses)):
                    analysis = i
                    analysis.as_target_objects = i[0]
                    analysis.as_reference_objects = i[1]
                    analysis.analyses_objects = i[2]

            if analysis:
                existing.append(analysis)
            else:
                missing.append((targets, references, analyses))

        return existing, missing

    def _create_analysis(self, as_target, as_reference, analyses):
        """
        Create an analysis given a tuple of targets, references and analyses.

        Arguments:
            as_target (list): Workflow objects to be used as target.
            as_reference (list): Workflow objects to be used as reference.
            analyses (list): Analysis objects to be linked.

        Raises:
            exceptions.ValidationError: if tuple is not valid.

        Returns:
            Analysis: fresh analysis object.
        """
        assert self.validate_tuple(as_target, as_reference, analyses)

        data = dict(
            pipeline=self.pk,
            as_target=self._coerce_to_keys(as_target),
            as_reference=self._coerce_to_keys(as_reference),
            analyses=self._coerce_to_keys(analyses),
            )

        url = constants.ANALYSES_API_URL
        data = api.request("post", url=url, json=data).json()
        return Analysis(from_dict=data)

    def _build_cache(self, tuples):
        """
        Build a dict of candidate existing analyses.

        Arguments:
            tuples (list): list of (targets, references, analyses) tuples.

        Returns:
            dict: were keys are a tuple of targets and references primary keys
                and values are a list of analyses that match these workflows.
        """
        projects = set(chain(i[0].projects for i in tuples))
        filters = [("pipeline", self.pk), ("projects__pk__in", projects)]
        cache = defaultdict(list)

        for i in get_analyses(filters=filters):
            key = self._get_cache_key(i.as_target, i.as_reference)
            cache[key].append(i)

        return cache

    def _get_cache_key(self, targets, references):
        """Get a key to match chached analyses."""
        targets = tuple(sorted(self._coerce_to_keys(targets)))
        references = tuple(sorted(self._coerce_to_keys(references)))
        return targets, references

    @staticmethod
    def _coerce_to_objects(tuples):
        """Make sure that we have objects instead of just identifiers."""
        result = []
        workflows_objects = dict()
        analyses_objects = dict()
        workflows_ids = []
        analyses_ids = []

        for targets, references, analyses in tuples:
            for i in chain(targets, references):
                if not isinstance(i, Workflow):
                    workflows_ids.append(i)
                else:
                    workflows_objects[i] = i

            for i in analyses:
                if not isinstance(i, Analysis):
                    analyses_ids.append(i)
                else:
                    analyses_objects[i] = i

        for i in get_workflows(workflows_ids):
            workflows_objects[i.pk] = i
            workflows_objects[i.slug] = i
            workflows_objects[i.leukid] = i
            workflows_objects[str(i.pk)] = i

        for i in get_analyses(analyses_ids):
            analyses_objects[i.pk] = i
            analyses_objects[i.slug] = i
            analyses_objects[str(i.pk)] = i

        for targets, references, analyses in tuples:
            result.append((
                list({workflows_objects[i] for i in targets}),
                list({workflows_objects[i] for i in references}),
                list({analyses_objects[i] for i in analyses}),
                ))

        return result

    @staticmethod
    def _coerce_to_keys(objects):
        """Get list of primary keys from `objects`."""
        return list(set(api.get_pks(objects)))
