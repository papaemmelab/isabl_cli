"""Engine validator."""
# pylint: disable=R0201,C0103

import click


class Validator():

    """Class with reusable tuples validators."""

    def validate_tuple_ispair(self, targets, references, analyses):
        """Validate targets, references tuple is a pair."""
        if len(targets) != 1 or len(references) != 1:
            raise click.UsageError("Pipeline requires tumor normal pairs")
        return True

    def validate_onetarget_noreferences(self, targets, references, analyses):
        """Test only one sample is passed targets and none on references."""
        if len(targets) != 1 and references:
            raise click.UsageError("Pipeline requires no references.")
        return True

    def validate_atleast_onetarget_onereference(
            self, targets, references, analyses):
        """Validate that at least one reference and target are passed."""
        if not references or not targets:
            raise click.UsageError(
                "Pipeline requires at least one reference and one target.")
        return True

    def validate_targets_notin_references(self, targets, references, analyses):
        """Make sure targets are not passed as references also."""
        refset = set(i['pk'] for i in references)
        msg = []

        for i in targets:
            if i['pk'] in refset:
                msg.append(f"{i['system_id']} was also used as reference.")

        if msg:
            raise click.UsageError('\n'.join(msg))

        return True

    def validate_dna_tuples(self, targets, references, analyses):
        """Make sure all targets and references are DNA data."""
        msg = []

        for i in targets + references:
            if i['technique']['analyte'] != "DNA":
                msg .append(f"{i['system_id']} analyte is not DNA")

        if msg:
            raise click.UsageError('\n'.join(msg))

        return True

    def validate_dna_pairs(self, targets, references, analyses):
        """Validate targets, references tuples for base dna pipelines."""
        self.validate_tuple_ispair(targets, references, analyses)
        self.validate_dna_tuples(targets, references, analyses)
        return True

    def validate_same_technique(self, targets, references, analyses):
        """Validate targets and references have same bedfile."""
        techniques = {i['technique']['slug'] for i in references}

        if len(techniques) > 1:
            raise click.UsageError(
                f"References have different techniques: {techniques}")

        msg = []
        for i in targets:
            if i['technique']['slug'] not in techniques:
                msg.append(
                    f"References technique differ from {i['system_id']}: "
                    f"{i['technique']['slug']} =! "
                    f"{references[0]['technique']['slug']}.")

        if msg:
            raise click.UsageError('\n'.join(msg))

        return True
