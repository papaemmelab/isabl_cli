"""Engine validator."""

# pylint: disable=R0201,C0103

from leuktools import exceptions


class Validator():

    """Class with reusable tuples validators."""

    def validate_same_assembly(self, targets, references, analyses):
        """Validate targets, references tuple is a pair."""
        assemblies = {i.genome_assembly for i in targets + references}

        if len(assemblies) > 1:
            raise exceptions.UsageError(f"Multiple assemblies: {assemblies}")

        return True

    def validate_tuple_ispair(self, targets, references, analyses):
        """Validate targets, references tuple is a pair."""
        if len(targets) != 1 or len(references) != 1:
            raise exceptions.UsageError("Pipeline requires tumor normal pairs")

        return True

    def validate_onetarget_noreferences(self, targets, references, analyses):
        """Test only one sample is passed targets and none on references."""
        if len(targets) != 1 and references:
            raise exceptions.UsageError("Pipeline requires no references.")

        return True

    def validate_atleast_onetarget_onereference(
            self, targets, references, analyses):
        """Validate that at least one reference and target are passed."""
        if not references or not targets:
            raise exceptions.UsageError(
                "Pipeline requires at least one reference and one target."
                )

        return True

    def validate_has_indexed_bam(self, targets, references, analyses):
        """Make sure tuples have indexed bams."""
        try:
            paths = (i.get_bampath(indexed=True) for i in targets + references)
            assert all(paths)
        except (exceptions.PackageBaseException, AssertionError) as error:
            raise exceptions.UsageError(error)

        return True

    def validate_targets_notin_references(self, targets, references, analyses):
        """Make sure targets are not passed as references also."""
        refset = set(i.pk for i in references)
        msg = []

        for i in targets:
            if i.pk in refset:
                msg.append(f"{i} was also used as reference.")

        if msg:
            raise exceptions.UsageError("\n".join(msg))

        return True

    def validate_dna_tuples(self, targets, references, analyses):
        """Make sure all targets and references are DNA data."""
        msg = []

        for i in targets + references:
            if i.extraction.analyte != "DNA":
                msg .append(f"{i} analyte is not DNA")

        if msg:
            raise exceptions.UsageError("\n".join(msg))

        return True

    def validate_dna_pairs(self, targets, references, analyses):
        """Validate targets, references tuples for base dna pipelines."""
        self.validate_tuple_ispair(targets, references, analyses)
        self.validate_dna_tuples(targets, references, analyses)
        return True

    def validate_same_technique(self, targets, references, analyses):
        """Validate targets and references have same bedfile."""
        techniques = {i.technique.pk for i in references}

        if len(techniques) > 1:
            raise exceptions.UsageError(
                f"References have different techniques: {techniques}"
                )

        msg = []
        for target in targets:
            if target.technique.pk not in techniques:
                msg.append(
                    f"References technique differ from {target}: "
                    f"{target.technique} =! {references[0].technique}."
                    )

        if msg:
            raise exceptions.UsageError("\n".join(msg))

        return True
