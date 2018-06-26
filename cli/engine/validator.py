"""Engine validator."""
# pylint: disable=R0201,C0103

from os.path import isfile

from cli.exceptions import ValidationError


class Validator():

    """Class with reusable tuples validators."""

    NAME = None
    VERSION = None
    ASSEMBLY = None
    SPECIES = None

    def validate_reference_genome(self, reference):
        """Validate genome exists and has required indexes."""
        required = ".fai", ".amb", ".ann", ".bwt", ".pac", ".sa"

        if not all(isfile(reference + i) for i in required):
            raise ValidationError(
                f'Missing indexes please run:\n\n\t'
                f'bwa index {reference}\n\tsamtools faidx {reference}')

        if not isfile(reference + '.dict'):
            raise ValidationError(
                f'Please generate {reference + ".dict"}, e.g.\n\n\t'
                f'samtools dict -a {self.ASSEMBLY} -s {self.SPECIES} '
                f'{reference} > {reference + ".dict"}')

    def validate_has_raw_sequencing_data(self, targets, references, analyses):
        """Validate targets and references have sequencing data."""
        msg = []

        for i in targets + references:
            if not i['sequencing_data']:
                msg.append(f'{i["system_id"]} has no sequencing data...')

        if msg:
            raise ValidationError('\n'.join(msg))

    def validate_single_data_type(
            self, targets, references, analyses, dtype=None):
        """Validate targets and references have sequencing data."""
        self.validate_has_raw_sequencing_data(targets, references, analyses)
        types = set()

        for i in targets + references:
            for j in i['sequencing_data']:
                types.add(j['file_type'])

        if len(types) > 1:
            raise ValidationError(f'Multiple data types not supported: {types}')

        if dtype and dtype not in types:
            raise ValidationError(f'Only {dtype} supported, found: {types}')

    def validate_tuple_is_pair(self, targets, references, analyses):
        """Validate targets, references tuple is a pair."""
        if len(targets) != 1 or len(references) != 1:
            raise ValidationError('Target, reference pairs required.')

    def validate_one_target_no_references(self, targets, references, analyses):
        """Test only one sample is passed targets and none on references."""
        if len(targets) != 1 or references:
            raise ValidationError('References not allowed.')

    def validate_atleast_onetarget_onereference(
            self, targets, references, analyses):
        """Validate that at least one reference and target are passed."""
        if not references or not targets:
            raise ValidationError('References and targets required.')

    def validate_targets_not_in_references(self, targets, references, analyses):
        """Make sure targets are not passed as references also."""
        refset = set(i['pk'] for i in references)
        msg = []

        for i in targets:
            if i['pk'] in refset:
                msg.append(f"{i['system_id']} was also used as reference.")

        if msg:
            raise ValidationError('\n'.join(msg))

    def validate_dna_only(self, targets, references, analyses):
        """Make sure all targets and references are DNA data."""
        msg = []

        for i in targets + references:
            if i['technique']['analyte'] != 'DNA':
                msg.append(f"{i['system_id']} analyte is not DNA")

        if msg:
            raise ValidationError('\n'.join(msg))

    def validate_rna_only(self, targets, references, analyses):
        """Make sure all targets and references are DNA data."""
        msg = []

        for i in targets + references:
            if i['technique']['analyte'] != 'RNA':
                msg.append(f"{i['system_id']} analyte is not RNA")

        if msg:
            raise ValidationError('\n'.join(msg))

    def validate_dna_pairs(self, targets, references, analyses):
        """Validate targets, references tuples for base dna pipelines."""
        self.validate_tuple_is_pair(targets, references, analyses)
        self.validate_dna_only(targets, references, analyses)

    def validate_same_technique(self, targets, references, analyses):
        """Validate targets and references have same bedfile."""
        techniques = {i['technique']['slug'] for i in references}

        if len(techniques) > 1:
            raise ValidationError(
                f'Multiple references techniques: {techniques}')

        msg = []
        for i in targets:
            if i['technique']['slug'] not in techniques:
                msg.append(
                    f"References technique differ from {i['system_id']}: "
                    f"{i['technique']['slug']} =! "
                    f"{references[0]['technique']['slug']}.")

        if msg:
            raise ValidationError('\n'.join(msg))

    def validate_species(self, targets, references, analyses):
        """Validate targets and references have sequencing data."""
        msg = []

        for i in targets + references:
            if i['specimen']['individual']['species'] != self.SPECIES:
                msg.append(f'{i["system_id"]} species not supported')

        if msg:
            raise ValidationError('\n'.join(msg))
