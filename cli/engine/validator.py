# pylint: disable=R0201

from collections import defaultdict
from os.path import isdir
from os.path import isfile

from cli.exceptions import ValidationError


class Validator:

    ASSEMBLY = None
    SPECIES = None

    def validate_is_file(self, path):  # pragma: no cover
        """Validate path is file."""
        if not isfile(path):
            raise ValidationError(f'{path} is not a file.')

    def validate_is_dir(self, path):  # pragma: no cover
        """Validate path is directory."""
        if not isdir(path):
            raise ValidationError(f'{path} is not a file.')

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

    def validate_has_raw_sequencing_data(self, targets, references):
        """Validate targets and references have sequencing data."""
        msg = []

        for i in targets + references:
            if not i['sequencing_data']:
                msg.append(f'{i["system_id"]} has no sequencing data...')

        if msg:
            raise ValidationError('\n'.join(msg))

    def validate_single_data_type(self, targets, references):
        """Validate targets and references have sequencing data."""
        self.validate_has_raw_sequencing_data(targets, references)
        types = defaultdict(list)

        for i in targets + references:
            for j in i['sequencing_data']:
                file_type = j['file_type']

                if file_type.startswith('FASTQ_'):
                    file_type = 'FASTQ'

                types[file_type].append(i['system_id'])

        if len(types) > 1:
            raise ValidationError(
                f'Multiple data types not supported: {dict(types)}')

        return list(types.keys())[0]

    def validate_fastq_only(self, targets, references):
        """Validate sequencing data is only fastq."""
        dtype = self.validate_single_data_type(targets, references)

        if dtype != 'FASTQ':
            raise ValidationError(f'Only FASTQ supported, found: {dtype}')

    def validate_tuple_is_pair(self, targets, references):
        """Validate targets, references tuple is a pair."""
        if len(targets) != 1 or len(references) != 1:
            raise ValidationError('Target, reference pairs required.')

    def validate_one_target_no_references(self, targets, references):
        """Test only one sample is passed targets and none on references."""
        if len(targets) != 1 or references:
            raise ValidationError('References not allowed.')

    def validate_at_least_one_target_one_reference(self, targets, references):
        """Validate that at least one reference and target are passed."""
        if not references or not targets:
            raise ValidationError('References and targets required.')

    def validate_targets_not_in_references(self, targets, references):
        """Make sure targets are not passed as references also."""
        refset = set(i['pk'] for i in references)
        msg = []

        for i in targets:
            if i['pk'] in refset:
                msg.append(f"{i['system_id']} was also used as reference.")

        if msg:
            raise ValidationError('\n'.join(msg))

    def validate_methods(self, targets, references, methods):
        """Make sure all targets and references are come from PDX samples."""
        msg = []

        for i in targets + references:
            if i['technique']['method'] not in methods:
                msg.append(
                    f"Only '{methods}' sequencing method allowed, "
                    f"found {i['technique']['method']} for {i['system_id']}.")

        if msg:
            raise ValidationError('\n'.join(msg))

    def validate_pdx_only(self, targets, references):
        """Make sure all targets and references are come from PDX samples."""
        msg = []

        for i in targets + references:
            if not i['specimen']['is_pdx']:
                msg.append(f"{i['system_id']} specimen is not PDX derived")

        if msg:
            raise ValidationError('\n'.join(msg))

    def validate_dna_only(self, targets, references):
        """Make sure all targets and references are DNA data."""
        msg = []

        for i in targets + references:
            if i['technique']['analyte'] != 'DNA':
                msg.append(f"{i['system_id']} analyte is not DNA")

        if msg:
            raise ValidationError('\n'.join(msg))

    def validate_rna_only(self, targets, references):
        """Make sure all targets and references are DNA data."""
        msg = []

        for i in targets + references:
            if i['technique']['analyte'] != 'RNA':
                msg.append(f"{i['system_id']} analyte is not RNA")

        if msg:
            raise ValidationError('\n'.join(msg))

    def validate_dna_pairs(self, targets, references):
        """Validate targets, references tuples for base dna pipelines."""
        self.validate_tuple_is_pair(targets, references)
        self.validate_dna_only(targets, references)

    def validate_same_technique(self, targets, references):
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

    def validate_species(self, targets, references):
        """Validate targets and references have sequencing data."""
        msg = []

        for i in targets + references:
            if i['specimen']['individual']['species'] != self.SPECIES:
                msg.append(f'{i["system_id"]} species not supported')

        if msg:
            raise ValidationError('\n'.join(msg))
