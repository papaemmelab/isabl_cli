# pylint: disable=R0201

from collections import defaultdict
from os.path import isdir
from os.path import isfile


class Validator:

    ASSEMBLY = None
    SPECIES = None

    def validate_is_file(self, path):  # pragma: no cover
        """Validate path is file."""
        assert isfile(path), f'{path} is not a file.'

    def validate_is_dir(self, path):  # pragma: no cover
        """Validate path is directory."""
        assert isdir(path), f'{path} is not a file.'

    def validate_reference_genome(self, reference):
        """Validate genome exists and has required indexes."""
        required = ".fai", ".amb", ".ann", ".bwt", ".pac", ".sa"

        assert all(isfile(reference + i) for i in required), (
            f'Missing indexes please run:\n\n\t'
            f'bwa index {reference}\n\tsamtools faidx {reference}')

        assert isfile(reference + '.dict'), (
            f'Please generate {reference + ".dict"}, e.g.\n\n\t'
            f'samtools dict -a {self.ASSEMBLY} -s {self.SPECIES} '
            f'{reference} > {reference + ".dict"}')

    def validate_has_raw_sequencing_data(self, workflows):
        """Validate workflows have sequencing data."""
        msg = []

        for i in workflows:
            if not i['sequencing_data']:
                msg.append(f'{i["system_id"]} has no sequencing data...')

        assert not msg, '\n'.join(msg)

    def validate_single_data_type(self, workflows):
        """Validate workflows have only one type of sequencing data."""
        self.validate_has_raw_sequencing_data(workflows)
        types = defaultdict(list)

        for i in workflows:
            for j in i['sequencing_data']:
                file_type = j['file_type']

                if file_type.startswith('FASTQ_'):
                    file_type = 'FASTQ'

                types[file_type].append(i['system_id'])

        assert len(types) == 1, f'Multiple types not supported: {dict(types)}'
        return list(types.keys())[0]

    def validate_fastq_only(self, workflows):
        """Validate sequencing data is only fastq."""
        dtype = self.validate_single_data_type(workflows)
        assert dtype == 'FASTQ', f'Only FASTQ supported, found: {dtype}'

    def validate_is_pair(self, targets, references):
        """Validate targets, references tuple is a pair."""
        assert len(targets) == 1 and len(references) == 1, 'Pairs only.'

    def validate_one_target(self, targets):
        """Validate only one target."""
        assert len(targets) == 1, f'Only 1 target allowed, got: {len(targets)}'

    def validate_one_target_no_references(self, targets, references):
        """Validate only one sample is passed targets and none on references."""
        self.validate_one_target(targets)
        assert not references, f'No reference workflows, got: {len(references)}'

    def validate_at_least_one_target_one_reference(self, targets, references):  # pylint: disable=C0103
        """Validate that at least one reference and target are passed."""
        assert targets and references, 'References and targets required.'

    def validate_targets_not_in_references(self, targets, references):
        """Make sure targets are not passed as references also."""
        refset = set(i['pk'] for i in references)
        template = "%s was also used as reference."
        msg = [template % i['system_id'] for i in targets if i['pk'] in refset]
        assert not msg, '\n'.join(msg)

    def validate_methods(self, workflows, methods):
        """Make sure all workflows methods are those expected."""
        msg = []

        for i in workflows:
            if i['technique']['method'] not in methods:
                msg.append(
                    f"Only '{methods}' sequencing method allowed, "
                    f"found {i['technique']['method']} for {i['system_id']}.")

        assert not msg, '\n'.join(msg)

    def validate_pdx_only(self, workflows):
        """Make sure workflows come from PDX samples."""
        msg = []

        for i in workflows:
            if not i['specimen']['is_pdx']:
                msg.append(f"{i['system_id']} specimen is not PDX derived")

        assert not msg, '\n'.join(msg)

    def validate_dna_only(self, workflows):
        """Make sure workflows are DNA data."""
        msg = []

        for i in workflows:
            if i['technique']['analyte'] != 'DNA':
                msg.append(f"{i['system_id']} analyte is not DNA")

        assert not msg, '\n'.join(msg)

    def validate_rna_only(self, workflows):
        """Make sure workflows are RNA data."""
        msg = []

        for i in workflows:
            if i['technique']['analyte'] != 'RNA':
                msg.append(f"{i['system_id']} analyte is not RNA")

        assert not msg, '\n'.join(msg)

    def validate_dna_pairs(self, targets, references):
        """Validate targets, references tuples for base dna pipelines."""
        self.validate_is_pair(targets, references)
        self.validate_dna_only(targets + references)

    def validate_same_technique(self, targets, references):
        """Validate targets and references have same bedfile."""
        ttec = {i['technique']['slug'] for i in targets}
        rtec = {i['technique']['slug'] for i in references}
        assert len(rtec) == 1, f'Expected one technique, got: {rtec}'
        assert len(ttec) == 1, f'Expected one technique, got: {ttec}'
        assert rtec == ttec, f'Same techniques required: {ttec}, {rtec}'

    def validate_species(self, workflows):
        """Validate workflows's species is same as pipeline's setting."""
        msg = []

        for i in workflows:
            if i['specimen']['individual']['species'] != self.SPECIES:
                msg.append(f'{i["system_id"]} species not supported')

        assert not msg, '\n'.join(msg)
