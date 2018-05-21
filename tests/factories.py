import string
import random

from factory import fuzzy
import factory

class BaseFactory(factory.DictFactory):
    created_by = 'admin'
    notes = fuzzy.FuzzyText(length=100, chars=string.ascii_letters)
    tags = [{'name': 'tag 1'}, {'name': 'tag 2'}]


class CenterFactory(BaseFactory):
    acronym = fuzzy.FuzzyText(length=3, chars=string.ascii_letters)
    internal = fuzzy.FuzzyChoice([True, False])
    name = fuzzy.FuzzyText(length=6, chars=string.ascii_letters)


class PlatformFactory(BaseFactory):
    manufacturer = fuzzy.FuzzyText(length=3, chars=string.ascii_letters)
    system = fuzzy.FuzzyText(length=3, chars=string.ascii_letters)
    version = fuzzy.FuzzyText(length=3, chars=string.ascii_letters)


class DiseaseFactory(BaseFactory):
    acronym = factory.Sequence(lambda n: f'DISEASE-{n}')
    name = fuzzy.FuzzyText(length=6, chars=string.ascii_letters)


class PipelineFactory(BaseFactory):
    description = fuzzy.FuzzyText(length=10, chars=string.ascii_letters)
    name = fuzzy.FuzzyText(length=10, chars=string.ascii_letters)
    version = fuzzy.FuzzyText(length=10, chars=string.ascii_letters)


class AnalysisFactory(BaseFactory):
    pipeline = factory.SubFactory(PipelineFactory)


class TechniqueFactory(BaseFactory):
    analyte = fuzzy.FuzzyChoice(['DNA', 'RNA'])
    name = fuzzy.FuzzyText(length=12, chars=string.hexdigits)

    @factory.lazy_attribute
    def method(self):
        if self.analyte == 'RNA':
            return random.choice(['TR', 'WT'])
        return random.choice(['CS', 'TD', 'WE', 'WG', 'MD'])


class IndividualFactory(BaseFactory):
    birth_year = fuzzy.FuzzyInteger(1800, 2100)
    center = factory.SubFactory(CenterFactory)
    gender = fuzzy.FuzzyChoice(['MALE', 'FEMALE', 'UNKNOWN'])
    species = fuzzy.FuzzyChoice(['HUMAN', 'MOUSE'])
    research_id = fuzzy.FuzzyText(length=12, chars=string.hexdigits)


class SpecimenFactory(BaseFactory):
    disease = factory.SubFactory(DiseaseFactory)
    individual = factory.SubFactory(IndividualFactory)
    pdx_id = fuzzy.FuzzyText(length=12, chars=string.hexdigits)
    specimen_class = fuzzy.FuzzyChoice(['TUMOR', 'NORMAL'])
    research_id = fuzzy.FuzzyText(length=12, chars=string.hexdigits)


class WorkflowFactory(BaseFactory):
    cell_type = fuzzy.FuzzyChoice(['SINGLE', 'BULK'])
    center_id = fuzzy.FuzzyText(length=12, chars=string.hexdigits)
    portion_id = fuzzy.FuzzyText(length=12, chars=string.hexdigits)
    raw_data = fuzzy.FuzzyChoice(['BAM', 'FASTQ'])
    read_length = fuzzy.FuzzyChoice(["100", "150"])
    read_type = fuzzy.FuzzyChoice(['PAIR-END', 'SINGLE-END'])
    sequencing_center = factory.SubFactory(CenterFactory)
    sequencing_platform = factory.SubFactory(PlatformFactory)
    specimen = factory.SubFactory(SpecimenFactory)
    technique = factory.SubFactory(TechniqueFactory)
    research_id = fuzzy.FuzzyText(length=12, chars=string.hexdigits)
