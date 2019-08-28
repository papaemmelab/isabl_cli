"""Instance factories used for testing."""

import random
import string

from factory import fuzzy
import factory


class BaseFactory(factory.DictFactory):
    notes = fuzzy.FuzzyText(length=100, chars=string.ascii_letters)
    tags = [{"name": "tag 1"}, {"name": "tag 2"}]


class ProjectFactory(BaseFactory):
    analyst = factory.Sequence(lambda n: f"analyst-{n}@test.com")
    coordinator = factory.Sequence(lambda n: f"coordinator-{n}@test.com")
    description = fuzzy.FuzzyText(length=20, chars=string.ascii_lowercase + " ")
    owner = factory.Sequence(lambda n: f"owner-{n}@test.com")
    principal_investigator = factory.Sequence(lambda n: f"pi-{n}@test.com")
    title = fuzzy.FuzzyText(length=20, chars=string.ascii_lowercase + " ")


class CenterFactory(BaseFactory):
    acronym = fuzzy.FuzzyText(length=3, chars=string.ascii_letters)
    name = fuzzy.FuzzyText(length=6, chars=string.ascii_letters)


class PlatformFactory(BaseFactory):
    manufacturer = fuzzy.FuzzyText(length=3, chars=string.ascii_letters)
    system = fuzzy.FuzzyText(length=3, chars=string.ascii_letters)
    version = fuzzy.FuzzyText(length=3, chars=string.ascii_letters)


class DiseaseFactory(BaseFactory):
    acronym = factory.Sequence(lambda n: f"D-{n}")
    name = fuzzy.FuzzyText(length=6, chars=string.ascii_letters)


class AssemblyFactory(BaseFactory):
    name = "GRCh37"
    reference_data = {}
    species = fuzzy.FuzzyChoice(["HUMAN", "MOUSE"])


class ApplicationFactory(BaseFactory):
    description = fuzzy.FuzzyText(length=10, chars=string.ascii_letters)
    name = fuzzy.FuzzyText(length=10, chars=string.ascii_letters)
    version = fuzzy.FuzzyText(length=10, chars=string.ascii_letters)
    assembly = factory.SubFactory(AssemblyFactory)


class AnalysisFactory(BaseFactory):
    ran_by = "admin"
    application = factory.SubFactory(ApplicationFactory)
    targets = []
    references = []


class TechniqueFactory(BaseFactory):
    bed_files = factory.SubFactory(factory.DictFactory)
    name = fuzzy.FuzzyText(length=12, chars=string.hexdigits)
    method = fuzzy.FuzzyChoice(choices=["CS", "TD", "WE", "WG", "MD", "TR", "WT"])


class IndividualFactory(BaseFactory):
    birth_year = fuzzy.FuzzyInteger(1800, 2100)
    center = factory.SubFactory(CenterFactory)
    gender = fuzzy.FuzzyChoice(["MALE", "FEMALE", "UNKNOWN"])
    species = fuzzy.FuzzyChoice(["HUMAN", "MOUSE"])
    identifier = fuzzy.FuzzyText(length=12, chars=string.hexdigits)


class SampleFactory(BaseFactory):
    disease = factory.SubFactory(DiseaseFactory)
    individual = factory.SubFactory(IndividualFactory)
    pdx_id = fuzzy.FuzzyText(length=12, chars=string.hexdigits)
    sample_class = fuzzy.FuzzyChoice(["TUMOR", "NORMAL"])
    identifier = fuzzy.FuzzyText(length=12, chars=string.hexdigits)


class ExperimentFactory(BaseFactory):
    bam_files = factory.SubFactory(factory.DictFactory)
    cell_type = fuzzy.FuzzyChoice(["SINGLE", "BULK"])
    center_id = fuzzy.FuzzyText(length=12, chars=string.hexdigits)
    aliquot_id = None
    analyte_id = fuzzy.FuzzyText(length=12, chars=string.hexdigits)
    read_length = fuzzy.FuzzyChoice(["100", "150"])
    read_type = fuzzy.FuzzyChoice(["PAIR-END", "SINGLE-END"])
    identifier = fuzzy.FuzzyText(length=12, chars=string.hexdigits)
    generating_center = factory.SubFactory(CenterFactory)
    raw_data = None
    generating_platform = factory.SubFactory(PlatformFactory)
    sample = factory.SubFactory(SampleFactory)
    technique = factory.SubFactory(TechniqueFactory)

    @factory.lazy_attribute
    def projects(self):
        return [ProjectFactory()]
