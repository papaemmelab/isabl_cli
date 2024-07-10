"""Instance factories used for testing."""

import string

from factory import fuzzy
import factory


class BaseFactory(factory.DictFactory):
    notes = fuzzy.FuzzyText(length=100, chars=string.ascii_letters)
    tags = [{"name": "tag 1"}, {"name": "tag 2"}]


class UserFactory(BaseFactory):
    username = fuzzy.FuzzyText(length=6, chars=string.ascii_letters)


class ProjectFactory(BaseFactory):
    analyst = factory.Sequence(lambda n: f"analyst-{n}@test.com")
    coordinator = factory.Sequence(lambda n: f"coordinator-{n}@test.com")
    description = fuzzy.FuzzyText(length=20, chars=string.ascii_lowercase + " ")
    owner = factory.Sequence(lambda n: f"owner-{n}@test.com")
    principal_investigator = factory.Sequence(lambda n: f"pi-{n}@test.com")
    title = fuzzy.FuzzyText(length=20, chars=string.ascii_lowercase + " ")
    sharing = {"can_share": [], "can_read": [], "is_public": True}


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
    reference_data = factory.SubFactory(factory.DictFactory)
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
    category = fuzzy.FuzzyChoice(["TUMOR", "NORMAL"])
    identifier = fuzzy.FuzzyText(length=12, chars=string.hexdigits)


class ExperimentFactory(BaseFactory):
    bam_files = factory.SubFactory(factory.DictFactory)
    aliquot_id = None
    identifier = fuzzy.FuzzyText(length=12, chars=string.hexdigits)
    center = factory.SubFactory(CenterFactory)
    raw_data = None
    platform = factory.SubFactory(PlatformFactory)
    sample = factory.SubFactory(SampleFactory)
    technique = factory.SubFactory(TechniqueFactory)

    @factory.lazy_attribute
    def projects(self):  # pylint: disable=no-self-use
        """Lazy method to return list of projects."""
        return [ProjectFactory()]
