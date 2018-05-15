import string

from factory import fuzzy
import factory


class DiseaseFactory(factory.DictFactory):
    acronym = factory.Sequence(lambda n: f'DISEASE-{n}')
    name = fuzzy.FuzzyText(length=6, chars=string.ascii_letters)
    created_by = 'admin'
