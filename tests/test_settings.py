from os import environ
from cli import system_settings

def test_assert_read_settings_from_yaml(tmpdir):
    settings_yml = tmpdir.join('test.yml')
    settings_yml.write(f'BASE_STORAGE_DIRECTORY: {tmpdir.strpath}')
    environ['CLI_SETTINGS_PATH'] = settings_yml.strpath
    assert system_settings.BASE_STORAGE_DIRECTORY == tmpdir.strpath
