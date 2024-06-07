import plasticparcels as pp


def test_settings():
    settings = pp.utils.load_settings('tests/test_data/test_settings.json')
    assert settings is not None


def test_download_script(tmpdir):
    settings = pp.utils.load_settings('tests/test_data/test_settings.json')
    pp.utils.download_plasticparcels_dataset("NEMO0083", settings=settings, data_home=tmpdir)
    assert True
