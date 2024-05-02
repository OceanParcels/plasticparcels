import plasticparcels as pp


def test_settings():
    settings = pp.utils.load_settings('docs/examples/default_settings.json')
    assert settings is not None
