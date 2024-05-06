import plasticparcels as pp
import parcels
import pytest

def make_simple_fieldset():
    fieldset = parcels.FieldSet.from_data({'U': 0, 'V': 0}, {'lon': 0, 'lat': 0}, mesh='spherical')
    return fieldset


@pytest.mark.parametrize('maptype', ['coastal', 'rivers', 'fisheries', 'global_concentrations'])
def test_maptypes(maptype):
    settings = pp.utils.load_settings('docs/examples/example_Italy_coast_settings.json')
    for key in settings["release_maps"]:
        settings["release_maps"][key] = settings["release_maps"][key].replace("input_data", "docs/examples/input_data")
    fieldset = make_simple_fieldset()

    settings['release'] = {
        'initialisation_type': maptype,
        'country': 'Italy',
    }

    settings['plastictype'] = {  # TODO these should be defaults?
        'wind_coefficient' : 0.01, # Percentage of wind to apply to particles
        'plastic_diameter' : 0.001, # Plastic particle diameter (m)
        'plastic_density' : 1030., # Plastic particle density (kg/m^3)
    }

    pset = pp.constructors.create_particleset_from_map(fieldset, settings)
