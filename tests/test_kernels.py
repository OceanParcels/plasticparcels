import plasticparcels as pp
import parcels

def make_simple_fieldset():
    fieldset = parcels.FieldSet.from_data({'U': 0, 'V': 0}, {'lon': 0, 'lat': 0}, mesh='spherical')
    return fieldset

def test_kernels_options():
    settings = pp.utils.load_settings('docs/examples/example_Italy_coast_settings.json')
    fieldset = make_simple_fieldset()
