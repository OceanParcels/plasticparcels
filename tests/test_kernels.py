import plasticparcels as pp
import parcels


def make_simple_fieldset(settings):
    fieldset = parcels.FieldSet.from_data({'U': 0, 'V': 0}, {'lon': 0, 'lat': 0}, mesh='spherical')
    fieldset.use_3D = settings['use_3D']
    fieldset.use_biofouling = settings['use_biofouling']
    fieldset.use_stokes = settings['use_stokes']
    fieldset.use_wind = settings['use_wind']
    fieldset.use_mixing = settings['use_mixing']
    return fieldset


def test_kernels_options():
    settings = pp.utils.load_settings('docs/examples/example_Italy_coast_settings.json')
    settings['use_3D'] = False
    settings['use_biofouling'] = False
    settings['use_stokes'] = True
    settings['use_wind'] = True
    settings['use_mixing'] = True
    fieldset = make_simple_fieldset(settings)
    kernels = pp.constructors.create_kernel(fieldset)

    assert len(kernels) > 0
