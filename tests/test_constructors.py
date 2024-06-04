import parcels.application_kernels
import plasticparcels as pp
import parcels
from datetime import datetime, timedelta
import pytest


def make_simple_fieldset():
    fieldset = parcels.FieldSet.from_data({'U': 0, 'V': 0},
                                          {'lon': 0, 'lat': 0},
                                          mesh='spherical')
    return fieldset


# Test the create_hydrodynamic_fieldset() function
def test_create_hydrodynamic_fieldset():
    settings_file = 'tests/test_data/test_settings.json'
    settings = pp.utils.load_settings(settings_file)
    settings['simulation'] = {'startdate': datetime.strptime('2020-01-04-00:00:00', '%Y-%m-%d-%H:%M:%S'),
                                'runtime': timedelta(days=1),
                                'outputdt': timedelta(hours=1),
                                'dt': timedelta(minutes=20),
                                }

    fieldset = pp.constructors.create_hydrodynamic_fieldset(settings)

    assert isinstance(fieldset, parcels.FieldSet)

# Test the create_fieldset() function
def test_create_fieldset():
    settings_file = 'tests/test_data/test_settings.json'
    settings = pp.utils.load_settings(settings_file)
    settings = pp.utils.download_plasticparcels_dataset('NEMO0083', settings, 'input_data')

    settings['simulation'] = {'startdate': datetime.strptime('2020-01-04-00:00:00', '%Y-%m-%d-%H:%M:%S'),
                                'runtime': timedelta(days=1),
                                'outputdt': timedelta(hours=1),
                                'dt': timedelta(minutes=20),
                                }
    
    fieldset = pp.constructors.create_fieldset(settings)

    assert isinstance(fieldset, parcels.FieldSet)

# Test three different initialisation maps
@pytest.mark.parametrize('initialisation_map', ['coastal', 'fisheries', 'rivers'])
def test_create_particleset_from_map(initialisation_map):
    settings_file = 'tests/test_data/test_settings.json'
    settings = pp.utils.load_settings(settings_file)
    settings = pp.utils.download_plasticparcels_dataset('NEMO0083', settings, 'input_data')

    settings['simulation'] = {'startdate': datetime.strptime('2020-01-04-00:00:00', '%Y-%m-%d-%H:%M:%S'),
                                'runtime': timedelta(days=1),
                                'outputdt': timedelta(hours=1),
                                'dt': timedelta(minutes=20),
                                }
 
    settings['plastictype'] = {'wind_coefficient' : 0.01,  # Percentage of wind to apply to particles
                                'plastic_diameter' : 0.001, # Plastic particle diameter (m)
                                'plastic_density' : 1030.,  # Plastic particle density (kg/m^3)
                                }
    
    settings['release'] = {'initialisation_type': initialisation_map,
                            'country': 'Italy',
                            }
    
    fieldset = make_simple_fieldset()
    pset = pp.constructors.create_particleset_from_map(fieldset, settings)

    assert isinstance(pset, parcels.ParticleSet)

# Test the two global concentration release map types
@pytest.mark.parametrize('concentration_type', ['Beach', 'Ocean'])
def test_create_particleset_from_map_concentrations(concentration_type):
    settings_file = 'tests/test_data/test_settings.json'
    settings = pp.utils.load_settings(settings_file)
    settings = pp.utils.download_plasticparcels_dataset('NEMO0083', settings, 'input_data')

    settings['simulation'] = {'startdate': datetime.strptime('2020-01-04-00:00:00', '%Y-%m-%d-%H:%M:%S'),
                                'runtime': timedelta(days=1),
                                'outputdt': timedelta(hours=1),
                                'dt': timedelta(minutes=20),
                                }
 
    settings['plastictype'] = {'wind_coefficient' : 0.01,  # Percentage of wind to apply to particles
                                'plastic_diameter' : 0.001, # Plastic particle diameter (m)
                                'plastic_density' : 1030.,  # Plastic particle density (kg/m^3)
                                }
    
    settings['release'] = {'initialisation_type': 'global_concentrations',
                            'concentration_type': concentration_type,
                            }
    
    fieldset = make_simple_fieldset()
    pset = pp.constructors.create_particleset_from_map(fieldset, settings)

    assert isinstance(pset, parcels.ParticleSet)

# Test three different initialisation maps
@pytest.mark.parametrize('use_3D', [True, False])
@pytest.mark.parametrize('use_biofouling', [True, False])
@pytest.mark.parametrize('use_stokes', [True, False])
@pytest.mark.parametrize('use_wind', [True, False])
@pytest.mark.parametrize('use_mixing', [True, False])
def test_create_kernel(use_3D, use_biofouling, use_stokes, use_wind, use_mixing):
    settings_file = 'tests/test_data/test_settings.json'
    settings = pp.utils.load_settings(settings_file)
    settings = pp.utils.download_plasticparcels_dataset('NEMO0083', settings, 'input_data')
    settings['simulation'] = {'startdate': datetime.strptime('2020-01-04-00:00:00', '%Y-%m-%d-%H:%M:%S'),
                              'runtime': timedelta(days=1)}

    settings['use_3D'] = use_3D
    settings['use_biofouling'] = use_biofouling
    settings['use_stokes'] = use_stokes
    settings['use_wind'] = use_wind
    settings['use_mixing'] = use_mixing
    
    fieldset = pp.constructors.create_fieldset(settings)
    kernels = pp.constructors.create_kernel(fieldset)

    if use_3D:
        assert parcels.application_kernels.AdvectionRK4_3D in kernels
        assert pp.kernels.checkThroughBathymetry in kernels
        assert pp.kernels.checkErrorThroughSurface in kernels
    else:
        assert parcels.application_kernels.AdvectionRK4 in kernels

    if use_3D and use_biofouling:
        assert pp.kernels.Biofouling in kernels
    elif use_3D and not use_biofouling:
        assert pp.kernels.SettlingVelocity in kernels
    
    if use_mixing:
        assert pp.kernels.VerticalMixing in kernels

    if use_stokes:
        assert pp.kernels.StokesDrift in kernels
        assert pp.kernels.unbeaching in kernels
    
    if use_wind:
        assert pp.kernels.WindageDrift in kernels
        assert pp.kernels.unbeaching in kernels
    
    assert pp.kernels.deleteParticle in kernels
    assert pp.kernels.periodicBC in kernels
    assert pp.kernels.PolyTEOS10_bsq in kernels