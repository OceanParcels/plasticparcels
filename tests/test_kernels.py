import parcels.application_kernels
import plasticparcels as pp
import parcels
import numpy as np
from datetime import datetime, timedelta
import pytest


def make_standard_simulation_settings():
    simulation_settings = {'startdate': datetime.strptime('2020-01-04-00:00:00', '%Y-%m-%d-%H:%M:%S'),
                           'runtime': timedelta(days=2),
                           'outputdt': timedelta(hours=1),
                           'dt': timedelta(minutes=20),
                           }
    return simulation_settings


def make_standard_plastictype_settings():
    # Use tiny wind percentage because test data set is not large and wind speeds are quick!
    plastictype_settings = {'wind_coefficient': 0.0001,     # Percentage of wind to apply to particles
                            'plastic_diameter': 0.001,  # Plastic particle diameter (m)
                            'plastic_density': 1030.,   # Plastic particle density (kg/m^3)
                            }
    return plastictype_settings


def make_standard_particleset(fieldset, settings):
    # Generate a particleset that has particles in the test domain
    release_locations = {'lons': [18, 18.25, 18.5], 'lats': [35, 35, 35],
                         'plastic_amount': [1, 1, 1]}
    pset = pp.constructors.create_particleset(fieldset, settings, release_locations)

    return pset


@pytest.mark.parametrize('use_3D', [True, False])
def test_advection_only(use_3D):
    settings_file = 'tests/test_data/test_settings.json'
    settings = pp.utils.load_settings(settings_file)

    settings['simulation'] = make_standard_simulation_settings()
    settings['plastictype'] = make_standard_plastictype_settings()

    # Turn on/off 3D advection
    settings['use_3D'] = use_3D

    # Turn off all other processes
    settings['use_biofouling'] = False
    settings['use_stokes'] = False
    settings['use_wind'] = False
    settings['use_mixing'] = False

    fieldset = pp.constructors.create_fieldset(settings)
    if use_3D:
        kernels = [parcels.application_kernels.AdvectionRK4_3D, pp.kernels.checkThroughBathymetry,
                   pp.kernels.checkErrorThroughSurface, pp.kernels.deleteParticle]
    else:
        kernels = [parcels.application_kernels.AdvectionRK4, pp.kernels.deleteParticle]

    pset = make_standard_particleset(fieldset, settings)

    start_lons = pset.lon.copy()
    start_lats = pset.lat.copy()

    pset.execute(kernels, runtime=settings['simulation']['runtime'], dt=settings['simulation']['dt'])

    # Assert that the particles move from their initial location
    assert (np.sum(np.abs(pset.lon - start_lons)) > 0.) & (np.sum(np.abs(pset.lat - start_lats)) > 0.)


def test_settling_velocity():
    settings_file = 'tests/test_data/test_settings.json'
    settings = pp.utils.load_settings(settings_file)

    settings['simulation'] = make_standard_simulation_settings()
    settings['plastictype'] = make_standard_plastictype_settings()

    # Turn on 3D advection
    settings['use_3D'] = True

    # Turn off all other processes
    settings['use_biofouling'] = False
    settings['use_stokes'] = False
    settings['use_wind'] = False
    settings['use_mixing'] = False

    fieldset = pp.constructors.create_fieldset(settings)

    kernels = pp.constructors.create_kernel(fieldset)

    pset = make_standard_particleset(fieldset, settings)

    start_lons = pset.lon.copy()
    start_lats = pset.lat.copy()
    start_depths = pset.depth.copy()

    pset.execute(kernels, runtime=settings['simulation']['runtime'], dt=settings['simulation']['dt'])

    # Assert that the particles move from their initial location
    assert (np.sum(np.abs(pset.lon - start_lons)) > 0.) & (np.sum(np.abs(pset.lat - start_lats)) > 0.) & (np.sum(np.abs(pset.depth - start_depths)) > 0.)


def test_biofouling():
    settings_file = 'tests/test_data/test_settings.json'
    settings = pp.utils.load_settings(settings_file)

    settings['simulation'] = make_standard_simulation_settings()
    settings['plastictype'] = make_standard_plastictype_settings()

    # Turn on biofouling
    settings['use_biofouling'] = True

    # Turn off all other processes
    settings['use_3D'] = False
    settings['use_stokes'] = False
    settings['use_wind'] = False
    settings['use_mixing'] = False

    fieldset = pp.constructors.create_fieldset(settings)

    kernels = [pp.kernels.PolyTEOS10_bsq, pp.kernels.Biofouling,
               pp.kernels.checkThroughBathymetry, pp.kernels.checkErrorThroughSurface,
               pp.kernels.deleteParticle]

    pset = make_standard_particleset(fieldset, settings)

    start_depths = pset.depth.copy()

    pset.execute(kernels, runtime=settings['simulation']['runtime'], dt=settings['simulation']['dt'])

    # Assert that the particles move from their initial location
    assert (np.sum(np.abs(pset.depth - start_depths)) > 0.)


def test_Stokes():
    settings_file = 'tests/test_data/test_settings.json'
    settings = pp.utils.load_settings(settings_file)

    # Required for the unbeaching kernel
    settings = pp.utils.download_plasticparcels_dataset('NEMO0083', settings, 'input_data')

    settings['simulation'] = make_standard_simulation_settings()
    settings['plastictype'] = make_standard_plastictype_settings()

    # Turn on Stokes Drift
    settings['use_stokes'] = True

    # Turn off all other processes
    settings['use_3D'] = False
    settings['use_biofouling'] = False
    settings['use_wind'] = False
    settings['use_mixing'] = False

    fieldset = pp.constructors.create_fieldset(settings)

    kernels = [pp.kernels.StokesDrift, pp.kernels.unbeaching, pp.kernels.periodicBC, pp.kernels.deleteParticle]

    pset = make_standard_particleset(fieldset, settings)

    start_lons = pset.lon.copy()
    start_lats = pset.lat.copy()
    print(pset[0].dt)
    pset.execute(kernels, runtime=settings['simulation']['runtime'], dt=settings['simulation']['dt'])

    # Assert that the particles move from their initial location
    assert (np.sum(np.abs(pset.lon - start_lons)) > 0.) & (np.sum(np.abs(pset.lat - start_lats)) > 0.)


def test_wind():
    settings_file = 'tests/test_data/test_settings.json'
    settings = pp.utils.load_settings(settings_file)

    # Required for the unbeaching kernel
    settings = pp.utils.download_plasticparcels_dataset('NEMO0083', settings, 'input_data')

    settings['simulation'] = make_standard_simulation_settings()
    settings['plastictype'] = make_standard_plastictype_settings()

    # Turn on Stokes Drift
    settings['use_wind'] = True

    # Turn off all other processes
    settings['use_3D'] = False
    settings['use_biofouling'] = False
    settings['use_stokes'] = False
    settings['use_mixing'] = False

    fieldset = pp.constructors.create_fieldset(settings)

    kernels = [pp.kernels.WindageDrift, pp.kernels.unbeaching, pp.kernels.periodicBC, pp.kernels.deleteParticle]

    pset = make_standard_particleset(fieldset, settings)

    start_lons = pset.lon.copy()
    start_lats = pset.lat.copy()

    pset.execute(kernels, runtime=settings['simulation']['runtime'], dt=settings['simulation']['dt'])

    # Assert that the particles move from their initial location
    assert (np.sum(np.abs(pset.lon - start_lons)) > 0.) & (np.sum(np.abs(pset.lat - start_lats)) > 0.)


def test_mixing():
    settings_file = 'tests/test_data/test_settings.json'
    settings = pp.utils.load_settings(settings_file)

    settings['simulation'] = make_standard_simulation_settings()
    settings['plastictype'] = make_standard_plastictype_settings()

    settings['use_3D'] = True
    settings['use_mixing'] = True

    # Turn off all other processes
    settings['use_wind'] = False
    settings['use_biofouling'] = False
    settings['use_stokes'] = False

    fieldset = pp.constructors.create_fieldset(settings)
    kernels = [parcels.application_kernels.AdvectionRK4_3D, pp.kernels.checkThroughBathymetry,
               pp.kernels.checkErrorThroughSurface, pp.kernels.deleteParticle]

    kernels_mixing = [parcels.application_kernels.AdvectionRK4_3D, pp.kernels.VerticalMixing,
                      pp.kernels.checkThroughBathymetry, pp.kernels.checkErrorThroughSurface,
                      pp.kernels.deleteParticle]

    pset = make_standard_particleset(fieldset, settings)
    pset_mixing = make_standard_particleset(fieldset, settings)

    pset.execute(kernels, runtime=settings['simulation']['runtime'], dt=settings['simulation']['dt'])
    pset.execute(kernels_mixing, runtime=settings['simulation']['runtime'], dt=settings['simulation']['dt'])

    # Assert that the particles move from their initial location
    assert (np.sum(np.abs(pset.lon - pset_mixing.lon)) > 0.) & (np.sum(np.abs(pset.lat - pset_mixing.lat)) > 0.)
