import plasticparcels as pp
import parcels
from datetime import datetime, timedelta

def test_create_fieldset():
    settings_file = 'tests/test_data/test_settings.json'
    settings = pp.utils.load_settings(settings_file)
    settings = pp.utils.download_plasticparcels_dataset('NEMO0083', settings, 'input_data')

    settings['simulation'] = {'startdate': datetime.strptime('2020-01-04-00:00:00', '%Y-%m-%d-%H:%M:%S'), # Start date of simulation
                                'runtime': timedelta(days=2),        # Runtime of simulation
                                'outputdt': timedelta(hours=1),      # Timestep of output
                                'dt': timedelta(minutes=20),          # Timestep of advection
                                }
    
    fieldset = pp.constructors.create_fieldset(settings)

    assert isinstance(fieldset, parcels.FieldSet)

