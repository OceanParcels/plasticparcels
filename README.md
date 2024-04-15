# PlasticParcels

`PlasticParcels` is a python tool for simulating the dispersion of plastics in the global ocean. The tool is based on the `Parcels` computational Lagrangian framework [@Lange2017, @Delandmeter2019], providing a modular and customizable python notebook for advecting virtual plastic particles with a wide range of properties, applying a collection of physical processes such as Stokes drift, wind-induced drift, biofouling, and turbulent mixing. Additionally, `PlasticParcels` provides particle initialisation maps that represent the best estimates for plastic pollution emmissions along coastlines [@Jambeck2015], from river sources [@Meijer2021], and in the open-ocean from fishing-related activities [@Kroodsma2018].

## Installation

The latest version of `PlasticParcels` can be installed directly from github, via:
```
git clone https://github.com/OceanParcels/PlasticParcels.git
cd PlasticParcels; pip install -r requirements.txt
python PlasticParcels/run_initialisation.py
export PYTHONPATH="$PYTHONPATH:$PWD"
```
Ensure that you have updated `settings.txt` with the required directories and filenames.
