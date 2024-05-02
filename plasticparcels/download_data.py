import os
from pathlib import Path
from typing import List
from urllib.request import urlretrieve

import platformdirs

__all__ = ["download_plasticparcels_dataset", "get_data_home", "list_plasticparcels_datasets"]

plasticparcels_data_files = {
    "NEMO0083": [
        (("release_maps", "coastal"), "coastal_population_MPW_NEMO0083.csv"),
        (("release_maps", "rivers"), "river_emissions_NEMO0083.csv"),
        (("release_maps", "fisheries"), "agg_data_fisheries_info.csv"),
        (("unbeaching", "filename"), "land_current_NEMO0083.nc"),
    ],
}


plasticparcels_data_url = "https://plasticadrift.science.uu.nl/plasticparcels/"


def get_data_home(data_home=None):
    """Return a path to the cache directory for example datasets.

    This directory is used by :func:`load_dataset`.

    If the ``data_home`` argument is not provided, it will use a directory
    specified by the ``PLASTICPARCELS_DATA`` environment variable (if it exists)
    or otherwise default to an OS-appropriate user cache location.
    """
    if data_home is None:
        data_home = os.environ.get("PLASTICPARCELS_DATA", platformdirs.user_cache_dir("plasticsparcels"))
    data_home = os.path.expanduser(data_home)
    if not os.path.exists(data_home):
        os.makedirs(data_home)
    return data_home


def list_plasticparcels_datasets() -> List[str]:
    """List the available example datasets.

    Use :func:`download_example_dataset` to download one of the datasets.

    Returns
    -------
    datasets : list of str
        The names of the available example datasets.
    """
    return list(plasticparcels_data_files.keys())


def download_plasticparcels_dataset(dataset: str, settings, data_home=None):
    """Load an example dataset from the parcels website.

    This function provides quick access to a small number of example datasets
    that are useful in documentation and testing in parcels.

    Parameters
    ----------
    dataset : str
        Name of the dataset to load.
    data_home : pathlike, optional
        The directory in which to cache data. If not specified, the value
        of the ``PLASTICPARCELS_DATA`` environment variable, if any, is used.
        Otherwise the default location is assigned by :func:`get_data_home`.

    Returns
    -------
    dataset_folder : Path
        Path to the folder containing the downloaded dataset files.
    """
    # Dev note: `dataset` is assumed to be a folder name with netcdf files
    if dataset not in plasticparcels_data_files:
        raise ValueError(
            f"Dataset {dataset!r} not found. Available datasets are: "
            ", ".join(plasticparcels_data_files.keys())
        )

    cache_folder = get_data_home(data_home)
    dataset_folder = Path(cache_folder) / dataset

    if not dataset_folder.exists():
        dataset_folder.mkdir(parents=True)

    for settings_path, filename in plasticparcels_data_files[dataset]:
        filepath = dataset_folder / filename
        settings[settings_path[0]][settings_path[1]] = filepath
        if not filepath.exists():
            url = f"{plasticparcels_data_url}/{dataset}/{filename}"
            urlretrieve(url, str(filepath))

    return settings
