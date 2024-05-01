import numpy as np
from glob import glob
import pandas as pd
from datetime import timedelta
import os
import math
import shapely as sh
import json

def getclosest_ij(lats, lons, latpt, lonpt):
    """Function to find the index of the closest point to a certain lon/lat value."""
    dist_sq = (lats - latpt) ** 2 + (lons - lonpt) ** 2  # find squared distance of every point on grid
    minindex_flattened = np.nanargmin(dist_sq)  # 1D index of minimum dist_sq element
    return np.unravel_index(minindex_flattened,
                            lats.shape)  # Get 2D index for latvals and lonvals arrays from 1D index

    
def select_files(dirread,string_,date_current,dt_sim,i_date_s=-13,i_date_e=-3, dt_margin=8):
    # set dt_margin to e.g. 32 when dealing with monthly data, or 8 when dealing with weekly data
    yr0 = date_current.year
    
    time_start = date_current - timedelta(days=dt_margin)
    time_end = date_current + timedelta(days=dt_sim) + timedelta(days=dt_margin)
    yrEnd = time_end.year

    ragged_files = [glob(dirread + string_ % (yr0+i)) for i in range(-1,yrEnd-yr0+1)]
    possible_files = []

    for i in range(len(ragged_files)):
        for file_ in ragged_files[i]:
            possible_files.append(file_)
    
    possible_files = sorted(possible_files)
    
    indices_use = []
    for i1,file_ in enumerate(possible_files):
        str_date = os.path.basename(file_)[i_date_s:i_date_e]
        date_ = pd.Timestamp(str_date) 

        if date_ > time_start and date_ < time_end:
            indices_use.append(i1)
            
    files_use = list( possible_files[i] for i in indices_use)
    return files_use


def create_directory(directory):
    if not os.path.exists(directory):
        print('Creating directory %s' % directory)
        os.makedirs(directory)


def distance(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance between two points
    on the earth (specified in decimal degrees)
    """

    #Convert decimal degrees to Radians:
    lon1r = np.radians(lon1)
    lat1r = np.radians(lat1)
    lon2r = np.radians(lon2)
    lat2r = np.radians(lat2)

    #Implementing Haversine Formula:
    dlon = np.subtract(lon2r, lon1r)
    dlat = np.subtract(lat2r, lat1r)

    a = np.add(np.power(np.sin(np.divide(dlat, 2)), 2),
                          np.multiply(np.cos(lat1r),
                                      np.multiply(np.cos(lat2r),
                                                  np.power(np.sin(np.divide(dlon, 2)), 2))))

    c = np.multiply(2, np.arcsin(np.sqrt(a)))
    r = 6371

    return c*r

def get_coords_from_polygon(shape):
    """
    Get a list of coordinate points on a Polygon
    (or MultiPolygon) shape

    Based on: https://stackoverflow.com/questions/58844463/how-to-get-a-list-of-every-point-inside-a-multipolygon-using-shapely
    """
    coords = []

    if isinstance(shape, sh.geometry.Polygon):
        coords.append(np.asarray(shape.exterior.coords[:-1]))
        for linearring in shape.interiors:
            coords.append(np.asarray(linearring.coords[:-1]))
    elif isinstance(shape, sh.geometry.MultiPolygon):
        #polygons = list(shape)
        for polygon in shape.geoms:
            coords.append(get_coords_from_polygon(polygon))
    coords = np.concatenate(coords)
    return coords

def load_settings(filename):
    """ A function to load a settings file in json format"""
    with open(filename, "r") as file:
        settings = json.load(file)
    return settings