#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 30 17:27:39 2021

@author: kaandorp
"""
import numpy as np
from glob import glob
import pandas as pd
from datetime import timedelta
import os
import math

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
    
    possible_files = sorted( glob(dirread + string_ % (yr0-1)) + glob(dirread + string_ % yr0) + glob(dirread + string_ % (yr0+1)) )
    
    indices_use = []
    for i1,file_ in enumerate(possible_files):
        str_date = os.path.basename(file_)[i_date_s:i_date_e]
        date_ = pd.Timestamp(str_date) 

        if date_ > time_start and date_ < time_end:
            indices_use.append(i1)
            
    files_use = list( possible_files[i] for i in indices_use)
    return files_use


# Added 2023-08-31
def create_directory(directory):
    if not os.path.exists(directory):
        print('Creating directory %s' % directory)
        os.makedirs(directory)


