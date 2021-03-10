# ---------------------------------------------------------------------------------------------
# MODULES
# ---------------------------------------------------------------------------------------------
import netCDF4 as nc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import os
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import LogNorm, from_levels_and_colors, ListedColormap, LinearSegmentedColormap
from matplotlib import cm

# ---------------------------------------------------------------------------------------------
# LOAD DATA
# ---------------------------------------------------------------------------------------------
ds_lpi = Dataset('/staging/leuven/stg_00024/OUTPUT/michelb/nu-wrf-dev/Great_Slave_Lake/2015/LPI.nc','r')
time = ds_lpi['times'][:]
lat = ds_lpi['lat'][:]
lon = ds_lpi['lon'][:]
lpi = ds_lpi['lpi'][:]

print('done')