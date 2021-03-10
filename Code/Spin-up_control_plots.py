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
ds_spinup = Dataset('/scratch/leuven/336/vsc33651/nu-wrf-dev/lis_spinup/LIS_RST_NOAHMP36d01.nc','r')
month = ds_spinup['time'][:]
SSTC = ds_spinup['SSTC'][:,:,:]
SH2O = ds_spinup['SH2O'][:,:,:]
SSTC[SSTC == 0] = np.nan
SH2O[SH2O == 0] = np.nan
SSTC_pl = np.nanmean(SSTC, axis=1)
SH2O_pl = np.nanmean(SH2O, axis=1)

# ---------------------------------------------------------------------------------------------
# TEMPORAL PLOTS
# ---------------------------------------------------------------------------------------------
# SSTC
plt.plot(month, SSTC_pl[:,0], 'k')
plt.xlabel('Months since 9-2005', size = 12)
plt.ylabel('snow/soil temperature (K)', size = 12)
plt.title('SSTC layer 0')
plt.ticklabel_format(axis='y',style='sci', scilimits=(0,0))
plt.show()
# SH2O
plt.plot(month, SH2O_pl[:,0], 'k')
plt.xlabel('Months since 9-2005', size = 12)
plt.ylabel('volumetric liquid soil moisture (m³ m$^{-3}$)', size = 12)
plt.title('SH2O layer 0')
plt.ticklabel_format(axis='y',style='sci', scilimits=(0,0))
plt.show()
