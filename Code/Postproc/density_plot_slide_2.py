# Import packages
from scipy import stats
import pygrib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# from thesis_tools import *
import mpl_scatter_density
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.gridspec as gridspec
import matplotlib.axis as axis
from matplotlib.cm import get_cmap
from matplotlib.colors import from_levels_and_colors, ListedColormap, LinearSegmentedColormap
import xarray as xr
from numpy import meshgrid
from matplotlib import cm
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
from datetime import date, timedelta
import netCDF4 as nc
import datetime as dt
from datetime import timedelta
from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()
from netCDF4 import Dataset
from wrf import (getvar, latlon_coords)
import os
from matplotlib.colors import from_levels_and_colors, ListedColormap, LinearSegmentedColormap

# ---------------------------------------------------------------------------------------------
# OBSERVATIONS
# ---------------------------------------------------------------------------------------------
##### OBSERVATIONS
## Read in the observations (currently short version)
file = '/data/leuven/336/vsc33651/projects/nu-wrf-dev/Lightning_data/Slave_lake_full_0.1deg.nc'
ds = nc.Dataset(file, 'r')
print(ds)
## Extract data from NetCDF file
lats = ds.variables['lat'][:]
lons = ds.variables['lon'][:]
time = ds.variables['time'][:]
Flashdensity_CC = ds.variables['Flashdensity_CC'][:]
Flashdensity_CG = ds.variables['Flashdensity_CG'][:]

## Convert all times to datetimes
time_obs = pd.to_datetime(time)+timedelta(minutes=2) ## 2 minutes are added to get correct days.
## Without adding 1 minute, some data was 1 day too early (2 minutes). This is due to the conversion from datetime to int
## and back to datetime. Since the daily data is important, the exact minute of the day is not important.

## Read in the WRF files for d02 and d03 to extract coordinates
WRF_DIRECTORY = '/scratch/leuven/336/vsc33651/nu-wrf-dev/Great_Slave_Lake/noahmp36_modis_merra2'
WRF_FILES = ["wrfout_d02_2015-06-01_00:00:00", "wrfout_d03_2015-06-01_00:00:00"]

_WRF_FILES = [os.path.abspath(os.path.join(WRF_DIRECTORY,f)) for f in WRF_FILES]

## Check if the files exist
try:
    for f in _WRF_FILES:
        if not os.path.exists(f):
            raise ValueError("{} does not exist. "
                             "Check for typos or incorrect directory.".format(f))
except ValueError as e:
    raise e

# Create functions so that the WRF files only need to be specified using the WRF_FILES global above
def d02_wrf_file():
    global _WRF_FILES
    return _WRF_FILES[0]
def d03_wrf_file():
    global _WRF_FILES
    return _WRF_FILES[1]
print('All tests passed!')

file_path_d02 = d02_wrf_file()
d02_wrf_file = Dataset(file_path_d02)

file_path_d03 = d03_wrf_file()
d03_wrf_file = Dataset(file_path_d03)

## average daily Flashdensity for d02 and d03
# D02: Extract latitudes and longitudes from WRF output.
terrain_d02 = getvar(d02_wrf_file, "ter", timeidx=0)
lats_d02, lons_d02 = latlon_coords(terrain_d02)
min_lat_d02 = np.array(lats_d02.min(), dtype='f4')
max_lat_d02 = np.array(lats_d02.max(), dtype='f4')
min_lon_d02 = np.array(lons_d02.min(), dtype='f4')
max_lon_d02 = np.array(lons_d02.max(), dtype='f4')
lat_inds_d02 = np.where((lats >= min_lat_d02) & (lats <= max_lat_d02))[0]
lon_inds_d02 = np.where((lons >= min_lon_d02) & (lons <= max_lon_d02))[0]
FD_CC_D02 = ds.variables['Flashdensity_CC'][:,lat_inds_d02, lon_inds_d02]
FD_CG_D02 = ds.variables['Flashdensity_CG'][:,lat_inds_d02, lon_inds_d02]

# D03: Extract latitudes and longitudes from WRF output.
terrain_d03 = getvar(d03_wrf_file, "ter", timeidx=0)
lats_d03, lons_d03 = latlon_coords(terrain_d03)
min_lat_d03 = np.array(lats_d03.min(), dtype='f4')
max_lat_d03 = np.array(lats_d03.max(), dtype='f4')
min_lon_d03 = np.array(lons_d03.min(), dtype='f4')
max_lon_d03 = np.array(lons_d03.max(), dtype='f4')
lat_inds_d03 = np.where((lats >= min_lat_d03) & (lats <= max_lat_d03))[0]
lon_inds_d03 = np.where((lons >= min_lon_d03) & (lons <= max_lon_d03))[0]
FD_CC_D03 = ds.variables['Flashdensity_CC'][:,lat_inds_d03, lon_inds_d03]
FD_CG_D03 = ds.variables['Flashdensity_CG'][:,lat_inds_d03, lon_inds_d03]

## Group the total amount of flashes for each domain per day
FD_CC_D02_pd = FD_CC_D02.sum(axis=(1,2))
FD_CG_D02_pd = FD_CG_D02.sum(axis=(1,2))
FD_CC_D03_pd = FD_CC_D03.sum(axis=(1,2))
FD_CG_D03_pd = FD_CG_D03.sum(axis=(1,2))


## Calculate CG/(CC+CG)
CC_CG_d02 = np.add(FD_CC_D02_pd, FD_CG_D02_pd) # First calculate CC+CG for each domain
CC_CG_d03 = np.add(FD_CC_D03_pd, FD_CG_D03_pd) # First calculate CC+CG for each domain
CG_ratio_d02 = np.divide(FD_CG_D02_pd,CC_CG_d02)  # Then divide elementwise to get CG/(CC+CG)
CG_ratio_d03 = np.divide(FD_CG_D03_pd,CC_CG_d03)  # Then divide elementwise to get CG/(CC+CG)

## Calculate total amount of flashes per location per day
CC_CG_PL_d02 = np.add(FD_CC_D02, FD_CG_D02) # First calculate CC+CG for each domain
#CC_CG_d02 = np.where(np.isnan(CC_CG_d02_nan), 0, CC_CG_d02_nan)
CC_CG_PL_d03 = np.add(FD_CC_D03, FD_CG_D03) # First calculate CC+CG for each domain
#CC_CG_d03 = np.where(np.isnan(CC_CG_d03), 0, CC_CG_d03)

## Total flashes summer 2015
flashes_2015_PL_d02 = CC_CG_PL_d02[25:112,:,:]
flashes_2015_PL_d03 = CC_CG_PL_d03[25:112,:,:]
CC_CG_PL_d02_2015 = np.sum(flashes_2015_PL_d02, axis=0)
CC_CG_PL_d03_2015 = np.sum(flashes_2015_PL_d03, axis=0)
## Total flashes summer 2016
flashes_2016_PL_d02 = CC_CG_PL_d02[159:248,:,:]
flashes_2016_PL_d03 = CC_CG_PL_d03[159:248,:,:]
CC_CG_PL_d02_2016 = np.sum(flashes_2016_PL_d02, axis=0)
CC_CG_PL_d03_2016 = np.sum(flashes_2016_PL_d03, axis=0)
## Total flashes summer 2017
flashes_2017_PL_d02 = CC_CG_PL_d02[289:376,:,:]
flashes_2017_PL_d03 = CC_CG_PL_d03[289:376,:,:]
CC_CG_PL_d02_2017 = np.sum(flashes_2017_PL_d02, axis=0)
CC_CG_PL_d03_2017 = np.sum(flashes_2017_PL_d03, axis=0)
## Total flashes summer 2018
flashes_2018_PL_d02 = CC_CG_PL_d02[418:505,:,:]
flashes_2018_PL_d03 = CC_CG_PL_d03[418:505,:,:]
CC_CG_PL_d02_2018 = np.sum(flashes_2018_PL_d02, axis=0)
CC_CG_PL_d03_2018 = np.sum(flashes_2018_PL_d03, axis=0)
## Total flashes summer 2019
flashes_2019_PL_d02 = CC_CG_PL_d02[540:624,:,:]
flashes_2019_PL_d03 = CC_CG_PL_d03[540:624,:,:]
CC_CG_PL_d02_2019 = np.sum(flashes_2019_PL_d02, axis=0)
CC_CG_PL_d03_2019 = np.sum(flashes_2019_PL_d03, axis=0)
## Total flashes summer 2020
flashes_2020_PL_d02 = CC_CG_PL_d02[660:746,:,:]
flashes_2020_PL_d03 = CC_CG_PL_d03[660:746,:,:]
CC_CG_PL_d02_2020 = np.sum(flashes_2020_PL_d02, axis=0)
CC_CG_PL_d03_2020 = np.sum(flashes_2020_PL_d03, axis=0)

## Total flashes summer 2015
flashes_2015_d02 = CC_CG_d02[25:112]
flashes_2015_d03 = CC_CG_d03[25:112]
## Total flashes summer 2016
flashes_2016_d02 = CC_CG_d02[159:248]
flashes_2016_d03 = CC_CG_d03[159:248]
## Total flashes summer 2017
flashes_2017_d02 = CC_CG_d02[289:376]
flashes_2017_d03 = CC_CG_d03[289:376]
## Total flashes summer 2018
flashes_2018_d02 = CC_CG_d02[418:505]
flashes_2018_d03 = CC_CG_d03[418:505]
## Total flashes summer 2019
flashes_2019_d02 = CC_CG_d02[540:624]
flashes_2019_d03 = CC_CG_d03[540:624]
## Total flashes summer 2020
flashes_2020_d02 = CC_CG_d02[660:746]
flashes_2020_d03 = CC_CG_d03[660:746]

# Calculate total flashes per year (sum / 4)
Total_flashes_4_year_d02 = np.add(CC_CG_PL_d02_2015, CC_CG_PL_d02_2016)
Total_flashes_4_year_d02 = np.add(Total_flashes_4_year_d02, CC_CG_PL_d02_2017)
Total_flashes_4_year_d02 = np.add(Total_flashes_4_year_d02, CC_CG_PL_d02_2018)

Total_flashes_4_year_d03 = np.add(CC_CG_PL_d03_2015, CC_CG_PL_d03_2016)
Total_flashes_4_year_d03 = np.add(Total_flashes_4_year_d03, CC_CG_PL_d03_2017)
Total_flashes_4_year_d03 = np.add(Total_flashes_4_year_d03, CC_CG_PL_d03_2018)

TF_d02 = np.divide(Total_flashes_4_year_d02, 4)
TF_d03 = np.divide(Total_flashes_4_year_d03, 4)

## Only get the summer months
time_obs1 = time_obs.year
time_obs2 = np.where([(time_obs1 < 2019)])
time_obs3 = time_obs[time_obs2[1]]
months = time_obs3.month
summermonths = np.where([(months > 5) & (months < 9)])
lats_d02_obs = lats_d02
lons_d02_obs = lons_d02
lats_d03_obs = lats_d03
lons_d03_obs = lons_d03
# ---------------------------------------------------------------------------------------------
# MODEL
# ---------------------------------------------------------------------------------------------
f2015_d02 = '/scratch/leuven/336/vsc33651/nu-wrf-dev/WRFPRS_d02_2015.nc'
ds2015_d02 = Dataset(f2015_d02, 'r')
lats_2015_d02 = ds2015_d02.variables['lat'][:]
lons_2015_d02 = ds2015_d02.variables['lon'][:]
time_2015_d02 = ds2015_d02.variables['time'][:]
Cdim_2015_d02 = ds2015_d02.variables['Cdim'][:]
CAPE_2015_d02 = ds2015_d02.variables['CAPE'][:]
F_PH_2015_d02 = ds2015_d02.variables['F_PH'][:]

f2016_d02 = '/scratch/leuven/336/vsc33651/nu-wrf-dev/WRFPRS_d02_2016.nc'
ds2016_d02 = Dataset(f2016_d02, 'r')
lats_2016_d02 = ds2016_d02.variables['lat'][:]
lons_2016_d02 = ds2016_d02.variables['lon'][:]
time_2016_d02 = ds2016_d02.variables['time'][:]
Cdim_2016_d02 = ds2016_d02.variables['Cdim'][:]
CAPE_2016_d02 = ds2016_d02.variables['CAPE'][:]
F_PH_2016_d02 = ds2016_d02.variables['F_PH'][:]

f2017_d02 = '/scratch/leuven/336/vsc33651/nu-wrf-dev/WRFPRS_d02_2017.nc'
ds2017_d02 = Dataset(f2017_d02, 'r')
lats_2017_d02 = ds2017_d02.variables['lat'][:]
lons_2017_d02 = ds2017_d02.variables['lon'][:]
time_2017_d02 = ds2017_d02.variables['time'][:]
Cdim_2017_d02 = ds2017_d02.variables['Cdim'][:]
CAPE_2017_d02 = ds2017_d02.variables['CAPE'][:]
F_PH_2017_d02 = ds2017_d02.variables['F_PH'][:]

f2018_d02 = '/scratch/leuven/336/vsc33651/nu-wrf-dev/WRFPRS_d02_2018.nc'
ds2018_d02 = Dataset(f2018_d02, 'r')
lats_2018_d02 = ds2018_d02.variables['lat'][:]
lons_2018_d02 = ds2018_d02.variables['lon'][:]
time_2018_d02 = ds2018_d02.variables['time'][:]
Cdim_2018_d02 = ds2018_d02.variables['Cdim'][:]
CAPE_2018_d02 = ds2018_d02.variables['CAPE'][:]
F_PH_2018_d02 = ds2018_d02.variables['F_PH'][:]

f2015_d03 = '/scratch/leuven/336/vsc33651/nu-wrf-dev/WRFPRS_d03_2015.nc'
ds2015_d03 = Dataset(f2015_d03, 'r')
lats_2015_d03 = ds2015_d03.variables['lat'][:]
lons_2015_d03 = ds2015_d03.variables['lon'][:]
Cdim_2015_d03 = ds2015_d03.variables['Cdim'][:]
CAPE_2015_d03 = ds2015_d03.variables['CAPE'][:]
F_PH_2015_d03 = ds2015_d03.variables['F_PH'][:]

f2016_d03 = '/scratch/leuven/336/vsc33651/nu-wrf-dev/WRFPRS_d03_2016.nc'
ds2016_d03 = Dataset(f2016_d03, 'r')
lats_2016_d03 = ds2016_d03.variables['lat'][:]
lons_2016_d03 = ds2016_d03.variables['lon'][:]
Cdim_2016_d03 = ds2016_d03.variables['Cdim'][:]
CAPE_2016_d03 = ds2016_d03.variables['CAPE'][:]
F_PH_2016_d03 = ds2016_d03.variables['F_PH'][:]

f2017_d03 = '/scratch/leuven/336/vsc33651/nu-wrf-dev/WRFPRS_d03_2017.nc'
ds2017_d03= Dataset(f2017_d03, 'r')
lats_2017_d03 = ds2017_d03.variables['lat'][:]
lons_2017_d03 = ds2017_d03.variables['lon'][:]
Cdim_2017_d03 = ds2017_d03.variables['Cdim'][:]
CAPE_2017_d03 = ds2017_d03.variables['CAPE'][:]
F_PH_2017_d03 = ds2017_d03.variables['F_PH'][:]

f2018_d03 = '/scratch/leuven/336/vsc33651/nu-wrf-dev/WRFPRS_d03_2018.nc'
ds2018_d03= Dataset(f2018_d03, 'r')
lats_2018_d03 = ds2018_d03.variables['lat'][:]
lons_2018_d03 = ds2018_d03.variables['lon'][:]
time_2018_d03 = ds2018_d03.variables['time'][:]
Cdim_2018_d03 = ds2018_d03.variables['Cdim'][:]
CAPE_2018_d03 = ds2018_d03.variables['CAPE'][:]
F_PH_2018_d03 = ds2018_d03.variables['F_PH'][:]

# Add all years in one array (per domain)
lats_d02 = lats_2015_d02 #lats are for all years the same
lons_d02 = lons_2015_d02 #lons as well
time_4_yr = np.zeros((8832))
time_4_yr[0:2208] = time_2015_d02
time_4_yr[2208:4416] = time_2016_d02
time_4_yr[4416:6624] = time_2017_d02
time_4_yr[6624:8832] = time_2018_d02
F_PH_d02 = np.zeros((8832,216,363))
F_PH_d02[0:2208,:,:] = F_PH_2015_d02
F_PH_d02[2208:4416,:,:] = F_PH_2016_d02
F_PH_d02[4416:6624,:,:] = F_PH_2017_d02
F_PH_d02[6624:8832,:,:] = F_PH_2018_d02
CAPE_d02 = np.zeros((8832,216,363))
CAPE_d02[0:2208,:,:] = CAPE_2015_d02
CAPE_d02[2208:4416,:,:] = CAPE_2016_d02
CAPE_d02[4416:6624,:,:] = CAPE_2017_d02
CAPE_d02[6624:8832,:,:] = CAPE_2018_d02
Cdim_d02 = np.zeros((8832,216,363))
Cdim_d02[0:2208,:,:] = Cdim_2015_d02
Cdim_d02[2208:4416,:,:] = Cdim_2016_d02
Cdim_d02[4416:6624,:,:] = Cdim_2017_d02
Cdim_d02[6624:8832,:,:] = Cdim_2018_d02

lats_d03 = lats_2015_d03 #lats are for all years the same
lons_d03 = lons_2015_d03 #lons as well
F_PH_d03 = np.zeros((8832,111,114))
F_PH_d03[0:2208,:,:] = F_PH_2015_d03
F_PH_d03[2208:4416,:,:] = F_PH_2016_d03
F_PH_d03[4416:6624,:,:] = F_PH_2017_d03
F_PH_d03[6624:8832,:,:] = F_PH_2018_d03
CAPE_d03 = np.zeros((8832,111,114))
CAPE_d03[0:2208,:,:] = CAPE_2015_d03
CAPE_d03[2208:4416,:,:] = CAPE_2016_d03
CAPE_d03[4416:6624,:,:] = CAPE_2017_d03
CAPE_d03[6624:8832,:,:] = CAPE_2018_d03
Cdim_d03 = np.zeros((8832,111,114))
Cdim_d03[0:2208,:,:] = Cdim_2015_d03
Cdim_d03[2208:4416,:,:] = Cdim_2016_d03
Cdim_d03[4416:6624,:,:] = Cdim_2017_d03
Cdim_d03[6624:8832,:,:] = Cdim_2018_d03

# ---------------------------------------------------------------------------------------------
# DENSITY PLOTS
# ---------------------------------------------------------------------------------------------

# "Viridis-like" colormap with white background
white_viridis = LinearSegmentedColormap.from_list('white_viridis', [
    (0, '#ffffff'),
    (1e-20, '#440053'),
    (0.2, '#404388'),
    (0.4, '#2a788e'),
    (0.6, '#21a784'),
    (0.8, '#78d151'),
    (1, '#fde624'),
], N=256)

def using_mpl_scatter_density(fig, x, y):
    ax = fig.add_subplot(1, 1, 1, projection='scatter_density')
    density = ax.scatter_density(x, y, cmap=white_viridis)
    fig.colorbar(density, label='Number of points per pixel')

def density_plots(x, y, axis=[], title='', xlabel='', ylabel=''):
    x1 = x[np.logical_not(np.isnan(x))]
    y1 = y[np.logical_not(np.isnan(x))]
    x2 = x1[np.logical_not(np.isnan(y1))]
    y2 = y1[np.logical_not(np.isnan(y1))]
    x3 = x2[np.logical_not(np.isinf(x2))]
    y3 = y2[np.logical_not(np.isinf(x2))]
    x4 = x3[np.logical_not(np.isinf(y3))]
    y4 = y3[np.logical_not(np.isinf(y3))]
    x = x4
    y = y4
    fig = plt.figure()
    using_mpl_scatter_density(fig, x, y)
    if axis != []:
        plt.axis(axis)
    if title != '':
        plt.title(title)
    if xlabel != '':
        plt.xlabel(xlabel)
    if ylabel != '':
        plt.ylabel(ylabel)
    plt.show()

# Get daily averages
# d02
days = np.divide(8832, 24).astype(int)
daily_avg_F_PH_d02 = np.zeros((days,216,363))
F_PH_d02_tot = np.sum(F_PH_d02, axis=(1,2))
F_PH_d02_tot_pl = np.zeros((days,216,363))
daily_tot_F_PH_d02 = np.zeros((days))
daily_avg_CAPE_d02 = np.zeros((days,216,363))
daily_avg_Cdim_d02 = np.zeros((days,216,363))
time_4_yr_days = np.zeros((days))
for i in range(0,days):
    j = (i+1)*24
    k = i*24
    daily_avg_F_PH_d02[i,:,:] = np.mean(F_PH_d02[k:j,:,:], axis = 0)
    daily_tot_F_PH_d02[i] = np.sum(F_PH_d02_tot[k:j])
    daily_avg_CAPE_d02[i,:,:] = np.mean(CAPE_d02[k:j,:,:], axis=0)
    daily_avg_Cdim_d02[i,:,:] = np.mean(Cdim_d02[k:j,:,:], axis=0)
    F_PH_d02_tot_pl[i,:,:] = np.sum(F_PH_d02[k:j,:,:], axis = 0)
    time_4_yr_days[i] = time_4_yr[k]

max_f_pd_d02 = np.amax(daily_tot_F_PH_d02)
max_f_pl_d02 = np.amax(F_PH_d02_tot_pl, axis=(1,2))
stand_f_pd_d02 = np.divide(daily_tot_F_PH_d02, max_f_pd_d02)
stand_f_pl_d02 = np.zeros((368,216,363))
avg_CAPE_pl_d02 = np.mean(CAPE_d02, axis=0)
for i in range(0,368):
    stand_f_pl_d02[i,:,:] = np.divide(F_PH_d02_tot_pl[i,:,:], max_f_pl_d02[i])
mean_threat_d02 = np.mean(stand_f_pl_d02, axis = 0)
# d03
daily_avg_F_PH_d03 = np.zeros((days,111,114))
F_PH_d03_tot = np.sum(F_PH_d03, axis=(1,2))
F_PH_d03_tot_pl = np.zeros((days,111,114))
daily_tot_F_PH_d03 = np.zeros((days))
daily_avg_CAPE_d03 = np.zeros((days,111,114))
daily_avg_Cdim_d03 = np.zeros((days,111,114))
for i in range(0,days):
    j = (i+1)*24
    k = i*24
    daily_avg_F_PH_d03[i,:,:] = np.mean(F_PH_d03[k:j,:,:], axis = 0)
    daily_tot_F_PH_d03[i] = np.sum(F_PH_d03_tot[k:j])
    daily_avg_CAPE_d03[i,:,:] = np.mean(CAPE_d03[k:j,:,:], axis=0)
    daily_avg_Cdim_d03[i,:,:] = np.mean(Cdim_d03[k:j,:,:], axis=0)
    F_PH_d03_tot_pl[i,:,:] = np.sum(F_PH_d03[k:j,:,:], axis = 0)

max_f_pd_d03 = np.amax(daily_tot_F_PH_d03)
max_f_pl_d03 = np.amax(F_PH_d03_tot_pl, axis=(1,2))
stand_f_pd_d03 = np.divide(daily_tot_F_PH_d03, max_f_pd_d03)
stand_f_pl_d03 = np.zeros((368,111,114))

for i in range(0,368):
    stand_f_pl_d03[i,:,:] = np.divide(F_PH_d03_tot_pl[i,:,:], max_f_pl_d03[i])
mean_threat_d03 = np.nanmean(stand_f_pl_d03, axis = 0)
## Plot timeseries
# density_plots(daily_avg_F_PH_d02, daily_avg_CAPE_d02, title='Daily average CAPE versus lightning flashes per hour (domain 2)'
#               , xlabel = 'Total number of flashes [hour$^{-1}$]', ylabel = 'CAPE [Jkg$^{-1}$]')
# density_plots(daily_avg_F_PH_d02, daily_avg_Cdim_d02, title='Daily average cloud dimensions versus lightning flashes (domain 2)'
#               , xlabel = 'Total number of flashes [hour$^{-1}$]', ylabel = 'Cloud dimensions [km]')
# density_plots(daily_avg_CAPE_d02, daily_avg_Cdim_d02, title='Daily average cloud dimensions versus CAPE (domain 2)'
#               , xlabel = 'CAPE [Jkg$^{-1}$]', ylabel = 'Cloud dimensions [km]')

# Plot d02_CC + d02_CG together (so total in d02) (for 4 years, only in summer)
d_4_yr = pd.date_range(start='2015-06-01',end='2018-09-1',freq='d')
months= d_4_yr.month
summer = np.where([(months > 5) & (months < 9)])
daily_CAPE_d02 = np.sum(daily_avg_CAPE_d02, axis=(1,2))
# # 4 year time series CAPE domain 2 and observed lightning flashes
plt.figure(figsize=(20,20))
plt.subplot(2,1,1)
plt.plot(time_obs3[summermonths[1]], CC_CG_d02[summermonths[1]],'.')
plt.grid()
plt.ylabel("Total number of flashes [day$^{-1}$]")
plt.title("Domain 2")
plt.tick_params(labelbottom=False)
plt.subplot(2,1,2)
plt.plot(d_4_yr[summer[1]], stand_f_pd_d02,'.')
plt.grid()
plt.ylabel("Lightning threat [-]")
plt.show()

#
# plt.figure(figsize=(20,20))
# plt.plot(d_4_yr[summer[1]], stand_f_pd_d03,'.')
# plt.grid()
# plt.ylabel("lightning threat")
# plt.title("Domain 3")
# plt.show()

switch_interactive = 0
## Map Flashdensity in domain
f = plt.figure(num=None, dpi=600, facecolor='w', edgecolor='k')
m = Basemap(projection= 'lcc', lat_0 = lats_d03.mean(),
            lon_0 = lons_d03.mean(),
            llcrnrlon=(np.amin(lons_d03)), llcrnrlat=(np.amin(lats_d03)),
            urcrnrlon=(np.amax(lons_d03)),urcrnrlat = (np.amax(lats_d03)),
            resolution = 'f')
m.drawcoastlines(linewidth=0.3)
m.drawcountries()

m.drawparallels(np.arange(-50,70,5))
m.drawmeridians(np.arange(-130,-100,5))

seismic_mod = cm.get_cmap('seismic',256)
newcmp = ListedColormap(seismic_mod(np.linspace(0.5,1,256)))
#x,y = np.meshgrid(lons[lon_inds_d02],lats[lat_inds_d02])
m.pcolormesh(lons_2015_d03,lats_2015_d03,mean_threat_d03,latlon=True,cmap=newcmp)

if switch_interactive=="1":
    plt.show()
else:
    outpath = "/scratch/leuven/336/vsc33651/nu-wrf-dev/"
    fname = "threat_d03"
    fname_long = os.path.join(outpath, fname+'.png')
    plt.title('Average Lightning Threat (June - August')
    plt.colorbar()
    plt.savefig(fname_long, dpi=f.dpi)
    plt.close()