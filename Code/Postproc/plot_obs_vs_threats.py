## Load the necessary modules
import numpy as np
import pandas as pd
import netCDF4 as nc
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.axis as axis
from matplotlib.cm import get_cmap
from matplotlib.colors import from_levels_and_colors, ListedColormap, LinearSegmentedColormap
from numpy import meshgrid
from matplotlib import cm
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
from datetime import date, timedelta
from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()
from netCDF4 import Dataset
from wrf import (getvar, latlon_coords)
import os

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
# THREATS
# ---------------------------------------------------------------------------------------------
f2015_d02 = '/scratch/leuven/336/vsc33651/nu-wrf-dev/wrfout_nc_files/wrfout_d02_2015_lightningtest.nc'
ds2015_d02 = Dataset(f2015_d02, 'r')
lats_2015_d02 = ds2015_d02.variables['lat'][:]
lons_2015_d02 = ds2015_d02.variables['lon'][:]
time_2015_d02 = ds2015_d02.variables['time'][:]
LTG1_MAX_2015_d02 = ds2015_d02.variables['LTG1_MAX'][:]
LTG2_MAX_2015_d02 = ds2015_d02.variables['LTG2_MAX'][:]
LTG3_MAX_2015_d02 = ds2015_d02.variables['LTG3_MAX'][:]

# ---------------------------------------------------------------------------------------------
# CALCULATING DAILY AVERAGE THREAT
# ---------------------------------------------------------------------------------------------
days = np.divide(216,24).astype(int)
daily_LTG1_MAX_d02 = np.zeros((days,216,363))
daily_LTG2_MAX_d02 = np.zeros((days,216,363))
daily_LTG3_MAX_d02 = np.zeros((days,216,363))
for i in range(0,days):
    j = (i+1)*24
    k = i*24
    daily_LTG1_MAX_d02[i,:,:] = np.mean(LTG1_MAX_2015_d02[k:j,:,:], axis = 0)
    daily_LTG2_MAX_d02[i, :, :] = np.mean(LTG2_MAX_2015_d02[k:j, :, :], axis=0)
    daily_LTG3_MAX_d02[i, :, :] = np.mean(LTG3_MAX_2015_d02[k:j, :, :], axis=0)

# ---------------------------------------------------------------------------------------------
# PLOTTING
# ---------------------------------------------------------------------------------------------
switch_interactive = 0
## Map Flashdensity in domain
f = plt.figure(num=None, dpi=600, facecolor='w', edgecolor='k')
m = Basemap(projection= 'lcc', lat_0 = lats_2015_d02.mean(),
            lon_0 = lons_2015_d02.mean(),
            llcrnrlon=(np.amin(lons_2015_d02)+1), llcrnrlat=(np.amin(lats_2015_d02)-0.5),
            urcrnrlon=(np.amax(lons_2015_d02)+0.5),urcrnrlat = (np.amax(lats_2015_d02)-0.5),
            resolution = 'f')
m.drawcoastlines(linewidth=0.3)
m.drawcountries()

m.drawparallels(np.arange(-50,70,5))
m.drawmeridians(np.arange(-130,-100,5))

seismic_mod = cm.get_cmap('seismic',256)
newcmp = ListedColormap(seismic_mod(np.linspace(0.5,1,256)))
#x,y = np.meshgrid(lons[lon_inds_d02],lats[lat_inds_d02])
m.pcolormesh(lons_2015_d02,lats_2015_d02,flashes_2015_PL_d02,latlon=True,cmap='viridis')

if switch_interactive=="1":
    plt.show()
else:
    outpath = "/scratch/leuven/336/vsc33651/nu-wrf-dev/Plots"
    fname = "observed_ltg_406"
    fname_long = os.path.join(outpath, fname+'.png')
    plt.title('observed LTG on 4-06-2015')
    plt.colorbar()
    plt.savefig(fname_long, dpi=f.dpi)
    plt.close()