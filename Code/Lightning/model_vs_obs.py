## Script to compare model results with the observations.
## The daily average of lightning flashes in d03 will be compared with the daily average of flashes in the same area
## of the observations. Later on, this will be done for the dry conditions (linked with precipitation; lightning flashes
## that happen without precipitation) as these will be the most important flashes for lightning induced fires.

## modules
import netCDF4 as nc
import datetime as dt
from datetime import timedelta
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
import pandas as pd
from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()
from netCDF4 import Dataset
from wrf import (getvar, interplevel, vertcross, vinterp, ALL_TIMES, latlon_coords, get_cartopy, to_np,
                cartopy_xlim, cartopy_ylim, CoordPair)
import os
from matplotlib.cm import get_cmap
from matplotlib.colors import from_levels_and_colors
import xarray as xr
from numpy import meshgrid
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
Flashdensity_CC = np.where(np.isnan(Flashdensity_CC), 0,Flashdensity_CC) # To change NaN to 0
Flashdensity_CG = np.where(np.isnan(Flashdensity_CG), 0,Flashdensity_CG) # To change NaN to 0
## Convert all times to datetimes
time = pd.to_datetime(time)+timedelta(minutes=1) ## One minute is added to get correct days.
## Without adding 1 minute, some data was 1 day too early (1 minute). This is due to the conversion from datetime to int
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
min_lat_d02 = np.array(lats_d02.min(), dtype='int')
max_lat_d02 = np.array(lats_d02.max(), dtype='int')
min_lon_d02 = np.array(lons_d02.min(), dtype='int')
max_lon_d02 = np.array(lons_d02.max(), dtype='int')
lat_inds_d02 = np.where((lats >= min_lat_d02) & (lats <= max_lat_d02))[0]
lon_inds_d02 = np.where((lons >= min_lon_d02) & (lons <= max_lon_d02))[0]
FD_CC_D02 = ds.variables['Flashdensity_CC'][:,lat_inds_d02, lon_inds_d02]
FD_CG_D02 = ds.variables['Flashdensity_CG'][:,lat_inds_d02, lon_inds_d02]
FD_CC_D02 = np.where(np.isnan(FD_CC_D02), 0, FD_CC_D02) # To change NaN to 0
FD_CG_D02 = np.where(np.isnan(FD_CG_D02), 0, FD_CG_D02) # To change NaN to 0

# D03: Extract latitudes and longitudes from WRF output.
terrain_d03 = getvar(d03_wrf_file, "ter", timeidx=0)
lats_d03, lons_d03 = latlon_coords(terrain_d03)
min_lat_d03 = np.array(lats_d03.min(), dtype='int')
max_lat_d03 = np.array(lats_d03.max(), dtype='int')
min_lon_d03 = np.array(lons_d03.min(), dtype='int')
max_lon_d03 = np.array(lons_d03.max(), dtype='int')
lat_inds_d03 = np.where((lats >= min_lat_d03) & (lats <= max_lat_d03))[0]
lon_inds_d03 = np.where((lons >= min_lon_d03) & (lons <= max_lon_d03))[0]
FD_CC_D03 = ds.variables['Flashdensity_CC'][:,lat_inds_d03, lon_inds_d03]
FD_CG_D03 = ds.variables['Flashdensity_CG'][:,lat_inds_d03, lon_inds_d03]
FD_CC_D03 = np.where(np.isnan(FD_CC_D03), 0, FD_CC_D03) # To change NaN to 0
FD_CG_D03 = np.where(np.isnan(FD_CG_D03), 0, FD_CG_D03) # To change NaN to 0


## Group the total amount of flashes for each domain per day
FD_CC_D02_pd = FD_CC_D02.sum(axis=(1,2))
FD_CG_D02_pd = FD_CG_D02.sum(axis=(1,2))
FD_CC_D03_pd = FD_CC_D03.sum(axis=(1,2))
FD_CG_D03_pd = FD_CG_D03.sum(axis=(1,2))


## Calculate CG/(CC+CG)
CC_CG_d02 = np.add(FD_CC_D02_pd, FD_CG_D02_pd) # First calculate CC+CG for each domain
CC_CG_d02 = np.where(np.isnan(CC_CG_d02), 0, CC_CG_d02)
CC_CG_d03 = np.add(FD_CC_D03_pd, FD_CG_D03_pd) # First calculate CC+CG for each domain
CC_CG_d03 = np.where(np.isnan(CC_CG_d03), 0, CC_CG_d03)
CG_ratio_d02 = np.divide(FD_CG_D02_pd,CC_CG_d02)  # Then divide elementwise to get CG/(CC+CG)
CG_ratio_d02 = np.where(np.isnan(CG_ratio_d02), 0, CG_ratio_d02)
CG_ratio_d03 = np.divide(FD_CG_D03_pd,CC_CG_d03)  # Then divide elementwise to get CG/(CC+CG)
CG_ratio_d03 = np.where(np.isnan(CG_ratio_d03), 0, CG_ratio_d03)

##### MODEL (UPP threats)

# ##### PLOTS
# ## d02
# # Plot d02_CC
# plt.plot(time, FD_CC_D02_pd, label='CC')
# plt.grid()
# plt.ylabel("Number of flashes")
# plt.title('Number of CC flashes per day in d02')
# plt.show()
# # Plot d02_CG
# plt.plot(time, FD_CG_D02_pd, label='CG')
# plt.grid()
# plt.ylabel("Number of flashes")
# plt.title('Number of CG flashes per day in d02')
# plt.show()
# # Plot d02 CG/(CC+CG)
# plt.plot(time, CG_ratio_d02,'.')
# plt.grid()
# plt.ylabel("CG/(CC+CG)")
# plt.show()
#
# ## d03
# # Plot d03_CC
# plt.plot(time, FD_CC_D03_pd, label='CC')
# plt.grid()
# plt.ylabel("Number of flashes")
# plt.title('Number of CC flashes per day in d03')
# plt.show()
# # Plot d03_CG
# plt.plot(time, FD_CG_D03_pd, label='CG')
# plt.grid()
# plt.ylabel("number of flashes")
# plt.title('Number of CG flashes per day in d03')
# plt.show()
# # Plot d03 CG/(CC+CG)
# plt.plot(time, CG_ratio_d03,'.')
# plt.grid()
# plt.ylabel("CG/(CC+CG)")
# plt.show()

## Map Flashdensity in domain
m = Basemap(projection= 'lcc', lat_0 = ((min_lat_d02 + max_lat_d02)/2), lon_0 = ((min_lon_d02 + max_lon_d02)/2),
            llcrnrlon=min_lon_d02, llcrnrlat=min_lat_d02, urcrnrlon=max_lon_d02,urcrnrlat = max_lat_d02,
            resolution = 'h')
# m.fillcontinents(color='coral', lake_color='aqua')
# m.drawmapboundary(fill_color='aqua')
m.drawparallels(np.arange(-50,70,5))
m.drawmeridians(np.arange(-127,-100,5))
x,y = np.meshgrid(lons[lon_inds_d02],lats[lat_inds_d02])
m.pcolormesh(x,y,FD_CC_D02[1,:,:])

plt.show() #doesn't show any results