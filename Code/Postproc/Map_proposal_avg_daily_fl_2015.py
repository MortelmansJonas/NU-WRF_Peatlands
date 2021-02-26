## Load the necessary modules
import numpy as np
import pandas as pd
import netCDF4 as nc
from matplotlib.colors import LogNorm
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
import scipy.stats as st

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
WRF_FILES = ["wrfout_d02_2015-06-01_00:00:00"]

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
print('All tests passed!')

file_path_d02 = d02_wrf_file()
d02_wrf_file = Dataset(file_path_d02)

## average daily Flashdensity for d02
# D02: Extract latitudes and longitudes from WRF output.
terrain_d02 = getvar(d02_wrf_file, "ter", timeidx=0)
lats_d02, lons_d02 = latlon_coords(terrain_d02)
print(terrain_d02)
min_lat_d02 = np.array(lats_d02.min(), dtype='f4')
max_lat_d02 = np.array(lats_d02.max(), dtype='f4')
min_lon_d02 = np.array(lons_d02.min(), dtype='f4')
max_lon_d02 = np.array(lons_d02.max(), dtype='f4')
lat_inds_d02 = np.where((lats >= min_lat_d02) & (lats <= max_lat_d02))[0]
lon_inds_d02 = np.where((lons >= min_lon_d02) & (lons <= max_lon_d02))[0]
FD_CC_D02 = ds.variables['Flashdensity_CC'][:,lat_inds_d02, lon_inds_d02]
FD_CG_D02 = ds.variables['Flashdensity_CG'][:,lat_inds_d02, lon_inds_d02]

## Calculate CC+CG
CC_CG_d02 = np.add(FD_CC_D02, FD_CG_D02)
## Total flashes summer 2015
flashes_2015_PL_d02 = CC_CG_d02[25:112,:,:]
# CC_CG_PL_d02_2015 = np.sum(flashes_2015_PL_d02, axis=0)

## Only get the summer months
time_obs1 = time_obs.year
time_obs2 = np.where([(time_obs1 == 2015)])
time_obs3 = time_obs[time_obs2[1]]
months = time_obs3.month
summermonths = np.where([(months > 5) & (months < 9)])
lats_d02_obs = lats_d02
lons_d02_obs = lons_d02


# Calculate average daily flashes 2015
# Avg_obs_flashes_pd_2015 = np.divide(CC_CG_PL_d02_2015,len(flashes_2015_PL_d02[0]))
Avg_obs_flashes_pd_2015 = np.mean(flashes_2015_PL_d02,axis=0)
# ---------------------------------------------------------------------------------------------
# MODEL
# ---------------------------------------------------------------------------------------------
#d02
f2015_d02 = '/scratch/leuven/336/vsc33651/nu-wrf-dev/wrfout_nc_files/wrfout_d02_2015.nc'
ds2015_d02 = Dataset(f2015_d02, 'r')
lats_2015_d02 = ds2015_d02.variables['lat'][:]
lons_2015_d02 = ds2015_d02.variables['lon'][:]
time_2015_d02 = ds2015_d02.variables['time'][:]
LTG1_MAX_2015_d02 = ds2015_d02.variables['LTG1_MAX'][:]
LTG2_MAX_2015_d02 = ds2015_d02.variables['LTG2_MAX'][:]
LTG3_MAX_2015_d02 = ds2015_d02.variables['LTG3_MAX'][:]
CTOP2D_2015_d02 = ds2015_d02.variables['CTOP2D'][:]
COD2D_2015_d02 = ds2015_d02.variables['COD2D'][:]

#d01
f2015_d01 = '/scratch/leuven/336/vsc33651/nu-wrf-dev/wrfout_nc_files/wrfout_d01_2015.nc'
ds2015_d01 = Dataset(f2015_d01, 'r')
lats_2015_d01 = ds2015_d01.variables['lat'][:]
lons_2015_d01 = ds2015_d01.variables['lon'][:]
maxlatd02 = np.amax(lats_2015_d02)
minlatd02 = np.amin(lats_2015_d02)
maxlond02 = np.amax(lons_2015_d02)
minlond02 = np.amin(lons_2015_d02)
lat_inds_d02_d01 = np.where((lats_2015_d01 > minlatd02) & (lats_2015_d01 < maxlatd02))
lon_inds_d02_d01 = np.where((lons_2015_d01 > minlond02) & (lons_2015_d01 < maxlond02))
latinds_d02 = np.unique(lat_inds_d02_d01[0])
loninds_d02 = np.unique(lon_inds_d02_d01[1])
LTG1_MAX_2015_d01 = ds2015_d01.variables['LTG1_MAX'][:, latinds_d02, loninds_d02]
LTG2_MAX_2015_d01 = ds2015_d01.variables['LTG2_MAX'][:, latinds_d02, loninds_d02]
LTG3_MAX_2015_d01 = ds2015_d01.variables['LTG3_MAX'][:, latinds_d02, loninds_d02]
CTOP2D_2015_d01 = ds2015_d01.variables['CTOP2D'][:, latinds_d02, loninds_d02]
COD2D_2015_d01 = ds2015_d01.variables['COD2D'][:, latinds_d02, loninds_d02]

# ---------------------------------------------------------------------------------------------
# CALCULATING DAILY AVERAGE THREAT
# ---------------------------------------------------------------------------------------------
#d01
days = np.divide(2208,24).astype(int)
daily_LTG1_MAX_d01 = np.zeros((days,90,157))
daily_LTG2_MAX_d01 = np.zeros((days,90,157))
daily_LTG3_MAX_d01 = np.zeros((days,90,157))
for i in range(0,days):
    j = (i+1)*24
    k = i*24
    daily_LTG1_MAX_d01[i,:,:] = np.mean(LTG1_MAX_2015_d01[k:j,:,:], axis = 0)
    daily_LTG2_MAX_d01[i, :, :] = np.mean(LTG2_MAX_2015_d01[k:j, :, :], axis=0)
    daily_LTG3_MAX_d01[i, :, :] = np.mean(LTG3_MAX_2015_d01[k:j, :, :], axis=0)

#d02
daily_LTG1_MAX_d02 = np.zeros((days,216,363))
daily_LTG2_MAX_d02 = np.zeros((days,216,363))
daily_LTG3_MAX_d02 = np.zeros((days,216,363))
for i in range(0,days):
    j = (i+1)*24
    k = i*24
    daily_LTG1_MAX_d02[i,:,:] = np.mean(LTG1_MAX_2015_d02[k:j,:,:], axis = 0)
    daily_LTG2_MAX_d02[i, :, :] = np.mean(LTG2_MAX_2015_d02[k:j, :, :], axis=0)
    daily_LTG3_MAX_d02[i, :, :] = np.mean(LTG3_MAX_2015_d02[k:j, :, :], axis=0)

# Calculate daily average of all summer
avg_daily_threat3_d01 = np.mean(daily_LTG3_MAX_d01, axis = 0)
avg_daily_threat3_d02 = np.mean(daily_LTG3_MAX_d02, axis = 0)

# Derive lightning flashes from threats
obs_lightning_mean = np.mean(Avg_obs_flashes_pd_2015)
mod_threat_mean_d01 = np.mean(avg_daily_threat3_d01)
mod_threat_mean_d02 = np.mean(avg_daily_threat3_d02)
mod_flashes_d01 = np.divide(np.multiply(obs_lightning_mean, avg_daily_threat3_d01),mod_threat_mean_d01)
mod_flashes_d02 = np.divide(np.multiply(obs_lightning_mean, avg_daily_threat3_d02),mod_threat_mean_d02)
print(np.nanmin(mod_flashes_d01))
# ---------------------------------------------------------------------------------------------
# MAPS
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

switch_interactive = 0
## Map Observed Flashdensity in domain
f = plt.figure(num=None, dpi=700, facecolor='w', edgecolor='k')
m = Basemap(projection= 'lcc', lat_0 = lats_d02.mean(),
            lon_0 = lons_d02.mean(),
            llcrnrlon=(-125.5), llcrnrlat=(57.2),
            urcrnrlon=(-105.4),urcrnrlat = (63),
            resolution = 'f')
m.drawcoastlines(linewidth=0.3)
m.drawcountries()
m.drawparallels(np.arange(-50,70,5))
m.drawmeridians(np.arange(-130,-100,5))
m.drawmapscale(np.amin(lons_2015_d02)+2.5, np.amin(lats_2015_d02)+0.6, lons_d02.mean(), lats_d02.mean(), 100)
m.drawparallels(np.arange(-50,70,5), linewidth=0.5,labels = [True,False,False,False])
m.drawmeridians(np.arange(-130,-100,5), linewidth=0.5,labels = [False,False,False, True])

Avg_obs_flashes_pd_2015[Avg_obs_flashes_pd_2015<=0]=10**(-1)

seismic_mod = cm.get_cmap('seismic',256)
newcmp = ListedColormap(seismic_mod(np.linspace(0.5,1,256)))
x,y = np.meshgrid(lons[lon_inds_d02],lats[lat_inds_d02])
m.pcolormesh(x,y,Avg_obs_flashes_pd_2015,latlon=True,cmap=newcmp)

if switch_interactive=="1":
    plt.show()
else:
    outpath = "/scratch/leuven/336/vsc33651/nu-wrf-dev/Plots"
    fname = "avg_lf_pd_obs_2015_d02"
    fname_long = os.path.join(outpath, fname+'.png')
    plt.title('Observations',size=12)
    plt.colorbar(label='Number of flashes (day$^{-1}$ km$^{-2}$)',shrink = 0.6)
    # plt.clim(10**(-1),10**2)
    plt.savefig(fname_long, dpi=f.dpi)
    plt.close()

## Map LTG threat 3 in domain 1
f = plt.figure(num=None, dpi=700, facecolor='w', edgecolor='k')
m = Basemap(projection= 'lcc', lat_0 = lats_d02.mean(),
            lon_0 = lons_d02.mean(),
            llcrnrlon=(-125.5), llcrnrlat=(57.2),
            urcrnrlon=(-105.4),urcrnrlat = (63),
            resolution = 'f')

m.drawcoastlines(linewidth=0.3)
m.drawcountries()
m.drawparallels(np.arange(-50,70,5))
m.drawmeridians(np.arange(-130,-100,5))
m.drawparallels(np.arange(-50,70,5), linewidth=0.5,labels = [True,False,False,False])
m.drawmeridians(np.arange(-130,-100,5), linewidth=0.5,labels = [False,False,False, True])

seismic_mod = cm.get_cmap('seismic',256)
latitudes = np.arange(minlatd02, maxlatd02, (maxlatd02 - minlatd02)/90)
longitudes = np.arange(minlond02, maxlond02, (maxlond02 - minlond02)/157)
newcmp = ListedColormap(seismic_mod(np.linspace(0.5,1,256)))
x,y = np.meshgrid(longitudes,latitudes)
m.pcolormesh(x,y, mod_flashes_d01,latlon=True,cmap=newcmp)

if switch_interactive=="1":
    plt.show()
else:
    outpath = "/scratch/leuven/336/vsc33651/nu-wrf-dev/Plots"
    fname = "avg_lf_pd_2015_d02_d01"
    fname_long = os.path.join(outpath, fname+'.png')
    plt.title('Convection parameterized (9-km)', size=12)
    plt.colorbar(label='Number of flashes (day$^{-1}$ km$^{-2}$)', extend='max',shrink = 0.6)
    # plt.clim(10**(-1),10**2)
    # plt.clim(0,20)
    plt.savefig(fname_long, dpi=f.dpi)
    plt.close()

## Domain 2
f = plt.figure(num=None, dpi=700, facecolor='w', edgecolor='k')
m = Basemap(projection= 'lcc', lat_0 = lats_d02.mean(),
            lon_0 = lons_d02.mean(),
            llcrnrlon=(-125.5), llcrnrlat=(57.2),
            urcrnrlon=(-105.4),urcrnrlat = (63),
            resolution = 'f')
m.drawcoastlines(linewidth=0.3)
m.drawcountries()
m.drawparallels(np.arange(-50,70,5), linewidth=0.5,labels = [True,False,False,False])
m.drawmeridians(np.arange(-130,-100,5), linewidth=0.5,labels = [False,False,False, True])

seismic_mod = cm.get_cmap('seismic',256)
newcmp = ListedColormap(seismic_mod(np.linspace(0.5,1,256)))
m.pcolormesh(lons_2015_d02,lats_2015_d02, mod_flashes_d02,latlon=True,cmap=newcmp)

if switch_interactive=="1":
    plt.show()
else:
    outpath = "/scratch/leuven/336/vsc33651/nu-wrf-dev/Plots"
    fname = "avg_lf_pd_2015_d02"
    fname_long = os.path.join(outpath, fname+'.png')
    plt.title('Convection-permitting (3-km)', size = 12)
    plt.colorbar(label='Number of flashes (day$^{-1}$ km$^{-2}$)', extend='max', shrink = 0.6)
  #  plt.clim(10**(-1),10**2)
   # plt.clim(0,4)
    plt.savefig(fname_long, dpi=f.dpi)
    plt.close()

# ---------------------------------------------------------------------------------------------
# PLOTS DIURNAL CYCLE
# ---------------------------------------------------------------------------------------------

# CALCULATE AVERAGE FOR EACH HOUR OF THE DAY
## OBSERVATIONS
# nf= '/scratch/leuven/336/vsc33651/nu-wrf-dev/Observations_diurnal.nc'
# ds_hourly = Dataset(nf, 'r')
# hour_of_day = ds_hourly['time'][:]
# LF_diurnal = ds_hourly['LF_diurnal'][:]
# plt.plot(hour_of_day, np.nanmean(LF_diurnal,axis=(1,2)), 'k')
# plt.xlabel('Hour of the day', size = 12)
# plt.ylabel('Number of flashes (day$^{-1}$ km$^{-2}$)', size = 12)
# plt.ticklabel_format(axis='y',style='sci', scilimits=(0,0))
# plt.show()
#
# ## DOMAIN 1
# ds_diurnal_d01 = Dataset('/scratch/leuven/336/vsc33651/nu-wrf-dev/wrfout_nc_files/2015_d01_diurnal.nc', mode='r')
# diurnal_LTG3_d01 = ds_diurnal_d01['LTG3_MAX'][:]
# hour_of_day = ds_diurnal_d01['time'][:]
# plt.plot(hour_of_day, np.mean(diurnal_LTG3_d01,axis=(1,2)), 'k')
# # plt.fill_between(hour_of_day, np.mean(diurnal_LTG3_d01,axis=(1,2))-np.std(diurnal_LTG3_d01, axis=(1,2)),
# #                  np.mean(diurnal_LTG3_d01,axis=(1,2))+np.std(diurnal_LTG3_d01, axis=(1,2)))
# plt.xlabel('Hour of the day', size = 12)
# plt.ylabel('Lightning threat (day$^{-1}$ (9km)$^{-2}$)', size = 12)
# plt.ticklabel_format(axis='y',style='sci', scilimits=(0,0))
# plt.show()
#
# # # DOMAIN 2
#
# ds_diurnal_d02 = Dataset('/scratch/leuven/336/vsc33651/nu-wrf-dev/wrfout_nc_files/2015_d02_diurnal.nc', mode='r')
# diurnal_LTG3_d02 = ds_diurnal_d02['LTG3_MAX'][:]
# plt.plot(hour_of_day, np.mean(diurnal_LTG3_d02,axis=(1,2)), 'k')
# #plt.fill_between(hour_of_day, np.mean(diurnal_LTG3_d02,axis=(1,2))-np.std(diurnal_LTG3_d02, axis=(1,2)),
# #                 np.mean(diurnal_LTG3_d02,axis=(1,2))+np.std(diurnal_LTG3_d02, axis=(1,2)))
# plt.xlabel('Hour of the day', size = 12)
# plt.ylabel('Lightning threat (day$^{-1}$ (9km)$^{-2}$)', size = 12)
# plt.ticklabel_format(axis='y',style='sci', scilimits=(0,0))
# plt.show()