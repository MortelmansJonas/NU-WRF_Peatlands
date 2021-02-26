# ---------------------------------------------------------------------------------------------
# MODULES
# ---------------------------------------------------------------------------------------------

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import netCDF4 as nc
from matplotlib.colors import LogNorm, from_levels_and_colors, ListedColormap, LinearSegmentedColormap
from matplotlib import cm
from mpl_toolkits.basemap import Basemap
from datetime import date, timedelta
from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()
from netCDF4 import Dataset
from wrf import (getvar, latlon_coords)
import os

# ---------------------------------------------------------------------------------------------
# MODEL OUTPUT
# ---------------------------------------------------------------------------------------------
# D02
f2015_d02 = '/scratch/leuven/336/vsc33651/nu-wrf-dev/wrfout_nc_files/wrfout_d02_2015.nc'
ds2015_d02 = Dataset(f2015_d02, 'r')
lats_2015_d02 = ds2015_d02.variables['lat'][:]
lons_2015_d02 = ds2015_d02.variables['lon'][:]
time_2015_d02 = ds2015_d02.variables['time'][:]
LTG1_MAX_2015_d02 = ds2015_d02.variables['LTG1_MAX'][:]
LTG2_MAX_2015_d02 = ds2015_d02.variables['LTG2_MAX'][:]
LTG3_MAX_2015_d02 = ds2015_d02.variables['LTG3_MAX'][:]

# D01
f2015_d01 = '/scratch/leuven/336/vsc33651/nu-wrf-dev/wrfout_nc_files/wrfout_d01_2015.nc'
ds2015_d01 = Dataset(f2015_d01, 'r')
lats_2015_d01 = ds2015_d01.variables['lat'][:]
lons_2015_d01 = ds2015_d01.variables['lon'][:]
maxlatd02 = np.amax(lats_2015_d02)
minlatd02 = np.amin(lats_2015_d02)
maxlond02 = np.amax(lons_2015_d02)
minlond02 = np.amin(lons_2015_d02)
latinds_d02= np.unique(np.where((lats_2015_d01 >= minlatd02) & (lats_2015_d01 <= maxlatd02))[0])
loninds_d02 = np.unique(np.where((lons_2015_d01 >= minlond02) & (lons_2015_d01 <= maxlond02))[1])
LTG1_MAX_2015_d01 = ds2015_d01.variables['LTG1_MAX'][:, latinds_d02, loninds_d02]
LTG2_MAX_2015_d01 = ds2015_d01.variables['LTG2_MAX'][:, latinds_d02, loninds_d02]
LTG3_MAX_2015_d01 = ds2015_d01.variables['LTG3_MAX'][:, latinds_d02, loninds_d02]

# ---------------------------------------------------------------------------------------------
# OBSERVATIONS
# ---------------------------------------------------------------------------------------------
# READ IN OBSERVATIONS
file = '/data/leuven/336/vsc33651/projects/nu-wrf-dev/Lightning_data/Slave_lake_full_0.1deg.nc'
ds = nc.Dataset(file, 'r')
# EXTRACT NECESSARY DATA
lats = ds.variables['lat'][:]
lons = ds.variables['lon'][:]
time = ds.variables['time'][:]
# CONVERT TIME TO DATETIME
time_obs = pd.to_datetime(time)+timedelta(minutes=2) ## 2 minutes are added to get correct days.
## Without adding 1 minute, some data was 1 day too early (2 minutes). This is due to the conversion from datetime to int
## and back to datetime. Since the daily data is important, the exact minute of the day is not important.
# SET DOMAIN TO FIT D02
latinds_d02_obs = np.where((lats >= minlatd02) & (lats <= maxlatd02))[0]
loninds_d02_obs = np.where((lons >= minlond02) & (lons <= maxlond02))[0]
# EXTRACT ONLY TIMES IN SUMMER 2015
year = time_obs.year
month = time_obs.month
summer_2015 = np.where([(year == 2015) & (month > 5) & (month < 9)])[1]
# CALCULATE TOTAL NUMBER OF FLASHES (CC+CG)
FD_CC_D02 = ds.variables['Flashdensity_CC'][summer_2015,latinds_d02_obs, loninds_d02_obs]
FD_CG_D02 = ds.variables['Flashdensity_CG'][summer_2015,latinds_d02_obs, loninds_d02_obs]
FD_CC_D02 = np.where(np.isnan(FD_CC_D02),0, FD_CC_D02)
FD_CG_D02 = np.where(np.isnan(FD_CG_D02),0, FD_CG_D02)
CC_CG_D02 = np.add(FD_CC_D02, FD_CG_D02)

# ---------------------------------------------------------------------------------------------
# CALCULATING DAILY AVERAGE THREAT AND DAILY AVERAGE OBSERVED FLASHES
# ---------------------------------------------------------------------------------------------
# OBSERVATIONS AVERAGE PER LOCATION PER DAY
Avg_obs_flashes_pd_2015 = np.mean(CC_CG_D02,axis=0)
# D01
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
# D02
daily_LTG1_MAX_d02 = np.zeros((days,216,363))
daily_LTG2_MAX_d02 = np.zeros((days,216,363))
daily_LTG3_MAX_d02 = np.zeros((days,216,363))
for i in range(0,days):
    j = (i+1)*24
    k = i*24
    daily_LTG1_MAX_d02[i,:,:] = np.mean(LTG1_MAX_2015_d02[k:j,:,:], axis = 0)
    daily_LTG2_MAX_d02[i, :, :] = np.mean(LTG2_MAX_2015_d02[k:j, :, :], axis=0)
    daily_LTG3_MAX_d02[i, :, :] = np.mean(LTG3_MAX_2015_d02[k:j, :, :], axis=0)
# CALCULATE DAILY AVERAGE OF ALL SUMMER
avg_daily_threat3_d01 = np.mean(daily_LTG3_MAX_d01, axis = 0)
avg_daily_threat3_d02 = np.mean(daily_LTG3_MAX_d02, axis = 0)

# ---------------------------------------------------------------------------------------------
# CALCULATING LIGHTNING FLASHES FROM AVERAGE THREAT
# ---------------------------------------------------------------------------------------------
obs_lightning_mean = np.mean(Avg_obs_flashes_pd_2015)
mod_threat_mean_d01 = np.mean(avg_daily_threat3_d01)
mod_threat_mean_d02 = np.mean(avg_daily_threat3_d02)
mod_flashes_d01 = np.divide(np.multiply(obs_lightning_mean, avg_daily_threat3_d01),mod_threat_mean_d01)
mod_flashes_d02 = np.divide(np.multiply(obs_lightning_mean, avg_daily_threat3_d02),mod_threat_mean_d02)

# ---------------------------------------------------------------------------------------------
# MAPS
# ---------------------------------------------------------------------------------------------
switch_interactive = 0
## Map Observed Flashdensity in domain
f = plt.figure(num=None, dpi=700, facecolor='w', edgecolor='k')
m = Basemap(projection= 'lcc', lat_0 = lats_2015_d02.mean(),
            lon_0 = lons_2015_d02.mean(),
            llcrnrlon=(-125.5), llcrnrlat=(57.2),
            urcrnrlon=(-105.4),urcrnrlat = (63),
            resolution = 'f')
m.drawcoastlines(linewidth=0.3)
m.drawcountries()
m.drawparallels(np.arange(-50,70,5))
m.drawmeridians(np.arange(-130,-100,5))
m.drawmapscale(np.amin(lons_2015_d02)+2.5, np.amin(lats_2015_d02)+0.6, lons_2015_d02.mean(), lats_2015_d02.mean(), 100)
m.drawparallels(np.arange(-50,70,5), linewidth=0.5,labels = [True,False,False,False])
m.drawmeridians(np.arange(-130,-100,5), linewidth=0.5,labels = [False,False,False, True])

seismic_mod = cm.get_cmap('seismic',256)
newcmp = ListedColormap(seismic_mod(np.linspace(0.5,1,256)))

# Avg_obs_flashes_pd_2015[Avg_obs_flashes_pd_2015<=0]=10**(-1) # ONLY NEEDED WHEN LOG SCALE

x,y = np.meshgrid(lons[loninds_d02_obs],lats[latinds_d02_obs])
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
m = Basemap(projection= 'lcc', lat_0 = lats_2015_d02.mean(),
            lon_0 = lons_2015_d02.mean(),
            llcrnrlon=(-125.5), llcrnrlat=(57.2),
            urcrnrlon=(-105.4),urcrnrlat = (63),
            resolution = 'f')

m.drawcoastlines(linewidth=0.3)
m.drawcountries()
m.drawparallels(np.arange(-50,70,5))
m.drawmeridians(np.arange(-130,-100,5))
m.drawparallels(np.arange(-50,70,5), linewidth=0.5,labels = [True,False,False,False])
m.drawmeridians(np.arange(-130,-100,5), linewidth=0.5,labels = [False,False,False, True])
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
m = Basemap(projection= 'lcc', lat_0 = lats_2015_d02.mean(),
            lon_0 = lons_2015_d02.mean(),
            llcrnrlon=(-125.5), llcrnrlat=(57.2),
            urcrnrlon=(-105.4),urcrnrlat = (63),
            resolution = 'f')
m.drawcoastlines(linewidth=0.3)
m.drawcountries()
m.drawparallels(np.arange(-50,70,5), linewidth=0.5,labels = [True,False,False,False])
m.drawmeridians(np.arange(-130,-100,5), linewidth=0.5,labels = [False,False,False, True])
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