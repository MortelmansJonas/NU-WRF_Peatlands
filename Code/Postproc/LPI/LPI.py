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
import matplotlib.dates as mdates

# ---------------------------------------------------------------------------------------------
# LOAD DATA OF D02
# ---------------------------------------------------------------------------------------------
# 2015
ds_lpi_2015_d02 = Dataset('/scratch/leuven/336/vsc33651/nu-wrf-dev/wrfout_nc_files/lpi_time_d02_2015_v03.nc','r')
time_2015 = ds_lpi_2015_d02['time'][:]
lat_d02 = ds_lpi_2015_d02['lat'][40:-40,40:-40]
lon_d02 = ds_lpi_2015_d02['lon'][40:-40,40:-40]
lpi_2015_d02 = ds_lpi_2015_d02['LPI'][:,40:-40,40:-40]
time_2015_pd = pd.to_datetime('2000010100',format='%Y%m%d%H') + pd.to_timedelta(time_2015, unit='h')

# 2016
ds_lpi_2016_d02 = Dataset('/scratch/leuven/336/vsc33651/nu-wrf-dev/wrfout_nc_files/lpi_time_d02_2016_v03.nc','r')
time_2016 = ds_lpi_2016_d02['time'][:]
lpi_2016_d02 = ds_lpi_2016_d02['LPI'][:,40:-40,40:-40]
time_2016_pd = pd.to_datetime('2000010100',format='%Y%m%d%H') + pd.to_timedelta(time_2016, unit='h')

# 2017
ds_lpi_2017_d02 = Dataset('/scratch/leuven/336/vsc33651/nu-wrf-dev/wrfout_nc_files/lpi_time_d02_2017_v03.nc','r')
time_2017 = ds_lpi_2017_d02['time'][:]
lpi_2017_d02 = ds_lpi_2017_d02['LPI'][:,40:-40,40:-40]
time_2017_pd = pd.to_datetime('2000010100',format='%Y%m%d%H') + pd.to_timedelta(time_2017, unit='h')

# 2018
ds_lpi_2018_d02 = Dataset('/scratch/leuven/336/vsc33651/nu-wrf-dev/wrfout_nc_files/lpi_time_d02_2018_v03.nc','r')
time_2018 = ds_lpi_2018_d02['time'][:]
lpi_2018_d02 = ds_lpi_2018_d02['LPI'][:,40:-40,40:-40]
time_2018_pd = pd.to_datetime('2000010100',format='%Y%m%d%H') + pd.to_timedelta(time_2018, unit='h')

# 2019
ds_lpi_2019_d02 = Dataset('/scratch/leuven/336/vsc33651/nu-wrf-dev/wrfout_nc_files/lpi_time_d02_2019_v03.nc','r')
time_2019 = ds_lpi_2019_d02['time'][:]
lpi_2019_d02 = ds_lpi_2019_d02['LPI'][:,40:-40,40:-40]
time_2019_pd = pd.to_datetime('2000010100',format='%Y%m%d%H') + pd.to_timedelta(time_2019, unit='h')

# 2020
ds_lpi_2020_d02 = Dataset('/scratch/leuven/336/vsc33651/nu-wrf-dev/wrfout_nc_files/lpi_time_d02_2020_v03.nc','r')
time_2020 = ds_lpi_2020_d02['time'][:]
lpi_2020_d02 = ds_lpi_2020_d02['LPI'][:,40:-40,40:-40]
time_2020_pd = pd.to_datetime('2000010100',format='%Y%m%d%H') + pd.to_timedelta(time_2020, unit='h')

# ---------------------------------------------------------------------------------------------
# LOAD DATA OF D01
# ---------------------------------------------------------------------------------------------
# 2015
ds_lpi_2015_d01 = Dataset('/scratch/leuven/336/vsc33651/nu-wrf-dev/wrfout_nc_files/lpi_time_d01_2015_v03.nc','r')
# Get only area of domain 2:
latinds = np.unique(np.where((ds_lpi_2015_d01['lat'] >= np.amin(lat_d02)) & (ds_lpi_2015_d01['lat'] <= np.amax(lat_d02)))[0])
loninds = np.unique(np.where((ds_lpi_2015_d01['lon'] >= np.amin(lon_d02)) & (ds_lpi_2015_d01['lon'] <= np.amax(lon_d02)))[1])
lat_d01 = ds_lpi_2015_d01['lat'][latinds,loninds]
lon_d01 = ds_lpi_2015_d01['lon'][latinds,loninds]
lpi_2015_d01 = ds_lpi_2015_d01['LPI'][:,latinds,loninds]

# 2016
ds_lpi_2016_d01 = Dataset('/scratch/leuven/336/vsc33651/nu-wrf-dev/wrfout_nc_files/lpi_time_d01_2016_v03.nc','r')
lpi_2016_d01 = ds_lpi_2016_d01['LPI'][:,latinds,loninds]

# 2017
ds_lpi_2017_d01 = Dataset('/scratch/leuven/336/vsc33651/nu-wrf-dev/wrfout_nc_files/lpi_time_d01_2017_v03.nc','r')
lpi_2017_d01 = ds_lpi_2017_d01['LPI'][:,latinds,loninds]

# 2018
ds_lpi_2018_d01 = Dataset('/scratch/leuven/336/vsc33651/nu-wrf-dev/wrfout_nc_files/lpi_time_d01_2018_v03.nc','r')
lpi_2018_d01 = ds_lpi_2018_d01['LPI'][:,latinds,loninds]

# 2019
ds_lpi_2019_d01 = Dataset('/scratch/leuven/336/vsc33651/nu-wrf-dev/wrfout_nc_files/lpi_time_d01_2019_v03.nc','r')
lpi_2019_d01 = ds_lpi_2019_d01['LPI'][:,latinds,loninds]

# 2020
ds_lpi_2020_d01 = Dataset('/scratch/leuven/336/vsc33651/nu-wrf-dev/wrfout_nc_files/lpi_time_d01_2020_v03.nc','r')
lpi_2020_d01 = ds_lpi_2020_d01['LPI'][:,latinds,loninds]

# ---------------------------------------------------------------------------------------------
# PUT EVERYTHING TOGETHER
# ---------------------------------------------------------------------------------------------
# TIME
time_all = np.zeros((13248))
time_all[0:2208] = time_2015_pd
time_all[2208:4416] = time_2016_pd
time_all[4416:6624] = time_2017_pd
time_all[6624:8832] = time_2018_pd
time_all[8832:11040] = time_2019_pd
time_all[11040:13248] = time_2020_pd

# LPI DOMAIN 2
LPI_d02 = np.zeros((13248,217,364))
LPI_d02[0:2208,:,:] = lpi_2015_d02
LPI_d02[2208:4416] = lpi_2016_d02
LPI_d02[4416:6624] = lpi_2017_d02
LPI_d02[6624:8832] = lpi_2018_d02
LPI_d02[8832:11040] = lpi_2019_d02
LPI_d02[11040:13248] = lpi_2020_d02
# LPI DOMAIN 1
LPI_d01 = np.zeros((13248,90,157))
LPI_d01[0:2208,:,:] = lpi_2015_d01
LPI_d01[2208:4416] = lpi_2016_d01
LPI_d01[4416:6624] = lpi_2017_d01
LPI_d01[6624:8832] = lpi_2018_d01
LPI_d01[8832:11040] = lpi_2019_d01
LPI_d01[11040:13248] = lpi_2020_d01
# ---------------------------------------------------------------------------------------------
# CALCULATE DAILY AVERAGE LPI
# ---------------------------------------------------------------------------------------------
days = np.divide(13248,24).astype(int)
daily_LPI_d02 = np.zeros((days,217,364))
daily_LPI_d01 = np.zeros((days,90,157))
for i in range(0,days):
    j = (i+1)*24
    k = i*24
    daily_LPI_d02[i,:,:] = np.mean(LPI_d02[k:j,:,:], axis = 0)
    daily_LPI_d01[i, :, :] = np.mean(LPI_d01[k:j, :, :], axis=0)
avg_daily_LPI_d02 = np.mean(daily_LPI_d02, axis = 0)
avg_daily_LPI_d01 = np.mean(daily_LPI_d01, axis = 0)

# ---------------------------------------------------------------------------------------------
# MAPS
# ---------------------------------------------------------------------------------------------
switch_interactive = 0
seismic_mod = cm.get_cmap('seismic',256)
newcmp = ListedColormap(seismic_mod(np.linspace(0.5,1,256)))
# MAP
# DOMAIN 1
f = plt.figure(num=None, dpi=700, facecolor='w', edgecolor='k')
m = Basemap(projection= 'lcc', lat_0 = lat_d01.mean(),
            lon_0 = lon_d01.mean(),
            llcrnrlon=(-125.5), llcrnrlat=(57.2),
            urcrnrlon=(-105.4),urcrnrlat = (63),
            resolution = 'f')
m.drawcoastlines(linewidth=0.3)
m.drawcountries()
m.drawparallels(np.arange(-50,70,5))
m.drawmeridians(np.arange(-130,-100,5))
m.drawmapscale(-123, 57.8, lon_d01.mean(), lat_d01.mean(), 100)
m.drawparallels(np.arange(-50,70,5), linewidth=0.5,labels = [True,False,False,False])
m.drawmeridians(np.arange(-130,-100,5), linewidth=0.5,labels = [False,False,False, True])
m.pcolormesh(lon_d01,lat_d01,avg_daily_LPI_d01,latlon=True,cmap=newcmp)

if switch_interactive=="1":
    plt.show()
else:
    outpath = "/scratch/leuven/317/vsc31786/nu-wrf-dev/Plots"
    fname = "LPI_ALL_d01_v03"
    fname_long = os.path.join(outpath, fname+'.png')
    plt.title('Convection-parameterized (9 km)',size=12)
    plt.colorbar(label='LPI (J kg$^{-1}$)',extend = 'max',shrink = 0.6)
    plt.savefig(fname_long, dpi=f.dpi)
    plt.close()

# DOMAIN 2
f = plt.figure(num=None, dpi=700, facecolor='w', edgecolor='k')
m = Basemap(projection= 'lcc', lat_0 = lat_d02.mean(),
            lon_0 = lon_d02.mean(),
            llcrnrlon=(-125.5), llcrnrlat=(57.2),
            urcrnrlon=(-105.4),urcrnrlat = (63),
            resolution = 'f')
m.drawcoastlines(linewidth=0.3)
m.drawcountries()
m.drawparallels(np.arange(-50,70,5))
m.drawmeridians(np.arange(-130,-100,5))
m.drawmapscale(-123, 57.8, lon_d02.mean(), lat_d02.mean(), 100)
m.drawparallels(np.arange(-50,70,5), linewidth=0.5,labels = [True,False,False,False])
m.drawmeridians(np.arange(-130,-100,5), linewidth=0.5,labels = [False,False,False, True])
m.pcolormesh(lon_d02,lat_d02,avg_daily_LPI_d02,latlon=True,cmap=newcmp)

if switch_interactive=="1":
    plt.show()
else:
    outpath = "/scratch/leuven/317/vsc31786/nu-wrf-dev/Plots"
    fname = "LPI_ALL_d02_v03"
    fname_long = os.path.join(outpath, fname+'.png')
    plt.title('Convection-permitting (3 km)',size=12)
    plt.colorbar(label='LPI (J kg$^{-1}$)',extend = 'max',shrink = 0.6)
    plt.savefig(fname_long, dpi=f.dpi)
    plt.close()

# ---------------------------------------------------------------------------------------------
# DIURNAL CYCLE
# ---------------------------------------------------------------------------------------------
# GET THE TIMES RIGHT
# time = pd.to_datetime('2000010100',format='%Y%m%d%H') + pd.to_timedelta(time_all, unit='h')
times = pd.to_datetime(time_all) - pd.to_timedelta(7, unit='h')
months = times.month
summer = np.where([(months > 5) & (months < 9)])
time_summer = times[summer[1]]

# EXTRACT THE RIGHT DATA
# lpi_summer_d02 = LPI_d02[summer[1],:,:]
# lpi_summer_d01 = LPI_d01[summer[1],:,:]

# lpi_diurnal_d02 = np.zeros((24,217,364))
# lpi_diurnal_d01 = np.zeros((24,90,157))
# for i in range(0,24):
#     time_select_vector = np.where(time_summer.hour == i)
#     lpi_diurnal_d02[i, :, :] = np.mean(lpi_summer_d02[time_select_vector[0], :,:], axis=0)
#     lpi_diurnal_d01[i, :, :] = np.mean(lpi_summer_d01[time_select_vector[0], :, :], axis=0)

# PLOT DOMAIN 1
# hour_of_day = np.arange(0,24,1)
# plt.plot(hour_of_day, np.mean(lpi_diurnal_d01,axis=(1,2)), 'k')
# plt.xlabel('Hour of the day', size = 12)
# plt.ylabel('LPI (J kg$^{-1}$)', size = 12)
# plt.title('Convection-parameterized (9 km)')
# plt.ticklabel_format(axis='y',style='sci', scilimits=(0,0))
# plt.show()

# PLOT DOMAIN 2
# plt.plot(hour_of_day, np.mean(lpi_diurnal_d02,axis=(1,2)), 'k')
# plt.xlabel('Hour of the day', size = 12)
# plt.ylabel('LPI (J kg$^{-1}$)', size = 12)
# plt.title('Convection-permitting (3 km)')
# plt.ticklabel_format(axis='y',style='sci', scilimits=(0,0))
# plt.show()

# ---------------------------------------------------------------------------------------------
# SEASONAL CLIMATOLOGY
# ---------------------------------------------------------------------------------------------
weeks = pd.date_range(time_2015_pd[0],time_2015_pd[2207],freq='W-MON')
seasonal_d02 = np.zeros((len(weeks),217,364))
seasonal_d01 = np.zeros((len(weeks),90,157))
for i in range(0,len(weeks)): # By far not the most efficient way, but the only way I could make it work.
    # MB: end of index 7 days times 24 hours
    j = i+(7*24)
    indices = np.array([np.arange(i,j), np.arange(i+2208,j+2208), np.arange(i+4416,j+4416), np.arange(i+6624,j+6624),
                        np.arange(i+8832,j+8832), np.arange(i+11040,j+11040)]).flatten()
    seasonal_d02[i,:,:] = np.mean(LPI_d02[indices], axis=(0))
    seasonal_d01[i,:,:] = np.mean(LPI_d01[indices], axis=(0))

# PLOT DOMAIN 1
plt.plot(weeks, np.mean(seasonal_d01, axis=(1,2)), 'k')
plt.xlabel('Date', size = 12)
plt.ylabel('LPI (J kg$^{-1}$)', size = 12)
plt.title('Convection-parameterized (9 km)')
plt.ticklabel_format(axis='y',style='sci', scilimits=(0,0))
locator = mdates.MonthLocator()
fmt = mdates.DateFormatter('%b')
X = plt.gca().xaxis
X.set_major_locator(locator)
X.set_major_formatter(fmt)
plt.grid(which='major', axis='both', color='lightgray')
plt.show()

# PLOT DOMAIN 2
plt.plot(weeks,np.mean(seasonal_d02,axis=(1,2)), 'k')
plt.xlabel('Date', size = 12)
plt.ylabel('LPI (J kg$^{-1}$)', size = 12)
plt.title('Convection-permitting (3 km)')
plt.ticklabel_format(axis='y',style='sci', scilimits=(0,0))
locator = mdates.MonthLocator()
fmt = mdates.DateFormatter('%b')
X = plt.gca().xaxis
X.set_major_locator(locator)
X.set_major_formatter(fmt)
plt.grid(which='major', axis='both', color='lightgray')
plt.show()

# ---------------------------------------------------------------------------------------------
# TIMESERIES
# ---------------------------------------------------------------------------------------------
time_all = pd.to_datetime(time_all)
# DOMAIN 1
plt.plot(time_all, np.mean(LPI_d01, axis = (1,2)),'.')
plt.xlabel('Date', size = 12)
plt.ylabel('LPI (J kg$^{-1}$)', size = 12)
plt.title('Convection-parameterized (9 km)')
plt.ticklabel_format(axis='y',style='sci', scilimits=(0,0))
plt.show()

# DOMAIN 2
plt.plot(time_all, np.mean(LPI_d02, axis = (1,2)),'.')
plt.xlabel('Date', size = 12)
plt.ylabel('LPI (J kg$^{-1}$)', size = 12)
plt.title('Convection-permitting (3 km)')
plt.ticklabel_format(axis='y',style='sci', scilimits=(0,0))
plt.show()
