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

print(np.amax(daily_LPI_d02))
print(np.amax(daily_LPI_d01))
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
    outpath = "/scratch/leuven/336/vsc33651/nu-wrf-dev/Plots"
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
    outpath = "/scratch/leuven/336/vsc33651/nu-wrf-dev/Plots"
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
# TIMESERIES
# ---------------------------------------------------------------------------------------------
# PLOT
timeseries = pd.date_range('2015-05-15', '2020-09-15', freq='D')
data_d01 = np.zeros((len(timeseries), 2))
data_d02 = np.zeros((len(timeseries), 2))
year = timeseries.year
month = timeseries.month
i2015 = np.where((year == 2015) & (month > 5) & (month < 9))
i2016 = np.where((year == 2016) & (month > 5) & (month < 9))
i2017 = np.where((year == 2017) & (month > 5) & (month < 9))
i2018 = np.where((year == 2018) & (month > 5) & (month < 9))
i2019 = np.where((year == 2019) & (month > 5) & (month < 9))
i2020 = np.where((year == 2020) & (month > 5) & (month < 9))
data_d01[:,0] = pd.to_datetime(timeseries)
data_d01[i2015,1] = np.mean(daily_LPI_d01[0:92], axis=(1,2))
data_d01[i2016,1] = np.mean(daily_LPI_d01[92:184], axis=(1,2))
data_d01[i2017,1] = np.mean(daily_LPI_d01[184:276], axis=(1,2))
data_d01[i2018,1] = np.mean(daily_LPI_d01[276:368], axis=(1,2))
data_d01[i2019,1] = np.mean(daily_LPI_d01[368:460], axis=(1,2))
data_d01[i2020,1] = np.mean(daily_LPI_d01[460:552], axis=(1,2))

data_d02[:,0] = pd.to_datetime(timeseries)
data_d02[i2015,1] = np.mean(daily_LPI_d02[0:92], axis=(1,2))
data_d02[i2016,1] = np.mean(daily_LPI_d02[92:184], axis=(1,2))
data_d02[i2017,1] = np.mean(daily_LPI_d02[184:276], axis=(1,2))
data_d02[i2018,1] = np.mean(daily_LPI_d02[276:368], axis=(1,2))
data_d02[i2019,1] = np.mean(daily_LPI_d02[368:460], axis=(1,2))
data_d02[i2020,1] = np.mean(daily_LPI_d02[460:552], axis=(1,2))

data_d01[:,1] = np.where(month == 6, data_d01[:,1],
                np.where(month == 7, data_d01[:,1],
                         np.where(month == 8, data_d01[:,1], np.nan)))
data_d02[:,1] = np.where(month == 6, data_d02[:,1],
                np.where(month == 7, data_d02[:,1],
                         np.where(month == 8, data_d02[:,1], np.nan)))
timeseries_array = timeseries.astype(int)
timeseries_array = np.where((month >5) & (month<9), timeseries_array, np.nan)
timeseries = pd.to_datetime(timeseries_array)
# 14 DAY AVERAGE
# biweeks = pd.date_range('2015-05-15', '2020-09-15', freq='2W')
# biweekly_d01 = np.zeros((len(timeseries), 2))
# biweekly_d02 = np.zeros((len(timeseries), 2))
# for i in range(0,biweeks):
#     j = (i+1)*(14*24)
#     k = i*(14*24)
#     biweekly_d02[i,:,:] = np.mean(daily_LPI_d02[k:j,:,:], axis = 0)
#     biweekly_d01[i, :, :] = np.mean(daily_LPI_d01[k:j, :, :], axis=0)
# avg_biweekly_LPI_d01 = np.mean(biweekly_d01, axis=(1,2))
# avg_biweekly_LPI_d02 = np.mean(biweekly_d02, axis=(1,2))
print(data_d01)
# DOMAIN 1 (from: https://www.semicolonworld.com/question/43500/python-matplotlib-is-there-a-way-to-make-a-discontinuous-axis)
# First: create a subplot for each year
fig, (ax, ax2, ax3, ax4, ax5, ax6) = plt.subplots(1, 6, sharey = 'row')
ax.plot(pd.to_datetime(timeseries), data_d01[:,1],'.',ms=3)
ax2.plot(pd.to_datetime(timeseries), data_d01[:,1],'.',ms=3)
ax3.plot(pd.to_datetime(timeseries), data_d01[:,1],'.',ms=3)
ax4.plot(pd.to_datetime(timeseries), data_d01[:,1],'.',ms=3)
ax5.plot(pd.to_datetime(timeseries), data_d01[:,1],'.',ms=3)
ax6.plot(pd.to_datetime(timeseries), data_d01[:,1],'.',ms=3)
# then set the limits of the x-axis to cover only one year
ax.set_xlim(pd.to_datetime('2015-05-15'), pd.to_datetime('2015-09-15'))
ax2.set_xlim(pd.to_datetime('2016-05-15'), pd.to_datetime('2016-09-15'))
ax3.set_xlim(pd.to_datetime('2017-05-15'), pd.to_datetime('2017-09-15'))
ax4.set_xlim(pd.to_datetime('2018-05-15'), pd.to_datetime('2018-09-15'))
ax5.set_xlim(pd.to_datetime('2019-05-15'), pd.to_datetime('2019-09-15'))
ax6.set_xlim(pd.to_datetime('2020-05-15'), pd.to_datetime('2020-09-15'))
# Set the spines invisible
ax.spines['right'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax3.spines['right'].set_visible(False)
ax4.spines['right'].set_visible(False)
ax5.spines['right'].set_visible(False)
ax2.spines['left'].set_visible(False)
ax3.spines['left'].set_visible(False)
ax4.spines['left'].set_visible(False)
ax5.spines['left'].set_visible(False)
ax6.spines['left'].set_visible(False)
ax.yaxis.tick_left()
# ylim = np.arange(0, np.amax(daily_LPI_d01), 0.25)
ax2.tick_params(axis='y', colors='w')
ax3.tick_params(axis='y', colors='w')
ax4.tick_params(axis='y', colors='w')
ax5.tick_params(axis='y', colors='w')
ax6.tick_params(axis='y', colors='w')
# ax.tick_parms(labeltop='off')
# Make spacing between two axes a bit smaller
plt.subplots_adjust(wspace=0.075)
# add broken axis lines
d = 0.015
kwargs = dict(transform=ax.transAxes, color='k', clip_on = False)
ax.plot((1-d,1+d),(-d,+d), **kwargs)
ax.plot((1-d,1+d),(1-d,1+d), **kwargs)
kwargs.update(transform=ax2.transAxes)
ax2.plot((1-d,1+d),(-d,+d), **kwargs)
ax2.plot((1-d,1+d),(1-d,1+d), **kwargs)
ax2.plot((-d,d),(-d,+d), **kwargs)
ax2.plot((-d,d),(1-d,1+d), **kwargs)
kwargs.update(transform=ax3.transAxes)
ax3.plot((1-d,1+d),(-d,+d), **kwargs)
ax3.plot((1-d,1+d),(1-d,1+d), **kwargs)
ax3.plot((-d,d),(-d,+d), **kwargs)
ax3.plot((-d,d),(1-d,1+d), **kwargs)
kwargs.update(transform=ax4.transAxes)
ax4.plot((1-d,1+d),(-d,+d), **kwargs)
ax4.plot((1-d,1+d),(1-d,1+d), **kwargs)
ax4.plot((-d,d),(-d,+d), **kwargs)
ax4.plot((-d,d),(1-d,1+d), **kwargs)
kwargs.update(transform=ax5.transAxes)
ax5.plot((1-d,1+d),(-d,+d), **kwargs)
ax5.plot((1-d,1+d),(1-d,1+d), **kwargs)
ax5.plot((-d,d),(-d,+d), **kwargs)
ax5.plot((-d,d),(1-d,1+d), **kwargs)
kwargs.update(transform=ax6.transAxes)
ax6.plot((-d,d),(-d,+d), **kwargs)
ax6.plot((-d,d),(1-d,1+d), **kwargs)
ax.set_ylabel('LPI (J kg$^{-1}$)', size = 12)
fig.suptitle('Convection-parameterized (9 km)')
ax.ticklabel_format(axis='y',style='sci', scilimits=(0,0))
ax.grid(which='major', axis='both', color='lightgray')
ax2.grid(which='major', axis='both', color='lightgray')
ax3.grid(which='major', axis='both', color='lightgray')
ax4.grid(which='major', axis='both', color='lightgray')
ax5.grid(which='major', axis='both', color='lightgray')
ax6.grid(which='major', axis='both', color='lightgray')

ax.set_xticks(ticks=[pd.to_datetime('2015-06', format='%Y-%m'), pd.to_datetime('2015-07', format='%Y-%m'),
                     pd.to_datetime('2015-08', format='%Y-%m'), pd.to_datetime('2015-09', format='%Y-%m')])
ax.set_xticklabels(['Jun \'15', 'Jul \'15', 'Aug \'15', 'Sep \'15'], rotation=30)
ax2.set_xticks(ticks=[pd.to_datetime('2016-06', format='%Y-%m'), pd.to_datetime('2016-07', format='%Y-%m'),
                      pd.to_datetime('2016-08', format='%Y-%m'), pd.to_datetime('2016-09', format='%Y-%m')])
ax2.set_xticklabels(['Jun \'16', 'Jul \'16', 'Aug \'16', 'Sep \'16'], rotation=30)
ax3.set_xticks(ticks=[pd.to_datetime('2017-06', format='%Y-%m'), pd.to_datetime('2017-07', format='%Y-%m'),
                      pd.to_datetime('2017-08', format='%Y-%m'), pd.to_datetime('2017-09', format='%Y-%m')])
ax3.set_xticklabels(['Jun \'17', 'Jul \'17', 'Aug \'17', 'Sep \'17'], rotation=30)
ax4.set_xticks(ticks=[pd.to_datetime('2018-06', format='%Y-%m'), pd.to_datetime('2018-07', format='%Y-%m'),
                      pd.to_datetime('2018-08', format='%Y-%m'), pd.to_datetime('2018-09', format='%Y-%m')])
ax4.set_xticklabels(['Jun \'18', 'Jul \'18', 'Aug \'18', 'Sep \'18'], rotation=30)
ax5.set_xticks(ticks=[pd.to_datetime('2019-06', format='%Y-%m'), pd.to_datetime('2019-07', format='%Y-%m'),
                      pd.to_datetime('2019-08', format='%Y-%m'), pd.to_datetime('2019-09', format='%Y-%m')])
ax5.set_xticklabels(['Jun \'19', 'Jul \'19', 'Aug \'19', 'Sep \'19'], rotation=30)
ax6.set_xticks(ticks=[pd.to_datetime('2020-06', format='%Y-%m'), pd.to_datetime('2020-07', format='%Y-%m'),
                      pd.to_datetime('2020-08', format='%Y-%m'), pd.to_datetime('2020-09', format='%Y-%m')])
ax6.set_xticklabels(['Jun \'20', 'Jul \'20', 'Aug \'20', 'Sep \'20'], rotation=30)
# plt.xticks(ticks=[pd.to_datetime('2015-06', format='%Y-%m'), pd.to_datetime('2015-07', format='%Y-%m'), pd.to_datetime('2015-08', format='%Y-%m'),
#            pd.to_datetime('2015-09', format='%Y-%m'), pd.to_datetime('2016-06', format='%Y-%m'), pd.to_datetime('2016-07', format='%Y-%m'),
#            pd.to_datetime('2016-08', format='%Y-%m'), pd.to_datetime('2016-09', format='%Y-%m'), pd.to_datetime('2017-06', format='%Y-%m'),
#            pd.to_datetime('2017-07', format='%Y-%m'), pd.to_datetime('2017-08', format='%Y-%m'), pd.to_datetime('2017-09', format='%Y-%m'),
#            pd.to_datetime('2018-06', format='%Y-%m'), pd.to_datetime('2018-07', format='%Y-%m'), pd.to_datetime('2018-08', format='%Y-%m'),
#            pd.to_datetime('2018-09', format='%Y-%m'), pd.to_datetime('2019-06', format='%Y-%m'), pd.to_datetime('2019-07', format='%Y-%m'),
#            pd.to_datetime('2019-08', format='%Y-%m'), pd.to_datetime('2020-06', format='%Y-%m'), pd.to_datetime('2020-07', format='%Y-%m'),
#            pd.to_datetime('2020-08', format='%Y-%m'), pd.to_datetime('2020-09', format='%Y-%m')],
#            labels = ['Jun \'15', 'Jul \'15', 'Aug \'15', 'Sep \'15', 'Jun \'16', 'Jul \'16', 'Aug \'16', 'Sep \'16',
#             'Jun \'17', 'Jul \'17', 'Aug \'17', 'Sep \'17', 'Jun \'18', 'Jul \'18', 'Aug \'18', 'Sep \'18',
#             'Jun \'19', 'Jul \'19', 'Aug \'19', 'Sep \'19', 'Jun \'20', 'Jul \'20', 'Aug \'20', 'Sep \'20'])
plt.show()

# DOMAIN 2
fig, (ax, ax2, ax3, ax4, ax5, ax6) = plt.subplots(1, 6, sharey = 'row')
ax.plot(pd.to_datetime(timeseries), data_d02[:,1],'.',ms=3)
ax2.plot(pd.to_datetime(timeseries), data_d02[:,1],'.',ms=3)
ax3.plot(pd.to_datetime(timeseries), data_d02[:,1],'.',ms=3)
ax4.plot(pd.to_datetime(timeseries), data_d02[:,1],'.',ms=3)
ax5.plot(pd.to_datetime(timeseries), data_d02[:,1],'.',ms=3)
ax6.plot(pd.to_datetime(timeseries), data_d02[:,1],'.',ms=3)
# then set the limits of the x-axis to cover only one year
ax.set_xlim(pd.to_datetime('2015-05-15'), pd.to_datetime('2015-09-15'))
ax2.set_xlim(pd.to_datetime('2016-05-15'), pd.to_datetime('2016-09-15'))
ax3.set_xlim(pd.to_datetime('2017-05-15'), pd.to_datetime('2017-09-15'))
ax4.set_xlim(pd.to_datetime('2018-05-15'), pd.to_datetime('2018-09-15'))
ax5.set_xlim(pd.to_datetime('2019-05-15'), pd.to_datetime('2019-09-15'))
ax6.set_xlim(pd.to_datetime('2020-05-15'), pd.to_datetime('2020-09-15'))
# Set the spines invisible
ax.spines['right'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax3.spines['right'].set_visible(False)
ax4.spines['right'].set_visible(False)
ax5.spines['right'].set_visible(False)
ax2.spines['left'].set_visible(False)
ax3.spines['left'].set_visible(False)
ax4.spines['left'].set_visible(False)
ax5.spines['left'].set_visible(False)
ax6.spines['left'].set_visible(False)
ax.yaxis.tick_left()
# ylim = np.arange(0, np.amax(daily_LPI_d01), 0.25)
ax2.tick_params(axis='y', colors='w')
ax3.tick_params(axis='y', colors='w')
ax4.tick_params(axis='y', colors='w')
ax5.tick_params(axis='y', colors='w')
ax6.tick_params(axis='y', colors='w')
# ax.tick_parms(labeltop='off')
# Make spacing between two axes a bit smaller
plt.subplots_adjust(wspace=0.075)
# add broken axis lines
d = 0.015
kwargs = dict(transform=ax.transAxes, color='k', clip_on = False)
ax.plot((1-d,1+d),(-d,+d), **kwargs)
ax.plot((1-d,1+d),(1-d,1+d), **kwargs)
kwargs.update(transform=ax2.transAxes)
ax2.plot((1-d,1+d),(-d,+d), **kwargs)
ax2.plot((1-d,1+d),(1-d,1+d), **kwargs)
ax2.plot((-d,d),(-d,+d), **kwargs)
ax2.plot((-d,d),(1-d,1+d), **kwargs)
kwargs.update(transform=ax3.transAxes)
ax3.plot((1-d,1+d),(-d,+d), **kwargs)
ax3.plot((1-d,1+d),(1-d,1+d), **kwargs)
ax3.plot((-d,d),(-d,+d), **kwargs)
ax3.plot((-d,d),(1-d,1+d), **kwargs)
kwargs.update(transform=ax4.transAxes)
ax4.plot((1-d,1+d),(-d,+d), **kwargs)
ax4.plot((1-d,1+d),(1-d,1+d), **kwargs)
ax4.plot((-d,d),(-d,+d), **kwargs)
ax4.plot((-d,d),(1-d,1+d), **kwargs)
kwargs.update(transform=ax5.transAxes)
ax5.plot((1-d,1+d),(-d,+d), **kwargs)
ax5.plot((1-d,1+d),(1-d,1+d), **kwargs)
ax5.plot((-d,d),(-d,+d), **kwargs)
ax5.plot((-d,d),(1-d,1+d), **kwargs)
kwargs.update(transform=ax6.transAxes)
ax6.plot((-d,d),(-d,+d), **kwargs)
ax6.plot((-d,d),(1-d,1+d), **kwargs)
ax.set_ylabel('LPI (J kg$^{-1}$)', size = 12)
fig.suptitle('Convection-permitting (3 km)')
ax.ticklabel_format(axis='y',style='sci', scilimits=(0,0))
ax.grid(which='major', axis='both', color='lightgray')
ax2.grid(which='major', axis='both', color='lightgray')
ax3.grid(which='major', axis='both', color='lightgray')
ax4.grid(which='major', axis='both', color='lightgray')
ax5.grid(which='major', axis='both', color='lightgray')
ax6.grid(which='major', axis='both', color='lightgray')

ax.set_xticks(ticks=[pd.to_datetime('2015-06', format='%Y-%m'), pd.to_datetime('2015-07', format='%Y-%m'),
                     pd.to_datetime('2015-08', format='%Y-%m'), pd.to_datetime('2015-09', format='%Y-%m')])
ax.set_xticklabels(['Jun \'15', 'Jul \'15', 'Aug \'15', 'Sep \'15'], rotation=30)
ax2.set_xticks(ticks=[pd.to_datetime('2016-06', format='%Y-%m'), pd.to_datetime('2016-07', format='%Y-%m'),
                      pd.to_datetime('2016-08', format='%Y-%m'), pd.to_datetime('2016-09', format='%Y-%m')])
ax2.set_xticklabels(['Jun \'16', 'Jul \'16', 'Aug \'16', 'Sep \'16'], rotation=30)
ax3.set_xticks(ticks=[pd.to_datetime('2017-06', format='%Y-%m'), pd.to_datetime('2017-07', format='%Y-%m'),
                      pd.to_datetime('2017-08', format='%Y-%m'), pd.to_datetime('2017-09', format='%Y-%m')])
ax3.set_xticklabels(['Jun \'17', 'Jul \'17', 'Aug \'17', 'Sep \'17'], rotation=30)
ax4.set_xticks(ticks=[pd.to_datetime('2018-06', format='%Y-%m'), pd.to_datetime('2018-07', format='%Y-%m'),
                      pd.to_datetime('2018-08', format='%Y-%m'), pd.to_datetime('2018-09', format='%Y-%m')])
ax4.set_xticklabels(['Jun \'18', 'Jul \'18', 'Aug \'18', 'Sep \'18'], rotation=30)
ax5.set_xticks(ticks=[pd.to_datetime('2019-06', format='%Y-%m'), pd.to_datetime('2019-07', format='%Y-%m'),
                      pd.to_datetime('2019-08', format='%Y-%m'), pd.to_datetime('2019-09', format='%Y-%m')])
ax5.set_xticklabels(['Jun \'19', 'Jul \'19', 'Aug \'19', 'Sep \'19'], rotation=30)
ax6.set_xticks(ticks=[pd.to_datetime('2020-06', format='%Y-%m'), pd.to_datetime('2020-07', format='%Y-%m'),
                      pd.to_datetime('2020-08', format='%Y-%m'), pd.to_datetime('2020-09', format='%Y-%m')])
ax6.set_xticklabels(['Jun \'20', 'Jul \'20', 'Aug \'20', 'Sep \'20'], rotation=30)
plt.show()

# ---------------------------------------------------------------------------------------------
# SEASONALITY 2.0
# ---------------------------------------------------------------------------------------------
new_data_d02 = np.zeros((546)) # Because 78 weeks in total over 6 years (13/year) -> 546 days (last day is dropped)
new_data_d01 = np.zeros((546)) # Because 78 weeks in total over 6 years (13/year) -> 546 days (last day is dropped)

new_data_d02[0:91] = data_d02[17:108,1]
new_data_d02[91:182] = data_d02[383:474,1]
new_data_d02[182:273] = data_d02[748:839,1]
new_data_d02[273:364] = data_d02[1113:1204,1]
new_data_d02[364:455] = data_d02[1478:1569,1]
new_data_d02[455:546] = data_d02[1844:1935,1]

new_data_d01[0:91] = data_d01[17:108,1]
new_data_d01[91:182] = data_d01[383:474,1]
new_data_d01[182:273] = data_d01[748:839,1]
new_data_d01[273:364] = data_d01[1113:1204,1]
new_data_d01[364:455] = data_d01[1478:1569,1]
new_data_d01[455:546] = data_d01[1844:1935,1]

weekly_LPI_d02 = np.zeros((78)) # Because 78 weeks in total over 6 years
weekly_LPI_d01 = np.zeros((78))
for i in range(0,78):
    j = (i+1)*7
    k = i*7
    weekly_LPI_d02[i] = np.nanmean(new_data_d02[k:j])
    weekly_LPI_d01[i] = np.nanmean(new_data_d01[k:j])

seasonal_LPI_d01 = np.zeros((14))
seasonal_LPI_d02 = np.zeros((14))
for i in range(0,13):
    indices = np.arange(i,78,13)
    seasonal_LPI_d02[i] = np.nanmean(weekly_LPI_d02[indices])
    seasonal_LPI_d01[i] = np.nanmean(weekly_LPI_d01[indices])

seasonal_LPI_d01[13] = np.nan
seasonal_LPI_d02[13] = np.nan
print(seasonal_LPI_d01)
weeks = pd.date_range('2015-06-01','2015-09-01',freq='W-MON')

plt.plot(weeks, seasonal_LPI_d01, 'k')
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
plt.ylim(0,0.0018)
plt.show()

plt.plot(weeks, seasonal_LPI_d02, 'k')
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
plt.ylim(0,0.010)
# plt.show()

# ---------------------------------------------------------------------------------------------
# FREQUENCY IFO HOURLY GRID FLASH DENSITY (SEE WONG ET AL. FIG. 5)
# ---------------------------------------------------------------------------------------------
# rounded_LPI_d02 = np.round(LPI_d02, 8)
unique,counts=np.unique(LPI_d02,return_counts=True)

# rolling_average = np.convolve(counts, np.ones(10), 'valid')/10
plt.plot(unique,counts,'.', ms=1)
plt.ylabel('Frequency')
plt.xlabel('Hourly gird flash density (#/km${2}$')
plt.xscale('log')
plt.yscale('log')
plt.show()
