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
from datetime import timedelta
from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()
from netCDF4 import Dataset
import os
from scipy.stats import sem
import numpy.matlib

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

# Select latitudes and longitudes that are at least 40 grid cells away from border
# un_lat = np.unique(lats_2015_d02)
# un_lon = np.unique(lons_2015_d02)
# un_lat = np.sort(un_lat)
# un_lon = np.sort(un_lon)
# lat_sel= numpy.matlib.repmat(un_lat,76777,1)
# lon_sel=numpy.matlib.repmat(un_lon,40088,1)

lat_sel= lats_2015_d02[40:-40,40:-40]
lon_sel= lons_2015_d02[40:-40,40:-40]
LTG3_MAX_2015_spinup = LTG3_MAX_2015_d02[:,40:-40,40:-40]
daily_LTG3_MAX_d02 = np.zeros((92,136,283))
for i in range(0,92):
    j = (i+1)*24
    k = i*24
    daily_LTG3_MAX_d02[i, :, :] = np.mean(LTG3_MAX_2015_spinup[k:j, :, :], axis=0)
# CALCULATE DAILY AVERAGE OF ALL SUMMER
avg_daily_threat3_d02 = np.mean(daily_LTG3_MAX_d02, axis = 0)
# ---------------------------------------------------------------------------------------------
# MAP MODEL OUTPUT
# ---------------------------------------------------------------------------------------------
# CREATE NEW COLORBAR WITH 0 = WHITE
seismic_mod = cm.get_cmap('seismic',256)
newcmp = ListedColormap(seismic_mod(np.linspace(0.5,1,256)))
# CREATE MAP
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
m.pcolormesh(lon_sel,lat_sel,avg_daily_threat3_d02 ,latlon=True,cmap=newcmp)
switch_interactive = 0

if switch_interactive=="1":
    plt.show()
else:
    outpath = "/scratch/leuven/336/vsc33651/nu-wrf-dev/Plots"
    fname = "avg_lf_pd_2015_d02_spinup"
    fname_long = os.path.join(outpath, fname+'.png')
    plt.title('Convection-permitting (3-km)', size = 12)
    plt.colorbar(label='Number of flashes (day$^{-1}$ km$^{-2}$)', extend='max', shrink = 0.6)
    # plt.clim(10**(-1),10**2) # Only needed if logarithmic scale is desired
    # plt.clim(0,4)
    plt.savefig(fname_long, dpi=f.dpi)
    plt.close()

# ---------------------------------------------------------------------------------------------
# PLOTS DIURNAL CYCLE
# ---------------------------------------------------------------------------------------------
# DOMAIN 2

ds_diurnal_d02 = Dataset('/scratch/leuven/336/vsc33651/nu-wrf-dev/wrfout_nc_files/2015_d02_diurnal.nc', mode='r')
hour_of_day = ds_diurnal_d02['time'][:]
diurnal_LTG3_d02 = ds_diurnal_d02['LTG3_MAX'][:]
plt.plot(hour_of_day, np.mean(diurnal_LTG3_d02,axis=(1,2)), 'k')
#plt.fill_between(hour_of_day, np.mean(diurnal_LTG3_d02,axis=(1,2))-np.std(diurnal_LTG3_d02, axis=(1,2)),
#                 np.mean(diurnal_LTG3_d02,axis=(1,2))+np.std(diurnal_LTG3_d02, axis=(1,2)))
plt.xlabel('Hour of the day', size = 12)
plt.ylabel('Lightning threat (day$^{-1}$ (3km)$^{-2}$)', size = 12)
plt.ticklabel_format(axis='y',style='sci', scilimits=(0,0))
plt.show()