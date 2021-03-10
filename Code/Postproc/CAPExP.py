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
import os

ds_CP = Dataset('/scratch/leuven/336/vsc33651/nu-wrf-dev/wrfout_nc_files/2015_d01_CP.nc', mode='r')
time = ds_CP['time'][:]
lats_2015_d02 = ds_CP['lat'][:]
lons_2015_d02 = ds_CP['lon'][:]
CAPExP_d02 = ds_CP['CAPExP'][:]
maxlatd02 = np.amax(lats_2015_d02)
minlatd02 = np.amin(lats_2015_d02)
maxlond02 = np.amax(lons_2015_d02)
minlond02 = np.amin(lons_2015_d02)
meanlatd02 = (maxlatd02 + minlatd02)/2
meanlond02 = (maxlond02 + minlond02)/2

daily_avg_CAPExP_d02 = np.zeros((92,140,200))
for i in range(0,92):
    j = (i+1)*24
    k = i*24
    daily_avg_CAPExP_d02[i,:,:] = np.mean(CAPExP_d02[k:j,:,:], axis=0)
avg_daily_CAPExP = np.mean(daily_avg_CAPExP_d02, axis = 0)

switch_interactive = 0
f = plt.figure(num=None, dpi=700, facecolor='w', edgecolor='k')
m = Basemap(projection= 'lcc', lat_0 = meanlatd02,
            lon_0 = meanlond02,
            llcrnrlon=(-125.5), llcrnrlat=(57.2),
            urcrnrlon=(-105.4),urcrnrlat = (63),
            resolution = 'f')
m.drawcoastlines(linewidth=0.3)
m.drawcountries()
m.drawparallels(np.arange(-50,70,5))
m.drawmeridians(np.arange(-130,-100,5))

seismic_mod = cm.get_cmap('seismic',256)
newcmp = ListedColormap(seismic_mod(np.linspace(0.5,1,256)))
m.pcolormesh(lons_2015_d02,lats_2015_d02,avg_daily_CAPExP,latlon=True,cmap=newcmp)

if switch_interactive=="1":
    plt.show()
else:
    outpath = "/scratch/leuven/336/vsc33651/nu-wrf-dev/Plots"
    fname = "avg_CAPExP_pd_2015_d01"
    fname_long = os.path.join(outpath, fname+'.png')
    #plt.title('Average daily Lightning Flashes 2015 (June - August')
    plt.colorbar()
    plt.savefig(fname_long, dpi=f.dpi)
    plt.close()