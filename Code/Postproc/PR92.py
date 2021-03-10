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
# FROM THE WRF FILES
# D02
ds_wrf_d02 = Dataset('/scratch/leuven/336/vsc33651/nu-wrf-dev/wrfout_nc_files/wrfout_d02_2015.nc', 'r')
lats_d02 = ds_wrf_d02.variables['lat'][:]
lons_d02 = ds_wrf_d02.variables['lon'][:]
time_2015= ds_wrf_d02.variables['time'][:]
CTOP2D_d02 = ds_wrf_d02['CTOP2D'][:]
COD2D_d02 = ds_wrf_d02['COD2D'][:]
W_UP_MAX_d02 = ds_wrf_d02['W_UP_MAX'][:]

# D01
ds_wrf_d01 = Dataset('/scratch/leuven/336/vsc33651/nu-wrf-dev/wrfout_nc_files/wrfout_d01_2015.nc', 'r')
lats_d01 = ds_wrf_d01.variables['lat'][:]
lons_d01 = ds_wrf_d01.variables['lon'][:]
latinds_d02 = np.where((lats_d01 >= np.min(lats_d02)) & (lats_d01 <= np.max(lats_d02)))[0]
loninds_d02 = np.where((lons_d01 >= np.min(lons_d02)) & (lons_d01 <= np.max(lons_d02)))[0]
CTOP2D_d01 = ds_wrf_d01['CTOP2D'][:, latinds_d02, loninds_d02]
COD2D_d01 = ds_wrf_d01['COD2D'][:, latinds_d02, loninds_d02]
W_UP_MAX_d01 = ds_wrf_d01['W_UP_MAX'][:, latinds_d02, loninds_d02]

# FROM THE UPP FILE
# D02
ds_upp_d02 = Dataset('/staging/leuven/stg_00024/OUTPUT/michelb/nu-wrf-dev/Great_Slave_Lake/2015/WRFPRS_d02_2015.nc', mode='r')
Cdim_d02 = ds_upp_d02['Cdim'][:]
Ctop_d02 = ds_upp_d02['Ctop'][:]
Cbot_d02 = ds_upp_d02['Cbot'][:]

# D01
ds_upp_d01 = Dataset('/staging/leuven/stg_00024/OUTPUT/michelb/nu-wrf-dev/Great_Slave_Lake/2015/WRFPRS_d01_2015.nc', mode='r')
Cdim_d01 = ds_upp_d01['Cdim'][:, latinds_d02, loninds_d02]
Ctop_d01 = ds_upp_d01['Ctop'][:, latinds_d02, loninds_d02]
Cbot_d01 = ds_upp_d01['Cbot'][:, latinds_d02, loninds_d02]

# ---------------------------------------------------------------------------------------------
# DETERMINE CELLS WITH CONVECTIVE CLOUDS
# ---------------------------------------------------------------------------------------------
# D02
CTOP_d02 = np.where(CTOP2D_d02 <440,1,0) # CONVECTIVE CLOUD: TOP PRESSURE < 440 mbar
COD_d02 = np.where(COD2D_d02>23, 1,0) # CONVECTIVE CLOUD: OPTICAL DENSITY > 23
Conv_cloud_d02 = np.multiply(CTOP_d02,COD_d02) # 1 = CONV. CLOUD, 0 = NO CONV. CLOUD

# D01
CTOP_d01 = np.where(CTOP2D_d01 < 440,1,0) # CONVECTIVE CLOUD: TOP PRESSURE < 440 mbar
COD_d01 = np.where(COD2D_d01 > 23, 1,0) # CONVECTIVE CLOUD: OPTICAL DENSITY > 23
Conv_cloud_d01 = np.multiply(CTOP_d01,COD_d01) # 1 = CONV. CLOUD, 0 = NO CONV. CLOUD

# ---------------------------------------------------------------------------------------------
# CALCULATE PR92 FOR CELLS WITH CONVECTIVE CLOUDS FOR CLOUD TOP HEIGHT
# ---------------------------------------------------------------------------------------------
# D02
power_d02 = np.power(Cdim_d02, 4.9)
PR92_per_min_d02 = np.multiply(0.0000344, power_d02)
PR92_per_hour_d02 = np.multiply(60, PR92_per_min_d02)
PR92_conv_d02 = np.multiply(PR92_per_hour_d02, Conv_cloud_d02)

# D01
power_d01 = np.power(Cdim_d01, 4.9)
PR92_per_min_d01 = np.multiply(0.0000344, power_d01)
PR92_per_hour_d01 = np.multiply(60, PR92_per_min_d01)
PR92_conv_d01 = np.multiply(PR92_per_hour_d01, Conv_cloud_d01)

# ---------------------------------------------------------------------------------------------
# CALCULATE PR92 FOR CELLS WITH CONVECTIVE CLOUDS FOR MAX VERTICAL UPDRAFT VELOCITY
# ---------------------------------------------------------------------------------------------
power_wu = np.power(W_UP_MAX, 5.45)
PR92_WU_per_min = np.multiply(0.000005, power_wu)
PR92_WU_ph = np.multiply(60,PR92_WU_per_min)

# ---------------------------------------------------------------------------------------------
# CALCULATE AVERAGE DAILY PR92
# ---------------------------------------------------------------------------------------------
days = np.divide(2208,24).astype(int)
daily_PR92 = np.zeros((days,216,363))
daily_PR92_WU = np.zeros((days,216,363))
for i in range(0,days):
    j = (i+1)*24
    k = i*24
    daily_PR92[i,:,:] = np.sum(PR92_conv[k:j,:,:], axis = 0)
    daily_PR92_WU[i,:,:] = np.sum(PR92_WU_ph[k:j,:,:], axis = 0)
avg_daily_PR92 = np.mean(daily_PR92, axis = 0)
avg_daily_PR92_WU = np.mean(daily_PR92_WU, axis = 0)

# ---------------------------------------------------------------------------------------------
# MAPS
# ---------------------------------------------------------------------------------------------
switch_interactive = 0
seismic_mod = cm.get_cmap('seismic',256)
newcmp = ListedColormap(seismic_mod(np.linspace(0.5,1,256)))
# MAP
# D02
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
m.drawmapscale(np.amin(lons_d02)+2.5, np.amin(lats_d02)+0.6, lons_d02.mean(), lats_d02.mean(), 100)
m.drawparallels(np.arange(-50,70,5), linewidth=0.5,labels = [True,False,False,False])
m.drawmeridians(np.arange(-130,-100,5), linewidth=0.5,labels = [False,False,False, True])
m.pcolormesh(lons_d02,lats_d02,avg_daily_PR92,latlon=True,cmap=newcmp)

if switch_interactive=="1":
    plt.show()
else:
    outpath = "/scratch/leuven/336/vsc33651/nu-wrf-dev/Plots"
    fname = "PR92_2015_d02"
    fname_long = os.path.join(outpath, fname+'.png')
    plt.title('Convection permitting (3-km)',size=12)
    plt.colorbar(label='Number of flashes (day$^{-1}$ (3km)$^{-2}$)',extend = 'max',shrink = 0.6)
    plt.savefig(fname_long, dpi=f.dpi)
    plt.close()

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
m.drawmapscale(np.amin(lons_d02)+2.5, np.amin(lats_d02)+0.6, lons_d02.mean(), lats_d02.mean(), 100)
m.drawparallels(np.arange(-50,70,5), linewidth=0.5,labels = [True,False,False,False])
m.drawmeridians(np.arange(-130,-100,5), linewidth=0.5,labels = [False,False,False, True])
m.pcolormesh(lons_d02,lats_d02,avg_daily_PR92_WU,latlon=True,cmap=newcmp)

if switch_interactive=="1":
    plt.show()
else:
    outpath = "/scratch/leuven/336/vsc33651/nu-wrf-dev/Plots"
    fname = "PR92_WU_2015_d02"
    fname_long = os.path.join(outpath, fname+'.png')
    plt.title('Convection permitting (3-km)',size=12)
    plt.colorbar(label='Number of flashes (day$^{-1}$ (3km)$^{-2}$)',extend = 'max',shrink = 0.6)
    plt.clim(0,100)
    plt.savefig(fname_long, dpi=f.dpi)
    plt.close()

# D01
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
m.drawmapscale(np.amin(lons_d02)+2.5, np.amin(lats_d02)+0.6, lons_d02.mean(), lats_d02.mean(), 100)
m.drawparallels(np.arange(-50,70,5), linewidth=0.5,labels = [True,False,False,False])
m.drawmeridians(np.arange(-130,-100,5), linewidth=0.5,labels = [False,False,False, True])
m.pcolormesh(lons_d02,lats_d02,avg_daily_PR92_d01,latlon=True,cmap=newcmp)

if switch_interactive=="1":
    plt.show()
else:
    outpath = "/scratch/leuven/336/vsc33651/nu-wrf-dev/Plots"
    fname = "PR92_2015_d01"
    fname_long = os.path.join(outpath, fname+'.png')
    plt.title('Convection parameterized (9-km)',size=12)
    plt.colorbar(label='Number of flashes (day$^{-1}$ (9km)$^{-2}$)',extend = 'max',shrink = 0.6)
    plt.savefig(fname_long, dpi=f.dpi)
    plt.close()

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
m.drawmapscale(np.amin(lons_d02)+2.5, np.amin(lats_d02)+0.6, lons_d02.mean(), lats_d02.mean(), 100)
m.drawparallels(np.arange(-50,70,5), linewidth=0.5,labels = [True,False,False,False])
m.drawmeridians(np.arange(-130,-100,5), linewidth=0.5,labels = [False,False,False, True])
m.pcolormesh(lons_d02,lats_d02,avg_daily_PR92_WU_d01,latlon=True,cmap=newcmp)

if switch_interactive=="1":
    plt.show()
else:
    outpath = "/scratch/leuven/336/vsc33651/nu-wrf-dev/Plots"
    fname = "PR92_WU_2015_d01"
    fname_long = os.path.join(outpath, fname+'.png')
    plt.title('Convection parameterized (9-km)',size=12)
    plt.colorbar(label='Number of flashes (day$^{-1}$ (9km)$^{-2}$)',extend = 'max',shrink = 0.6)
    plt.clim(0,100)
    plt.savefig(fname_long, dpi=f.dpi)
    plt.close()