# ---------------------------------------------------------------------------------------------
# MODULES
# ---------------------------------------------------------------------------------------------
import netCDF4 as nc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from netCDF4 import Dataset
import os
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import LogNorm, from_levels_and_colors, ListedColormap, LinearSegmentedColormap
from matplotlib import cm
import matplotlib.dates as mdates
from datetime import timedelta

ds_d02 = Dataset('/scratch/leuven/336/vsc33651/nu-wrf-dev/domain2_daily.nc', mode='r')
ds_d01 = Dataset('/scratch/leuven/336/vsc33651/nu-wrf-dev/domain1_daily.nc', mode='r')
# ds_obs_d01 = Dataset('/scratch/leuven/336/vsc33651/nu-wrf-dev/domain1_obs_summer_daily.nc')
ds_obs_d01 = Dataset('/scratch/leuven/336/vsc33651/nu-wrf-dev/wrfout_nc_files/Slave_lake_daily_0.1d.nc')
ds_obs_d02 = Dataset('/scratch/leuven/336/vsc33651/nu-wrf-dev/wrfout_nc_files/Slave_lake_daily_0.1d.nc')
lats_d02 = ds_d02['lat'][:]
lons_d02 = ds_d02['lon'][:]
lats_d01 = ds_d01['lat'][:]
lons_d01 = ds_d01['lon'][:]

time = pd.to_datetime(ds_obs_d01.variables['time'][:]) + timedelta(minutes=2) # Add 2 minutes again to get correct times
months = time.month
summer_2015 = np.where([(months > 5) & (months < 9)])
inds_obs_lat = np.unique(np.where((ds_obs_d01['lat'][:] > 58) & (ds_obs_d01['lat'][:] <63))[0])
inds_obs_lon = np.unique(np.where((ds_obs_d01['lon'][:] > -125) & (ds_obs_d01['lon'][:] <-108))[0])
CC_d01 = ds_obs_d01['Flashdensity_CC'][summer_2015[1],inds_obs_lat,inds_obs_lon] + ds_obs_d01['Flashdensity_CG'][summer_2015[1],inds_obs_lat,inds_obs_lon]
lats_obs_d01 = ds_obs_d01['lat'][inds_obs_lat]
lons_obs_d01 = ds_obs_d01['lon'][inds_obs_lon]

inds_obs_lat_d02 = np.unique(np.where((ds_obs_d02['lat'][:] > 58) & (ds_obs_d02['lat'][:] <63))[0])
inds_obs_lon_d02 = np.unique(np.where((ds_obs_d02['lon'][:] > -125) & (ds_obs_d02['lon'][:] <-108))[0])
CC_d02 = ds_obs_d02['Flashdensity_CC'][summer_2015[1],inds_obs_lat_d02,inds_obs_lon_d02] + ds_obs_d02['Flashdensity_CG'][summer_2015[1],inds_obs_lat_d02,inds_obs_lon_d02]
lats_obs_d02 = ds_obs_d02['lat'][inds_obs_lat_d02]
lons_obs_d02 = ds_obs_d02['lon'][inds_obs_lon_d02]
print(np.nanmax(ds_d01['PR92H'][:]))
print(np.nanmax(ds_d02['PR92H'][:]))
avg_daily_LPI_d02 = np.nanmean(ds_d02['LPI'][:], axis = 0)
avg_daily_LPI_d01 = np.nanmean(ds_d01['LPI'][:], axis = 0)
avg_daily_LTG3_d02 = np.nanmean(ds_d02['LTG3'][:], axis = 0)
avg_daily_LTG3_d01 = np.nanmean(ds_d01['LTG3'][:], axis = 0)
avg_daily_PR92_H_d02 = np.nanmean(ds_d02['PR92H'][:], axis = 0)
avg_daily_PR92_H_d01 = np.nanmean(ds_d01['PR92H'][:], axis = 0)
avg_daily_PR92_W_d02 = np.nanmean(ds_d02['PR92W'][:], axis = 0)
avg_daily_PR92_W_d01 = np.nanmean(ds_d01['PR92W'][:], axis = 0)
avg_daily_CAPExP_d01 = np.nanmean(ds_d01['CAPExP'][:], axis = 0)
avg_daily_obs_d01 = np.nanmean(CC_d01[:], axis = 0)
avg_daily_obs_d02 = np.nanmean(CC_d02[:], axis = 0)
print((avg_daily_PR92_H_d01))
print((avg_daily_PR92_H_d02))

switch_interactive = 0
# seismic_mod = cm.get_cmap('seismic',256)
# newcmp = ListedColormap(seismic_mod(np.linspace(0.5,1,256)))
#
# f = plt.figure(num=None, dpi=700, facecolor='w', edgecolor='k')
# m = Basemap(projection= 'lcc', lat_0 = lats_d01.mean(),
#             lon_0 = lons_d01.mean(),
#             llcrnrlon=(-125), llcrnrlat=(58),
#             urcrnrlon=(-108),urcrnrlat = (63),
#             resolution = 'f')
# m.drawcoastlines(linewidth=0.3)
# m.drawcountries()
# m.drawparallels(np.arange(-50,70,5))
# m.drawmeridians(np.arange(-130,-100,5))
# m.drawmapscale(-123, 58.7, lons_d01.mean(), lats_d01.mean(), 100)
# m.drawparallels(np.arange(-50,70,5), linewidth=0.5,labels = [True,False,False,False])
# m.drawmeridians(np.arange(-130,-100,5), linewidth=0.5,labels = [False,False,False, True])
# m.pcolormesh(lons_d01,lats_d01,avg_daily_PR92_H_d01,latlon=True,cmap=newcmp)
# plt.colorbar(label='Number of flashes (day$^{-1}$ (km)$^{-2}$)',extend = 'max',shrink = 0.6)
# plt.title('PR92H')
# plt.show()


# f = plt.figure(num=None, dpi=700, facecolor='w', edgecolor='k')
# m = Basemap(projection= 'lcc', lat_0 = lats_d02.mean(),
#             lon_0 = lons_d02.mean(),
#             llcrnrlon=(-125), llcrnrlat=(58),
#             urcrnrlon=(-108),urcrnrlat = (63),
#             resolution = 'f')
# m.drawcoastlines(linewidth=0.3)
# m.drawcountries()
# m.drawparallels(np.arange(-50,70,5))
# m.drawmeridians(np.arange(-130,-100,5))
# m.drawmapscale(-123, 58.7, lons_d02.mean(), lats_d02.mean(), 100)
# m.drawparallels(np.arange(-50,70,5), linewidth=0.5,labels = [True,False,False,False])
# m.drawmeridians(np.arange(-130,-100,5), linewidth=0.5,labels = [False,False,False, True])
# m.pcolormesh(lons_d02,lats_d02,avg_daily_PR92_H_d02,latlon=True,cmap=newcmp)
# plt.colorbar(label='Number of flashes (day$^{-1}$ (km)$^{-2}$)',extend = 'max',shrink = 0.6)
# plt.title('PR92H')
# plt.show()



# ---------------------------------------------------------------------------------------------
# MAP DOMAIN 1
# ---------------------------------------------------------------------------------------------
print('maps')
switch_interactive = 0
seismic_mod = cm.get_cmap('seismic',256)
newcmp = ListedColormap(seismic_mod(np.linspace(0.5,1,256)))

print('d01')
fig = plt.figure()
plt.title('Convection-parameterized (9 km)')
ax1 = plt.subplot2grid((3,2),(0,0))
ax2 = plt.subplot2grid((3,2),(0,1))
ax3 = plt.subplot2grid((3,2),(1,0))
ax4 = plt.subplot2grid((3,2),(1,1))
ax5 = plt.subplot2grid((3,2),(2,0))
ax6 = plt.subplot2grid((3,2),(2,1))

m1 = Basemap(projection= 'lcc', lat_0 = lats_obs_d01.mean(),
            lon_0 = lons_obs_d01.mean(),
            llcrnrlon=(-125), llcrnrlat=(58),
            urcrnrlon=(-108),urcrnrlat = (63),
            resolution = 'f',ax=ax1)
m1.drawcoastlines(linewidth=0.3)
m1.drawcountries()
m1.drawparallels(np.arange(-50,70,5))
m1.drawmeridians(np.arange(-130,-100,5))
m1.drawmapscale(-123, 58.7, lons_obs_d01.mean(), lats_obs_d01.mean(), 100)
m1.drawparallels(np.arange(-50,70,5), linewidth=0.5,labels = [True,False,False,False])
m1.drawmeridians(np.arange(-130,-100,5), linewidth=0.5,labels = [False,False,False, True])
x,y = np.meshgrid(lons_obs_d01,lats_obs_d01)
m1.pcolormesh(x,y,avg_daily_obs_d01,latlon=True,cmap=newcmp)
# m1.colorbar(label='Number of flashes (day$^{-1}$ km$^{-2}$)',extend = 'max',shrink = 0.6)
ax1.set_title('CLDN')

m2 = Basemap(projection= 'lcc', lat_0 = lats_d01.mean(),
            lon_0 = lons_d01.mean(),
            llcrnrlon=(-125), llcrnrlat=(58),
            urcrnrlon=(-108),urcrnrlat = (63),
            resolution = 'f',ax=ax2)
m2.drawcoastlines(linewidth=0.3)
m2.drawcountries()
m2.drawparallels(np.arange(-50,70,5))
m2.drawmeridians(np.arange(-130,-100,5))
m2.drawparallels(np.arange(-50,70,5), linewidth=0.5,labels = [True,False,False,False])
m2.drawmeridians(np.arange(-130,-100,5), linewidth=0.5,labels = [False,False,False, True])
m2.pcolormesh(lons_d01,lats_d01,avg_daily_CAPExP_d01,latlon=True,cmap=newcmp)
# m2.colorbar(label='Number of flashes (day$^{-1}$ km$^{-2}$)',extend = 'max',shrink = 0.6)
ax2.set_title('CAPExP')

m3 = Basemap(projection= 'lcc', lat_0 = lats_d01.mean(),
            lon_0 = lons_d01.mean(),
            llcrnrlon=(-125), llcrnrlat=(58),
            urcrnrlon=(-108),urcrnrlat = (63),
            resolution = 'f',ax=ax3)
m3.drawcoastlines(linewidth=0.3)
m3.drawcountries()
m3.drawparallels(np.arange(-50,70,5))
m3.drawmeridians(np.arange(-130,-100,5))
m3.drawparallels(np.arange(-50,70,5), linewidth=0.5,labels = [True,False,False,False])
m3.drawmeridians(np.arange(-130,-100,5), linewidth=0.5,labels = [False,False,False, True])
m3.pcolormesh(lons_d01,lats_d01,avg_daily_LPI_d01,latlon=True,cmap=newcmp)
# m3.colorbar(label='LPI (J kg$^{-1}$)',extend = 'max',shrink = 0.6)
ax3.set_title('LPI')

m4 = Basemap(projection= 'lcc', lat_0 = lats_d01.mean(),
            lon_0 = lons_d01.mean(),
            llcrnrlon=(-125), llcrnrlat=(58),
            urcrnrlon=(-108),urcrnrlat = (63),
            resolution = 'f',ax=ax4)
m4.drawcoastlines(linewidth=0.3)
m4.drawcountries()
m4.drawparallels(np.arange(-50,70,5))
m4.drawmeridians(np.arange(-130,-100,5))
m4.drawparallels(np.arange(-50,70,5), linewidth=0.5,labels = [True,False,False,False])
m4.drawmeridians(np.arange(-130,-100,5), linewidth=0.5,labels = [False,False,False, True])
m4.pcolormesh(lons_d01,lats_d01,avg_daily_LTG3_d01,latlon=True,cmap=newcmp)
# m4.colorbar(label='Number of flashes (day$^{-1}$ km$^{-2}$)',extend = 'max',shrink = 0.6)
ax4.set_title('LT3')

m5 = Basemap(projection= 'lcc', lat_0 = lats_d01.mean(),
            lon_0 = lons_d01.mean(),
            llcrnrlon=(-125), llcrnrlat=(58),
            urcrnrlon=(-108),urcrnrlat = (63),
            resolution = 'f')
m5.drawcoastlines(linewidth=0.3,ax=ax5)
m5.drawcountries()
m5.drawparallels(np.arange(-50,70,5))
m5.drawmeridians(np.arange(-130,-100,5))
m5.drawparallels(np.arange(-50,70,5), linewidth=0.5,labels = [True,False,False,False])
m5.drawmeridians(np.arange(-130,-100,5), linewidth=0.5,labels = [False,False,False, True])
m5.pcolormesh(lons_d01,lats_d01,avg_daily_PR92_H_d01,latlon=True,cmap=newcmp)
# m5.colorbar(label='Number of flashes (day$^{-1}$ km$^{-2}$)',extend = 'max',shrink = 0.6)
ax5.set_title('PR92H')

m6 = Basemap(projection= 'lcc', lat_0 = lats_d01.mean(),
            lon_0 = lons_d01.mean(),
            llcrnrlon=(-125), llcrnrlat=(58),
            urcrnrlon=(-108),urcrnrlat = (63),
            resolution = 'f', ax=ax6)
m6.drawcoastlines(linewidth=0.3)
m6.drawcountries()
m6.drawparallels(np.arange(-50,70,5))
m6.drawmeridians(np.arange(-130,-100,5))
m6.drawparallels(np.arange(-50,70,5), linewidth=0.5,labels = [True,False,False,False])
m6.drawmeridians(np.arange(-130,-100,5), linewidth=0.5,labels = [False,False,False, True])
m6.pcolormesh(lons_d01,lats_d01,avg_daily_PR92_W_d01,latlon=True,cmap=newcmp)
# m6.colorbar(label='Number of flashes (day$^{-1}$ km$^{-2}$)',extend = 'max',shrink = 0.6)
ax6.set_title('PR92W')

plt.subplots_adjust(top=0.85)
fig.tight_layout()
plt.show()
plt.savefig("/scratch/leuven/336/vsc33651/nu-wrf-dev/Plots/ALL_d01.png")

# ---------------------------------------------------------------------------------------------
# MAP DOMAIN 2
# ---------------------------------------------------------------------------------------------
print('d02')
fig = plt.figure()
plt.title('Convection-permitting (3 km)')
plt.subplots_adjust(top=0.85)
ax1 = plt.subplot2grid((3,2),(0,0))
ax2 = plt.subplot2grid((3,2),(0,1))
ax3 = plt.subplot2grid((3,2),(1,0))
ax4 = plt.subplot2grid((3,2),(1,1))
ax5 = plt.subplot2grid((3,2),(2,0))
ax6 = plt.subplot2grid((3,2),(2,1))

m1 = Basemap(projection= 'lcc', lat_0 = lats_obs_d02.mean(),
            lon_0 = lons_obs_d02.mean(),
            llcrnrlon=(-125), llcrnrlat=(58),
            urcrnrlon=(-108),urcrnrlat = (63),
            resolution = 'f',ax=ax1)
m1.drawcoastlines(linewidth=0.3)
m1.drawcountries()
m1.drawparallels(np.arange(-50,70,5))
m1.drawmeridians(np.arange(-130,-100,5))
m1.drawmapscale(-123, 58.7, lons_obs_d02.mean(), lats_obs_d02.mean(), 100)
m1.drawparallels(np.arange(-50,70,5), linewidth=0.5,labels = [True,False,False,False])
m1.drawmeridians(np.arange(-130,-100,5), linewidth=0.5,labels = [False,False,False, True])
x,y = np.meshgrid(lons_obs_d02,lats_obs_d02)
m1.pcolormesh(x,y,avg_daily_obs_d02,latlon=True,cmap=newcmp)
# m1.colorbar(label='Number of flashes (day$^{-1}$ km$^{-2}$)',extend = 'max',shrink = 0.6)
ax1.set_title('CLDN')

m3 = Basemap(projection= 'lcc', lat_0 = lats_d02.mean(),
            lon_0 = lons_d02.mean(),
            llcrnrlon=(-125), llcrnrlat=(58),
            urcrnrlon=(-108),urcrnrlat = (63),
            resolution = 'f',ax=ax3)
m3.drawcoastlines(linewidth=0.3)
m3.drawcountries()
m3.drawparallels(np.arange(-50,70,5))
m3.drawmeridians(np.arange(-130,-100,5))
m3.drawparallels(np.arange(-50,70,5), linewidth=0.5,labels = [True,False,False,False])
m3.drawmeridians(np.arange(-130,-100,5), linewidth=0.5,labels = [False,False,False, True])
m3.pcolormesh(lons_d02,lats_d02,avg_daily_LPI_d02,latlon=True,cmap=newcmp)
# m3.colorbar(label='LPI (J kg$^{-1}$)',extend = 'max',shrink = 0.6)
ax3.set_title('LPI')

m4 = Basemap(projection= 'lcc', lat_0 = lats_d02.mean(),
            lon_0 = lons_d02.mean(),
            llcrnrlon=(-125), llcrnrlat=(58),
            urcrnrlon=(-108),urcrnrlat = (63),
            resolution = 'f',ax=ax4)
m4.drawcoastlines(linewidth=0.3)
m4.drawcountries()
m4.drawparallels(np.arange(-50,70,5))
m4.drawmeridians(np.arange(-130,-100,5))
m4.drawparallels(np.arange(-50,70,5), linewidth=0.5,labels = [True,False,False,False])
m4.drawmeridians(np.arange(-130,-100,5), linewidth=0.5,labels = [False,False,False, True])
m4.pcolormesh(lons_d02,lats_d02,avg_daily_LTG3_d02,latlon=True,cmap=newcmp)
# m4.colorbar(label='Number of flashes (day$^{-1}$ km$^{-2}$)',extend = 'max',shrink = 0.6)
ax4.set_title('LT3')

m5 = Basemap(projection= 'lcc', lat_0 = lats_d02.mean(),
            lon_0 = lons_d02.mean(),
            llcrnrlon=(-125), llcrnrlat=(58),
            urcrnrlon=(-108),urcrnrlat = (63),
            resolution = 'f')
m5.drawcoastlines(linewidth=0.3,ax=ax5)
m5.drawcountries()
m5.drawparallels(np.arange(-50,70,5))
m5.drawmeridians(np.arange(-130,-100,5))
m5.drawparallels(np.arange(-50,70,5), linewidth=0.5,labels = [True,False,False,False])
m5.drawmeridians(np.arange(-130,-100,5), linewidth=0.5,labels = [False,False,False, True])
m5.pcolormesh(lons_d02,lats_d02,avg_daily_PR92_H_d02,latlon=True,cmap=newcmp)
# m5.colorbar(label='Number of flashes (day$^{-1}$ km$^{-2}$)',extend = 'max',shrink = 0.6)
ax5.set_title('PR92H')

m6 = Basemap(projection= 'lcc', lat_0 = lats_d02.mean(),
            lon_0 = lons_d02.mean(),
            llcrnrlon=(-125), llcrnrlat=(58),
            urcrnrlon=(-108),urcrnrlat = (63),
            resolution = 'f', ax=ax6)
m6.drawcoastlines(linewidth=0.3)
m6.drawcountries()
m6.drawparallels(np.arange(-50,70,5))
m6.drawmeridians(np.arange(-130,-100,5))
m6.drawparallels(np.arange(-50,70,5), linewidth=0.5,labels = [True,False,False,False])
m6.drawmeridians(np.arange(-130,-100,5), linewidth=0.5,labels = [False,False,False, True])
m6.pcolormesh(lons_d02,lats_d02,avg_daily_PR92_W_d02,latlon=True,cmap=newcmp)
# m6.colorbar(label='Number of flashes (day$^{-1}$ km$^{-2}$)',extend = 'max',shrink = 0.6)
ax6.set_title('PR92W')

plt.subplots_adjust(top=0.85)
fig.tight_layout(h_pad = 2)
plt.show()
# plt.savefig("/scratch/leuven/336/vsc33651/nu-wrf-dev/Plots/ALL_d02.png")
