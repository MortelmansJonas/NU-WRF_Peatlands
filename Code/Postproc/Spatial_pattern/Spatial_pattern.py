# ---------------------------------------------------------------------------------------------
# MODULES
# ---------------------------------------------------------------------------------------------
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import LogNorm, from_levels_and_colors, ListedColormap, LinearSegmentedColormap
from matplotlib import cm

# ---------------------------------------------------------------------------------------------
# LOAD DATA
# ---------------------------------------------------------------------------------------------
ds_d02 = Dataset('/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/wrfout_nc_files/domain2_daily_Thompson.nc', mode='r')
ds_d01 = Dataset('/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/wrfout_nc_files/domain1_daily_Thompson.nc', mode='r')
lats_d02 = ds_d02['lat'][:]
lons_d02 = ds_d02['lon'][:]
lats_d01 = ds_d01['lat'][:]
lons_d01 = ds_d01['lon'][:]

avg_daily_LPI_d02 = np.nanmean(ds_d02['LPI'][:], axis = 0)
avg_daily_LPI_d01 = np.nanmean(ds_d01['LPI'][:], axis = 0)
avg_daily_LTG3_d02 = np.nanmean(ds_d02['LTG3'][:], axis = 0)
avg_daily_LTG3_d01 = np.nanmean(ds_d01['LTG3'][:], axis = 0)
avg_daily_PR92W_d02 = np.nanmean(ds_d02['PR92W'][:], axis = 0)
avg_daily_PR92W_d01 = np.nanmean(ds_d01['PR92W'][:], axis = 0)
avg_daily_CAPExP_R_d01 = np.nanmean(ds_d01['CAPExP_R'][:], axis = 0)
avg_daily_CAPExP_CSI_d01 = np.nanmean(ds_d01['CAPExP_CSI'][:], axis = 0)
avg_daily_CAPExP_R_d02 = np.nanmean(ds_d02['CAPExP_R'][:], axis = 0)
avg_daily_CAPExP_CSI_d02 = np.nanmean(ds_d02['CAPExP_CSI'][:], axis = 0)
avg_daily_obs_d01 = np.nanmean(ds_d01['Obs'][:], axis = 0)
avg_daily_obs_d02 = np.nanmean(ds_d02['Obs'][:], axis = 0)

# ---------------------------------------------------------------------------------------------
# MAP DOMAIN 1
# ---------------------------------------------------------------------------------------------
print('maps')
# specify before all subplots
switch_interactive = 0

# SPECTRAL colorbar
spectral_mod = cm.get_cmap('Spectral',256)
newcolors = spectral_mod(np.linspace(0.5,1,5))
white = np.array([256/256, 256/256, 256/256, 1])
newcolors[0,:] = white
newcmp = ListedColormap(newcolors)

x1 = 0.0005
my_fontsize= 12

fig, axes = plt.subplots(3,2, figsize=(15.3, 8.27))
fig.suptitle('THOM - 9 km',fontsize=20)
m1 = Basemap(projection= 'lcc', lat_0 = lats_d01.mean(),
            lon_0 = lons_d01.mean(),
            llcrnrlon=(-123), llcrnrlat=(58),
            urcrnrlon=(-108),urcrnrlat = (62),
            resolution = 'f', ax=axes[0][0])
m1.drawcoastlines(linewidth=0.3)
m1.drawcountries()
m1.drawparallels(np.arange(-50,70,5))
m1.drawmeridians(np.arange(-130,-100,5))
m1.drawmapscale(-121, 58.7, lons_d01.mean(), lats_d01.mean(), 100)
m1.drawparallels(np.arange(-50,70,5), linewidth=0.5,labels = [True,False,False,False])
m1.drawmeridians(np.arange(-130,-100,5), linewidth=0.5,labels = [False,False,False, True])
map1= m1.pcolormesh(lons_d01,lats_d01,avg_daily_obs_d01,latlon=True,cmap=newcmp)
axes[0][0].set_title('CLDN\n', fontsize=16)
mstats = 'm = %.4f, s = %.4f' % (np.nanmean(avg_daily_obs_d01), np.nanstd(avg_daily_obs_d01))
plt.text(1.0,1.0, mstats,  horizontalalignment='right', verticalalignment='bottom', transform=axes[0][0].transAxes, fontsize=my_fontsize)
map1.set_clim(0,0.05)


m2 = Basemap(projection= 'lcc', lat_0 = lats_d01.mean(),
            lon_0 = lons_d01.mean(),
            llcrnrlon=(-123), llcrnrlat=(58),
            urcrnrlon=(-108),urcrnrlat = (62),
            resolution = 'f',ax=axes[2][0])
m2.drawcoastlines(linewidth=0.3)
m2.drawcountries()
m2.drawparallels(np.arange(-50,70,5))
m2.drawmeridians(np.arange(-130,-100,5))
m2.drawparallels(np.arange(-50,70,5), linewidth=0.5,labels = [True,False,False,False])
m2.drawmeridians(np.arange(-130,-100,5), linewidth=0.5,labels = [False,False,False, True])
map2 = m2.pcolormesh(lons_d01,lats_d01,avg_daily_CAPExP_R_d01,latlon=True,cmap=newcmp)
axes[2][0].set_title('CAPExP_R\n', fontsize=16)
mstats = 'm = %.4f, s = %.4f' % (np.nanmean(avg_daily_CAPExP_R_d01), np.nanstd(avg_daily_CAPExP_R_d01))
plt.text(1.0,1.0, mstats,  horizontalalignment='right', verticalalignment='bottom', transform=axes[2][0].transAxes, fontsize=my_fontsize)
map2.set_clim(0,0.05)

m3 = Basemap(projection= 'lcc', lat_0 = lats_d01.mean(),
            lon_0 = lons_d01.mean(),
            llcrnrlon=(-123), llcrnrlat=(58),
            urcrnrlon=(-108),urcrnrlat = (62),
            resolution = 'f',ax=axes[0][1])
m3.drawcoastlines(linewidth=0.3)
m3.drawcountries()
m3.drawparallels(np.arange(-50,70,5))
m3.drawmeridians(np.arange(-130,-100,5))
m3.drawparallels(np.arange(-50,70,5), linewidth=0.5,labels = [True,False,False,False])
m3.drawmeridians(np.arange(-130,-100,5), linewidth=0.5,labels = [False,False,False, True])
map3=m3.pcolormesh(lons_d01,lats_d01,avg_daily_LPI_d01,latlon=True,cmap=newcmp)
axes[0][1].set_title('LPI\n', fontsize=16)
mstats = 'm = %.4f, s = %.4f' % (np.nanmean(avg_daily_LPI_d01), np.nanstd(avg_daily_LPI_d01))
plt.text(1.0,1.0, mstats,  horizontalalignment='right', verticalalignment='bottom', transform=axes[0][1].transAxes, fontsize=my_fontsize)
map3.set_clim(0,0.05)

m4 = Basemap(projection= 'lcc', lat_0 = lats_d01.mean(),
            lon_0 = lons_d01.mean(),
            llcrnrlon=(-123), llcrnrlat=(58),
            urcrnrlon=(-108),urcrnrlat = (62),
            resolution = 'f',ax=axes[1][0])
m4.drawcoastlines(linewidth=0.3)
m4.drawcountries()
m4.drawparallels(np.arange(-50,70,5))
m4.drawmeridians(np.arange(-130,-100,5))
m4.drawparallels(np.arange(-50,70,5), linewidth=0.5,labels = [True,False,False,False])
m4.drawmeridians(np.arange(-130,-100,5), linewidth=0.5,labels = [False,False,False, True])
map4=m4.pcolormesh(lons_d01,lats_d01,avg_daily_LTG3_d01,latlon=True,cmap=newcmp)
axes[1][0].set_title('LT3\n', fontsize=16)
mstats = 'm = %.4f, s = %.4f' % (np.nanmean(avg_daily_LTG3_d01), np.nanstd(avg_daily_LTG3_d01))
plt.text(1.0,1.0, mstats,  horizontalalignment='right', verticalalignment='bottom', transform=axes[1][0].transAxes, fontsize=my_fontsize)
map4.set_clim(0,0.05)

m5 = Basemap(projection= 'lcc', lat_0 = lats_d01.mean(),
            lon_0 = lons_d01.mean(),
            llcrnrlon=(-123), llcrnrlat=(58),
            urcrnrlon=(-108),urcrnrlat = (62),
            resolution = 'f',ax=axes[2][1])
m5.drawcoastlines(linewidth=0.3)
m5.drawcountries()
m5.drawparallels(np.arange(-50,70,5))
m5.drawmeridians(np.arange(-130,-100,5))
m5.drawparallels(np.arange(-50,70,5), linewidth=0.5,labels = [True,False,False,False])
m5.drawmeridians(np.arange(-130,-100,5), linewidth=0.5,labels = [False,False,False, True])
map5 = m5.pcolormesh(lons_d01,lats_d01,avg_daily_CAPExP_CSI_d01,latlon=True,cmap=newcmp)
axes[2][1].set_title('CAPExP_CSI\n', fontsize=16)
mstats = 'm = %.4f, s = %.4f' % (np.nanmean(avg_daily_CAPExP_CSI_d01), np.nanstd(avg_daily_CAPExP_CSI_d01))
plt.text(1.0,1.0, mstats,  horizontalalignment='right', verticalalignment='bottom', transform=axes[2][1].transAxes, fontsize=my_fontsize)
map5.set_clim(0,0.05)


m6 = Basemap(projection= 'lcc', lat_0 = lats_d01.mean(),
            lon_0 = lons_d01.mean(),
            llcrnrlon=(-123), llcrnrlat=(58),
            urcrnrlon=(-108),urcrnrlat = (62),
            resolution = 'f', ax=axes[1][1])
m6.drawcoastlines(linewidth=0.3)
m6.drawcountries()
m6.drawparallels(np.arange(-50,70,5))
m6.drawmeridians(np.arange(-130,-100,5))
m6.drawparallels(np.arange(-50,70,5), linewidth=0.5,labels = [True,False,False,False])
m6.drawmeridians(np.arange(-130,-100,5), linewidth=0.5,labels = [False,False,False, True])
map6 =m6.pcolormesh(lons_d01,lats_d01,avg_daily_PR92W_d01,latlon=True,cmap=newcmp)
axes[1][1].set_title('PR92W\n', fontsize = 16)
mstats = 'm = %.4f, s = %.4f' % (np.nanmean(avg_daily_PR92W_d01), np.nanstd(avg_daily_PR92W_d01))
plt.text(1.0,1.0, mstats,  horizontalalignment='right', verticalalignment='bottom', transform=axes[1][1].transAxes, fontsize=my_fontsize)
map6.set_clim(0,0.05)

plt.subplots_adjust(left=0.173,
                    bottom=0.033,
                    right=0.808,
                    top=0.855,
                    wspace=0.0,
                    hspace=0.550)
plt.show()

## MAP DOMAIN 2
fig, axes = plt.subplots(3,2, figsize=(15.3, 8.27))
fig.suptitle('THOM - 3 km',fontsize=20)
m1 = Basemap(projection= 'lcc', lat_0 = lats_d02.mean(),
            lon_0 = lons_d02.mean(),
            llcrnrlon=(-123), llcrnrlat=(58),
            urcrnrlon=(-108),urcrnrlat = (62),
            resolution = 'f', ax=axes[0][0])
m1.drawcoastlines(linewidth=0.3)
m1.drawcountries()
m1.drawparallels(np.arange(-50,70,5))
m1.drawmeridians(np.arange(-130,-100,5))
m1.drawmapscale(-121, 58.7, lons_d02.mean(), lats_d02.mean(), 100)
m1.drawparallels(np.arange(-50,70,5), linewidth=0.5,labels = [True,False,False,False])
m1.drawmeridians(np.arange(-130,-100,5), linewidth=0.5,labels = [False,False,False, True])
map1= m1.pcolormesh(lons_d02,lats_d02,avg_daily_obs_d02,latlon=True,cmap=newcmp)
axes[0][0].set_title('CLDN\n', fontsize=16)
mstats = 'm = %.4f, s = %.4f' % (np.nanmean(avg_daily_obs_d02), np.nanstd(avg_daily_obs_d02))
plt.text(1.0,1.0, mstats,  horizontalalignment='right', verticalalignment='bottom', transform=axes[0][0].transAxes, fontsize=my_fontsize)

map1.set_clim(0,0.05)

m2 = Basemap(projection= 'lcc', lat_0 = lats_d02.mean(),
            lon_0 = lons_d02.mean(),
            llcrnrlon=(-123), llcrnrlat=(58),
            urcrnrlon=(-108),urcrnrlat = (62),
            resolution = 'f',ax=axes[2][0])
m2.drawcoastlines(linewidth=0.3)
m2.drawcountries()
m2.drawparallels(np.arange(-50,70,5))
m2.drawmeridians(np.arange(-130,-100,5))
m2.drawparallels(np.arange(-50,70,5), linewidth=0.5,labels = [True,False,False,False])
m2.drawmeridians(np.arange(-130,-100,5), linewidth=0.5,labels = [False,False,False, True])
map2 = m2.pcolormesh(lons_d02,lats_d02,avg_daily_CAPExP_R_d02,latlon=True,cmap=newcmp)
axes[2][0].set_title('CAPExP_R\n', fontsize=16)

mstats = 'm = %.4f, s = %.4f' % (np.nanmean(avg_daily_CAPExP_R_d02), np.nanstd(avg_daily_CAPExP_R_d02))
plt.text(1.0,1.0, mstats,  horizontalalignment='right', verticalalignment='bottom', transform=axes[2][0].transAxes, fontsize=my_fontsize)
map2.set_clim(0,0.05)

m3 = Basemap(projection= 'lcc', lat_0 = lats_d02.mean(),
            lon_0 = lons_d02.mean(),
            llcrnrlon=(-123), llcrnrlat=(58),
            urcrnrlon=(-108),urcrnrlat = (62),
            resolution = 'f',ax=axes[0][1])
m3.drawcoastlines(linewidth=0.3)
m3.drawcountries()
m3.drawparallels(np.arange(-50,70,5))
m3.drawmeridians(np.arange(-130,-100,5))
m3.drawparallels(np.arange(-50,70,5), linewidth=0.5,labels = [True,False,False,False])
m3.drawmeridians(np.arange(-130,-100,5), linewidth=0.5,labels = [False,False,False, True])
map3=m3.pcolormesh(lons_d02,lats_d02,avg_daily_LPI_d02,latlon=True,cmap=newcmp)
axes[0][1].set_title('LPI\n', fontsize=16)
mstats = 'm = %.4f, s = %.4f' % (np.nanmean(avg_daily_LPI_d02), np.nanstd(avg_daily_LPI_d02))
plt.text(1.0,1.0, mstats,  horizontalalignment='right', verticalalignment='bottom', transform=axes[0][1].transAxes, fontsize=my_fontsize)
map3.set_clim(0,0.05)

m4 = Basemap(projection= 'lcc', lat_0 = lats_d02.mean(),
            lon_0 = lons_d02.mean(),
            llcrnrlon=(-123), llcrnrlat=(58),
            urcrnrlon=(-108),urcrnrlat = (62),
            resolution = 'f',ax=axes[1][0])
m4.drawcoastlines(linewidth=0.3)
m4.drawcountries()
m4.drawparallels(np.arange(-50,70,5))
m4.drawmeridians(np.arange(-130,-100,5))
m4.drawparallels(np.arange(-50,70,5), linewidth=0.5,labels = [True,False,False,False])
m4.drawmeridians(np.arange(-130,-100,5), linewidth=0.5,labels = [False,False,False, True])
map4=m4.pcolormesh(lons_d02,lats_d02,avg_daily_LTG3_d02,latlon=True,cmap=newcmp)
axes[1][0].set_title('LT3\n', fontsize=16)
mstats = 'm = %.4f, s = %.4f' % (np.nanmean(avg_daily_LTG3_d02), np.nanstd(avg_daily_LTG3_d02))
plt.text(1.0,1.0, mstats,  horizontalalignment='right', verticalalignment='bottom', transform=axes[1][0].transAxes, fontsize=my_fontsize)
map4.set_clim(0,0.05)

m5 = Basemap(projection= 'lcc', lat_0 = lats_d02.mean(),
            lon_0 = lons_d02.mean(),
            llcrnrlon=(-123), llcrnrlat=(58),
            urcrnrlon=(-108),urcrnrlat = (62),
            resolution = 'f',ax=axes[2][1])
m5.drawcoastlines(linewidth=0.3)
m5.drawcountries()
m5.drawparallels(np.arange(-50,70,5))
m5.drawmeridians(np.arange(-130,-100,5))
m5.drawparallels(np.arange(-50,70,5), linewidth=0.5,labels = [True,False,False,False])
m5.drawmeridians(np.arange(-130,-100,5), linewidth=0.5,labels = [False,False,False, True])
map5 = m5.pcolormesh(lons_d02,lats_d02,avg_daily_CAPExP_CSI_d02,latlon=True,cmap=newcmp)
axes[2][1].set_title('CAPExP_CSI\n', fontsize=16)
mstats = 'm = %.4f, s = %.4f' % (np.nanmean(avg_daily_CAPExP_CSI_d02), np.nanstd(avg_daily_CAPExP_CSI_d02))
plt.text(1.0,1.0, mstats,  horizontalalignment='right', verticalalignment='bottom', transform=axes[2][1].transAxes, fontsize=my_fontsize)
map5.set_clim(0,0.05)

m6 = Basemap(projection= 'lcc', lat_0 = lats_d02.mean(),
            lon_0 = lons_d02.mean(),
            llcrnrlon=(-123), llcrnrlat=(58),
            urcrnrlon=(-108),urcrnrlat = (62),
            resolution = 'f', ax=axes[1][1])
m6.drawcoastlines(linewidth=0.3)
m6.drawcountries()
m6.drawparallels(np.arange(-50,70,5))
m6.drawmeridians(np.arange(-130,-100,5))
m6.drawparallels(np.arange(-50,70,5), linewidth=0.5,labels = [True,False,False,False])
m6.drawmeridians(np.arange(-130,-100,5), linewidth=0.5,labels = [False,False,False, True])
map6 =m6.pcolormesh(lons_d02,lats_d02,avg_daily_PR92W_d02,latlon=True,cmap=newcmp)
axes[1][1].set_title('PR92W\n', fontsize=16)
mstats = 'm = %.4f, s = %.4f' % (np.nanmean(avg_daily_PR92W_d02), np.nanstd(avg_daily_PR92W_d02))
plt.text(1.0,1.0, mstats,  horizontalalignment='right', verticalalignment='bottom', transform=axes[1][1].transAxes, fontsize=my_fontsize)
map6.set_clim(0,0.05)

plt.subplots_adjust(left=0.173,
                    bottom=0.033,
                    right=0.808,
                    top=0.855,
                    wspace=0.0,
                    hspace=0.550)
plt.show()
