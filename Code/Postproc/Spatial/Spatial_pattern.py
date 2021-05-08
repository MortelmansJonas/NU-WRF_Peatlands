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
ds_d02 = Dataset('/scratch/leuven/336/vsc33651/nu-wrf-dev/domain2_daily.nc', mode='r')
ds_d01 = Dataset('/scratch/leuven/336/vsc33651/nu-wrf-dev/domain1_daily.nc', mode='r')
lats_d02 = ds_d02['lat'][:]
lons_d02 = ds_d02['lon'][:]
lats_d01 = ds_d01['lat'][:]
lons_d01 = ds_d01['lon'][:]

avg_daily_LPI_d02 = np.nanmean(ds_d02['LPI'][:], axis = 0)
avg_daily_LPI_d01 = np.nanmean(ds_d01['LPI'][:], axis = 0)
avg_daily_LTG3_d02 = np.nanmean(ds_d02['LTG3'][:], axis = 0)
avg_daily_LTG3_d01 = np.nanmean(ds_d01['LTG3'][:], axis = 0)
avg_daily_PR92_H_d02 = np.nanmean(ds_d02['PR92H'][:], axis = 0)
avg_daily_PR92_H_d01 = np.nanmean(ds_d01['PR92H'][:], axis = 0)
avg_daily_PR92_W_d02 = np.nanmean(ds_d02['PR92W'][:], axis = 0)
avg_daily_PR92_W_d01 = np.nanmean(ds_d01['PR92W'][:], axis = 0)
avg_daily_CAPExP_d01 = np.nanmean(ds_d01['CAPExP'][:], axis = 0)
avg_daily_obs_d01 = np.nanmean(ds_d01['Obs'][:], axis = 0)
avg_daily_obs_d02 = np.nanmean(ds_d02['Obs'][:], axis = 0)

# ---------------------------------------------------------------------------------------------
# MAP DOMAIN 1
# ---------------------------------------------------------------------------------------------
print('maps')
# specify before all subplots
switch_interactive = 0
seismic_mod = cm.get_cmap('seismic',256)
newcmp = ListedColormap(seismic_mod(np.linspace(0.5,1,256)))
x1 = 0.00005
my_fontsize= 10

fig, axes = plt.subplots(3,2, figsize=(5.85, 8.27))
fig.suptitle('Convection-parameterized (9 km)',fontsize=16)
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
axes[0][0].set_title('CLDN\n')
axes[0][0].annotate('(a) \n ', xy=(x1,  axes[0][0].get_ylim()[1]),annotation_clip=False)
mstats = 'm = %.4f, s = %.4f' % (np.nanmean(avg_daily_obs_d01), np.nanstd(avg_daily_obs_d01))
plt.text(1.0,1.0, mstats,  horizontalalignment='right', verticalalignment='bottom', transform=axes[0][0].transAxes, fontsize=my_fontsize)
fig.colorbar(map1,ax=axes[0][0], label='Flash density \n (# day$^{-1}$ km$^{-2}$)',extend = 'max', shrink = 0.6)

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
map2 = m2.pcolormesh(lons_d01,lats_d01,avg_daily_CAPExP_d01,latlon=True,cmap=newcmp)
axes[2][0].set_title('CAPExP\n')
fig.colorbar(map2, ax=axes[2][0],label='Flash density \n (# day$^{-1}$ km$^{-2}$)',extend = 'max', shrink = 0.6)
axes[2][0].annotate('(e) \n ', xy=(x1,  axes[2][0].get_ylim()[1]),annotation_clip=False)
mstats = 'm = %.4f, s = %.4f' % (np.nanmean(avg_daily_CAPExP_d01), np.nanstd(avg_daily_CAPExP_d01))
plt.text(1.0,1.0, mstats,  horizontalalignment='right', verticalalignment='bottom', transform=axes[2][0].transAxes, fontsize=my_fontsize)

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
axes[0][1].set_title('LPI\n')
plt.colorbar(map3,ax=axes[0][1],label='LPI (J kg$^{-1}$)',extend = 'max', shrink = 0.6)
axes[0][1].annotate('(b) \n ', xy=(x1,  axes[0][1].get_ylim()[1]),annotation_clip=False)
mstats = 'm = %.4f, s = %.4f' % (np.nanmean(avg_daily_LPI_d01), np.nanstd(avg_daily_LPI_d01))
plt.text(1.0,1.0, mstats,  horizontalalignment='right', verticalalignment='bottom', transform=axes[0][1].transAxes, fontsize=my_fontsize)

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
axes[1][0].set_title('LT3\n')
plt.colorbar(map4,ax=axes[1][0],label='Flash density \n (# day$^{-1}$ km$^{-2}$)',extend = 'max', shrink = 0.6)
axes[1][0].annotate('(c) \n ', xy=(x1,  axes[1][0].get_ylim()[1]),annotation_clip=False)
mstats = 'm = %.4f, s = %.4f' % (np.nanmean(avg_daily_LTG3_d01), np.nanstd(avg_daily_LTG3_d01))
plt.text(1.0,1.0, mstats,  horizontalalignment='right', verticalalignment='bottom', transform=axes[1][0].transAxes, fontsize=my_fontsize)

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
map6 =m6.pcolormesh(lons_d01,lats_d01,avg_daily_PR92_W_d01,latlon=True,cmap=newcmp)
axes[1][1].set_title('PR92W\n')
plt.colorbar(map6,ax=axes[1][1], label='Flash density \n (# day$^{-1}$ km$^{-2}$)',extend = 'max', shrink = 0.6)
axes[1][1].annotate('(d) \n ', xy=(x1,  axes[1][1].get_ylim()[1]),annotation_clip=False)
mstats = 'm = %.4f, s = %.4f' % (np.nanmean(avg_daily_PR92_W_d01), np.nanstd(avg_daily_PR92_W_d01))
plt.text(1.0,1.0, mstats,  horizontalalignment='right', verticalalignment='bottom', transform=axes[1][1].transAxes, fontsize=my_fontsize)

plt.subplots_adjust(left=0.033,
                    bottom=0.033,
                    right=0.958,
                    top=0.881,
                    wspace=0.2,
                    hspace=0.397)
plt.show()
# ---------------------------------------------------------------------------------------------
# MAP DOMAIN 2
# ---------------------------------------------------------------------------------------------
print('d02')
fig = plt.figure(figsize=(5.85, 8.27))
fig.suptitle('Convection-permitting (3 km)', fontsize=16)
ax1 = plt.subplot2grid((2,2),(0,0))
ax3 = plt.subplot2grid((2,2),(0,1))
ax4 = plt.subplot2grid((2,2),(1,0))
ax6 = plt.subplot2grid((2,2),(1,1))

m1 = Basemap(projection= 'lcc', lat_0 = lats_d02.mean(),
            lon_0 = lons_d02.mean(),
            llcrnrlon=(-123), llcrnrlat=(58),
            urcrnrlon=(-108),urcrnrlat = (62),
            resolution = 'f',ax=ax1)
m1.drawcoastlines(linewidth=0.3)
m1.drawcountries()
m1.drawparallels(np.arange(-50,70,5))
m1.drawmeridians(np.arange(-130,-100,5))
m1.drawmapscale(-121, 58.7, lons_d02.mean(), lats_d02.mean(), 100)
m1.drawparallels(np.arange(-50,70,5), linewidth=0.5,labels = [True,False,False,False])
m1.drawmeridians(np.arange(-130,-100,5), linewidth=0.5,labels = [False,False,False, True])
map1= m1.pcolormesh(lons_d02,lats_d02,avg_daily_obs_d02,latlon=True,cmap=newcmp)
ax1.set_title('CLDN\n')
cbar = fig.colorbar(map1,ax=ax1, label='Flash density \n (# day$^{-1}$ km$^{-2}$)',extend = 'max', shrink = 0.6)
ax1.annotate('(a) \n ', xy=(x1,  ax1.get_ylim()[1]),annotation_clip=False)
mstats = 'm = %.4f, s = %.4f' % (np.nanmean(avg_daily_obs_d02), np.nanstd(avg_daily_obs_d02))
plt.text(1.0,1.0, mstats,  horizontalalignment='right', verticalalignment='bottom', transform=ax1.transAxes, fontsize=my_fontsize)

m3 = Basemap(projection= 'lcc', lat_0 = lats_d02.mean(),
            lon_0 = lons_d02.mean(),
            llcrnrlon=(-123), llcrnrlat=(58),
            urcrnrlon=(-108),urcrnrlat = (62),
            resolution = 'f',ax=ax3)
m3.drawcoastlines(linewidth=0.3)
m3.drawcountries()
m3.drawparallels(np.arange(-50,70,5))
m3.drawmeridians(np.arange(-130,-100,5))
m3.drawparallels(np.arange(-50,70,5), linewidth=0.5,labels = [True,False,False,False])
m3.drawmeridians(np.arange(-130,-100,5), linewidth=0.5,labels = [False,False,False, True])
map3=m3.pcolormesh(lons_d02,lats_d02,avg_daily_LPI_d02,latlon=True,cmap=newcmp)
ax3.set_title('LPI\n')
plt.colorbar(map3,ax=ax3,label='LPI (J kg$^{-1}$)',extend = 'max', shrink = 0.6)
ax3.annotate('(b) \n ', xy=(x1,  ax3.get_ylim()[1]),annotation_clip=False)
mstats = 'm = %.4f, s = %.4f' % (np.nanmean(avg_daily_LPI_d02), np.nanstd(avg_daily_LPI_d02))
plt.text(1.0,1.0, mstats,  horizontalalignment='right', verticalalignment='bottom', transform=ax3.transAxes, fontsize=my_fontsize)

m4 = Basemap(projection= 'lcc', lat_0 = lats_d02.mean(),
            lon_0 = lons_d02.mean(),
            llcrnrlon=(-123), llcrnrlat=(58),
            urcrnrlon=(-108),urcrnrlat = (62),
            resolution = 'f',ax=ax4)
m4.drawcoastlines(linewidth=0.3)
m4.drawcountries()
m4.drawparallels(np.arange(-50,70,5))
m4.drawmeridians(np.arange(-130,-100,5))
m4.drawparallels(np.arange(-50,70,5), linewidth=0.5,labels = [True,False,False,False])
m4.drawmeridians(np.arange(-130,-100,5), linewidth=0.5,labels = [False,False,False, True])
map4=m4.pcolormesh(lons_d02,lats_d02,avg_daily_LTG3_d02,latlon=True,cmap=newcmp)
ax4.set_title('LT3\n')
plt.colorbar(map4,ax=ax4,label='Flash density \n (# day$^{-1}$ km$^{-2}$)',extend = 'max', shrink = 0.6)
ax4.annotate('(c) \n ', xy=(x1,  ax4.get_ylim()[1]),annotation_clip=False)
mstats = 'm = %.4f, s = %.4f' % (np.nanmean(avg_daily_LTG3_d02), np.nanstd(avg_daily_LTG3_d02))
plt.text(1.0,1.0, mstats,  horizontalalignment='right', verticalalignment='bottom', transform=ax4.transAxes, fontsize=my_fontsize)

m6 = Basemap(projection= 'lcc', lat_0 = lats_d02.mean(),
            lon_0 = lons_d02.mean(),
            llcrnrlon=(-123), llcrnrlat=(58),
            urcrnrlon=(-108),urcrnrlat = (62),
            resolution = 'f', ax=ax6)
m6.drawcoastlines(linewidth=0.3)
m6.drawcountries()
m6.drawparallels(np.arange(-50,70,5))
m6.drawmeridians(np.arange(-130,-100,5))
m6.drawparallels(np.arange(-50,70,5), linewidth=0.5,labels = [True,False,False,False])
m6.drawmeridians(np.arange(-130,-100,5), linewidth=0.5,labels = [False,False,False, True])
map6 =m6.pcolormesh(lons_d02,lats_d02,avg_daily_PR92_W_d02,latlon=True,cmap=newcmp)
ax6.set_title('PR92W\n')
plt.colorbar(map6,ax=ax6, label='Flash density \n (# day$^{-1}$ km$^{-2}$)',extend = 'max', shrink = 0.6)
ax6.annotate('(d) \n ', xy=(x1,  ax6.get_ylim()[1]),annotation_clip=False)
mstats = 'm = %.4f, s = %.4f' % (np.nanmean(avg_daily_PR92_W_d02), np.nanstd(avg_daily_PR92_W_d02))
plt.text(1.0,1.0, mstats,  horizontalalignment='right', verticalalignment='bottom', transform=ax6.transAxes, fontsize=my_fontsize)

plt.subplots_adjust(left=0.125,
                    bottom=0.1,
                    right=0.9,
                    top=0.9,
                    wspace=0.2,
                    hspace=0.35)
plt.show()
