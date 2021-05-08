 ---------------------------------------------------------------------------------------------
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
ds_d01 = Dataset('/scratch/leuven/336/vsc33651/nu-wrf-dev/domain1_M2.nc', 'r')
ds_d02 = Dataset('/scratch/leuven/336/vsc33651/nu-wrf-dev/domain2_M2.nc', 'r')
lats_d01 = ds_d01['lat'][:]
lons_d01 = ds_d01['lon'][:]
lats_d02 = ds_d02['lat'][:]
lons_d02 = ds_d02['lon'][:]
LH_d01 = np.where(np.isnan(ds_d01['d01_LH'][:]), 0 ,ds_d01['d01_LH'][:])
LH_d02 = np.where(np.isnan(ds_d02['d02_LH'][:]), 0 ,ds_d02['d02_LH'][:])
SH_d01 = np.where(np.isnan(ds_d01['d01_SH'][:]), 0 ,ds_d01['d01_SH'][:])
SH_d02 = np.where(np.isnan(ds_d02['d02_SH'][:]), 0 ,ds_d02['d02_SH'][:])
# ---------------------------------------------------------------------------------------------
# MAPS
# ---------------------------------------------------------------------------------------------
switch_interactive = 0
seismic_mod = cm.get_cmap('seismic',256)
newcmp = ListedColormap(seismic_mod(np.linspace(0.5,1,256)))
x1 = 0.00005
# specify before all subplots
my_fontsize= 10

fig, axes = plt.subplots(3,3)
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
map1= m1.pcolormesh(lons_d01,lats_d01,ds_d01['M2_PREC'][0,:,:],latlon=True,cmap='jet')
axes[0][0].set_title('MERRA-2\n')
axes[0][0].annotate('(a) \n ', xy=(x1,  axes[0][0].get_ylim()[1]),annotation_clip=False)
mstats = 'm = %.2f, s = %.2f' % (np.nanmean(ds_d01['M2_PREC'][0,:,:]), np.nanstd(ds_d01['M2_PREC'][0,:,:]))
plt.text(1.0,1.0, mstats,  horizontalalignment='right', verticalalignment='bottom', transform=axes[0][0].transAxes, fontsize=my_fontsize)
fig.colorbar(map1,ax=axes[0][0], label='mm/day',extend = 'max', shrink = 0.6)
map1.set_clim(0,8)

m2 = Basemap(projection= 'lcc', lat_0 = lats_d01.mean(),
            lon_0 = lons_d01.mean(),
            llcrnrlon=(-123), llcrnrlat=(58),
            urcrnrlon=(-108),urcrnrlat = (62),
            resolution = 'f',ax=axes[0][1])
m2.drawcoastlines(linewidth=0.3)
m2.drawcountries()
m2.drawparallels(np.arange(-50,70,5))
m2.drawmeridians(np.arange(-130,-100,5))
m2.drawparallels(np.arange(-50,70,5), linewidth=0.5,labels = [True,False,False,False])
m2.drawmeridians(np.arange(-130,-100,5), linewidth=0.5,labels = [False,False,False, True])
map2 = m2.pcolormesh(lons_d01,lats_d01,ds_d01['d01_PREC'][0,:,:],latlon=True,cmap='jet')
axes[0][1].set_title('d01\n')
fig.colorbar(map2, ax=axes[0][1],label='mm/day',extend = 'max', shrink = 0.6)
map2.set_clim(0,8)
axes[0][1].annotate('(b) \n ', xy=(x1,  axes[0][1].get_ylim()[1]),annotation_clip=False)
mstats = 'm = %.2f, s = %.2f' % (np.nanmean(ds_d01['d01_PREC'][0,:,:]), np.nanstd(ds_d01['d01_PREC'][0,:,:]))
plt.text(1.0,1.0, mstats,  horizontalalignment='right', verticalalignment='bottom', transform=axes[0][1].transAxes, fontsize=my_fontsize)

m3 = Basemap(projection= 'lcc', lat_0 = lats_d02.mean(),
            lon_0 = lons_d02.mean(),
            llcrnrlon=(-123), llcrnrlat=(58),
            urcrnrlon=(-108),urcrnrlat = (62),
            resolution = 'f',ax=axes[0][2])
m3.drawcoastlines(linewidth=0.3)
m3.drawcountries()
m3.drawparallels(np.arange(-50,70,5))
m3.drawmeridians(np.arange(-130,-100,5))
m3.drawparallels(np.arange(-50,70,5), linewidth=0.5,labels = [True,False,False,False])
m3.drawmeridians(np.arange(-130,-100,5), linewidth=0.5,labels = [False,False,False, True])
map3=m3.pcolormesh(lons_d02,lats_d02,ds_d02['d02_PREC'][0,:,:],latlon=True,cmap='jet')
axes[0][2].set_title('d02\n')
plt.colorbar(map3,ax=axes[0][2],label='mm/day',extend = 'max', shrink = 0.6)
map3.set_clim(0,8)
axes[0][2].annotate('(c) \n ', xy=(x1,  axes[0][2].get_ylim()[1]),annotation_clip=False)
mstats = 'm = %.2f, s = %.2f' % (np.nanmean(ds_d02['d02_PREC'][0,:,:]), np.nanstd(ds_d02['d02_PREC'][0,:,:]))
plt.text(1.0,1.0, mstats,  horizontalalignment='right', verticalalignment='bottom', transform=axes[0][2].transAxes, fontsize=my_fontsize)


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
map4=m4.pcolormesh(lons_d01,lats_d01,ds_d01['M2_LH'][0,:,:],latlon=True,cmap='jet')
plt.colorbar(map4,ax=axes[1][0],label='W/m²',extend = 'max', shrink = 0.6)
map4.set_clim(0,100)
axes[1][0].annotate('(d) \n ', xy=(x1,  axes[1][0].get_ylim()[1]),annotation_clip=False)
mstats = 'm = %.2f, s = %.2f' % (np.nanmean(ds_d01['M2_LH'][0,:,:]), np.nanstd(ds_d01['M2_LH'][0,:,:]))
plt.text(1.0,1.0, mstats,  horizontalalignment='right', verticalalignment='bottom', transform=axes[1][0].transAxes, fontsize=my_fontsize)

m5 = Basemap(projection= 'lcc', lat_0 = lats_d01.mean(),
            lon_0 = lons_d01.mean(),
            llcrnrlon=(-123), llcrnrlat=(58),
            urcrnrlon=(-108),urcrnrlat = (62),
            resolution = 'f', ax=axes[1][1])
m5.drawcoastlines(linewidth=0.3)
m5.drawcountries()
m5.drawparallels(np.arange(-50,70,5))
m5.drawmeridians(np.arange(-130,-100,5))
m5.drawparallels(np.arange(-50,70,5), linewidth=0.5,labels = [True,False,False,False])
m5.drawmeridians(np.arange(-130,-100,5), linewidth=0.5,labels = [False,False,False, True])
map5=m5.pcolormesh(lons_d01,lats_d01,LH_d01[0,:,:],latlon=True,cmap='jet')
plt.colorbar(map5, ax=axes[1][1], label='W/m²',extend = 'both', shrink = 0.6)
map5.set_clim(0,100)
axes[1][1].annotate('(e) \n ', xy=(x1,  axes[1][1].get_ylim()[1]),annotation_clip=False)
mstats = 'm = %.2f, s = %.2f' % (np.nanmean(LH_d01[0,:,:]), np.nanstd(LH_d01[0,:,:]))
plt.text(1.0,1.0, mstats,  horizontalalignment='right', verticalalignment='bottom', transform=axes[1][1].transAxes, fontsize=my_fontsize)

m6 = Basemap(projection= 'lcc', lat_0 = lats_d02.mean(),
            lon_0 = lons_d02.mean(),
            llcrnrlon=(-123), llcrnrlat=(58),
            urcrnrlon=(-108),urcrnrlat = (62),
            resolution = 'f', ax=axes[1][2])
m6.drawcoastlines(linewidth=0.3)
m6.drawcountries()
m6.drawparallels(np.arange(-50,70,5))
m6.drawmeridians(np.arange(-130,-100,5))
m6.drawparallels(np.arange(-50,70,5), linewidth=0.5,labels = [True,False,False,False])
m6.drawmeridians(np.arange(-130,-100,5), linewidth=0.5,labels = [False,False,False, True])
map6 =m6.pcolormesh(lons_d02,lats_d02,LH_d02[0,:,:],latlon=True,cmap='jet')
plt.colorbar(map6,ax=axes[1][2], label='W/m²',extend = 'both', shrink = 0.6)
map6.set_clim(0,100)
axes[1][2].annotate('(f) \n ', xy=(x1,  axes[1][2].get_ylim()[1]),annotation_clip=False)
mstats = 'm = %.2f, s = %.2f' % (np.nanmean(LH_d02[0,:,:]), np.nanstd(LH_d02[0,:,:]))
plt.text(1.0,1.0, mstats,  horizontalalignment='right', verticalalignment='bottom', transform=axes[1][2].transAxes, fontsize=my_fontsize)

m7 = Basemap(projection= 'lcc', lat_0 = lats_d01.mean(),
            lon_0 = lons_d01.mean(),
            llcrnrlon=(-123), llcrnrlat=(58),
            urcrnrlon=(-108),urcrnrlat = (62),
            resolution = 'f',ax=axes[2][0])
m7.drawcoastlines(linewidth=0.3)
m7.drawcountries()
m7.drawparallels(np.arange(-50,70,5))
m7.drawmeridians(np.arange(-130,-100,5))
m7.drawparallels(np.arange(-50,70,5), linewidth=0.5,labels = [True,False,False,False])
m7.drawmeridians(np.arange(-130,-100,5), linewidth=0.5,labels = [False,False,False, True])
map7=m7.pcolormesh(lons_d01,lats_d01,ds_d01['M2_SH'][0,:,:],latlon=True,cmap='jet')
plt.colorbar(map7,ax=axes[2][0],label='W/m²',extend = 'max', shrink = 0.6)
map7.set_clim(0,50)
axes[2][0].annotate('(d) \n ', xy=(x1,  axes[2][0].get_ylim()[1]),annotation_clip=False)
mstats = 'm = %.2f, s = %.2f' % (np.nanmean(ds_d01['M2_SH'][0,:,:]), np.nanstd(ds_d01['M2_SH'][0,:,:]))
plt.text(1.0,1.0, mstats,  horizontalalignment='right', verticalalignment='bottom', transform=axes[2][0].transAxes, fontsize=my_fontsize)

m8 = Basemap(projection= 'lcc', lat_0 = lats_d01.mean(),
            lon_0 = lons_d01.mean(),
            llcrnrlon=(-123), llcrnrlat=(58),
            urcrnrlon=(-108),urcrnrlat = (62),
            resolution = 'f', ax=axes[2][1])
m8.drawcoastlines(linewidth=0.3)
m8.drawcountries()
m8.drawparallels(np.arange(-50,70,5))
m8.drawmeridians(np.arange(-130,-100,5))
m8.drawparallels(np.arange(-50,70,5), linewidth=0.5,labels = [True,False,False,False])
m8.drawmeridians(np.arange(-130,-100,5), linewidth=0.5,labels = [False,False,False, True])
map8=m8.pcolormesh(lons_d01,lats_d01,SH_d01[0,:,:],latlon=True,cmap='jet')
plt.colorbar(map8, ax=axes[2][1], label='W/m²',extend = 'both', shrink = 0.6)
map8.set_clim(0,50)
axes[2][1].annotate('(e) \n ', xy=(x1,  axes[2][1].get_ylim()[1]),annotation_clip=False)
mstats = 'm = %.2f, s = %.2f' % (np.nanmean(SH_d01[0,:,:]), np.nanstd(SH_d01[0,:,:]))
plt.text(1.0,1.0, mstats,  horizontalalignment='right', verticalalignment='bottom', transform=axes[2][1].transAxes, fontsize=my_fontsize)

m9 = Basemap(projection= 'lcc', lat_0 = lats_d02.mean(),
            lon_0 = lons_d02.mean(),
            llcrnrlon=(-123), llcrnrlat=(58),
            urcrnrlon=(-108),urcrnrlat = (62),
            resolution = 'f', ax=axes[2][2])
m9.drawcoastlines(linewidth=0.3)
m9.drawcountries()
m9.drawparallels(np.arange(-50,70,5))
m9.drawmeridians(np.arange(-130,-100,5))
m9.drawparallels(np.arange(-50,70,5), linewidth=0.5,labels = [True,False,False,False])
m9.drawmeridians(np.arange(-130,-100,5), linewidth=0.5,labels = [False,False,False, True])
map9 =m9.pcolormesh(lons_d02,lats_d02,SH_d02[0,:,:],latlon=True,cmap='jet')
plt.colorbar(map9,ax=axes[2][2], label='W/m²',extend = 'both', shrink = 0.6)
map9.set_clim(0,50)
axes[2][2].annotate('(f) \n ', xy=(x1,  axes[2][2].get_ylim()[1]),annotation_clip=False)
mstats = 'm = %.2f, s = %.2f' % (np.nanmean(SH_d02[0,:,:]), np.nanstd(SH_d02[0,:,:]))
plt.text(1.0,1.0, mstats,  horizontalalignment='right', verticalalignment='bottom', transform=axes[2][2].transAxes, fontsize=my_fontsize)

plt.show()
