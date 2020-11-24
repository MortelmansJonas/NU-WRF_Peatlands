### First import the necessary libraries
import numpy as np
import cartopy
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from wrf import (getvar, interplevel, vertcross, vinterp, ALL_TIMES, latlon_coords, get_cartopy, to_np,
                cartopy_xlim, cartopy_ylim, CoordPair)
import os
import cartopy.feature as cfeature
from cartopy.feature import NaturalEarthFeature, COLORS
from matplotlib.cm import get_cmap
from cartopy import crs
from matplotlib.colors import from_levels_and_colors

### Obtain the WRF Output files
WRF_DIRECTORY = '/scratch/leuven/336/vsc33651/nu-wrf-dev/testcases/wrflis/noahmp36_modis_merra2'
WRF_FILES = ["wrfout_d01_2005-07-01_12:00:00"]

_WRF_FILES = [os.path.abspath(os.path.join(WRF_DIRECTORY,f)) for f in WRF_FILES]

## Check if the files exist
try:
    for f in _WRF_FILES:
        if not os.path.exists(f):
            raise ValueError("{} does not exist. "
                             "Check for typos or incorrect directory.".format(f))
except ValueError as e:
    raise e

# Create functions so that the WRF files only need to be specified using the WRF_FILES global above
def single_wrf_file():
    global _WRF_FILES
    return _WRF_FILES[0]
def multiple_wrf_files():
    global _WRF_FILES
    return _WRF_FILES
print('All tests passed!')


file_path = single_wrf_file()
wrf_file = Dataset(file_path)

file_path = single_wrf_file()
wrf_file = Dataset(file_path)
print('\n')
print('netcdf4: wrf file')
print(wrf_file)

# ## Create a topography map
# # Get the terrain height
# terrain = getvar(wrf_file, "ter", timeidx=0)
#
# # Get the cartopy object and the lat,lon coords
# cart_proj = get_cartopy(terrain)
# lats, lons = latlon_coords(terrain)
#
# # Create a figure and get the GetAxes object
# fig = plt.figure(figsize=(10, 7.5))
# ax = plt.axes(projection=cart_proj)
#
# # Set the contour levels
# levels = np.arange(0., 1250, 40.)
#
# plt.contourf(to_np(lons), to_np(lats), to_np(terrain), levels=levels, transform=crs.PlateCarree(), cmap=get_cmap("terrain"))
#
# # Add a color bar.
# plt.colorbar(ax=ax, shrink=.99)
# plt.title('Topography (height in m)')
# # Download and add the stateborders and coastlines
# provinces_50m = cfeature.NaturalEarthFeature('cultural', 'admin_1_states_provinces_lines', '50m', facecolor='none')
# ax.gridlines()
# ax.add_feature(cfeature.LAKES)
# ax.add_feature(cfeature.BORDERS)
# ax.add_feature(cfeature.COASTLINE)
# ax.add_feature(cfeature.RIVERS)
# ax.add_feature(provinces_50m)
#
# plt.show()
# plt.savefig('Topography.png',dpi=300)

# ## Create full plot of Sea Level Pressure
# slp = getvar(wrf_file, "slp", timeidx=0)
# # Get the cartopy object and the lat,lon coords
# cart_proj = get_cartopy(slp)
# lats, lons = latlon_coords(slp)
# # Create a figure and get the GetAxes object
# fig = plt.figure(figsize=(10, 7.5))
# geo_axes = plt.axes(projection=cart_proj)
#
# geo_axes.gridlines()
# geo_axes.add_feature(cfeature.LAKES)
# geo_axes.add_feature(cfeature.BORDERS)
# geo_axes.add_feature(cfeature.COASTLINE)
# geo_axes.add_feature(cfeature.RIVERS)
# geo_axes.add_feature(provinces_50m)
#
# # Set the contour levels so that all plots match
# levels = np.arange(980., 1030., 2.5)
# # Make the contour lines and fill them.
# plt.contour(to_np(lons), to_np(lats),
#                to_np(slp), levels=levels, colors="black",
#                transform=crs.PlateCarree())
# plt.contourf(to_np(lons), to_np(lats),
#                 to_np(slp), levels=levels,
#                 transform=crs.PlateCarree(),
#                 cmap=get_cmap("jet"))
# # Add a color bar. The shrink often needs to be set by trial and error.
# plt.colorbar(ax=geo_axes, shrink=.86)
# plt.title('Sea Level Presssure (hPa)')
# plt.show()
#
#
# ## Create a map showing the wind barbs
# # Extract the pressure, geopotential height, and wind variables
# p = getvar(wrf_file, "pressure")
# z = getvar(wrf_file, "z", units="dm")
# ua = getvar(wrf_file, "ua", units="kt")
# va = getvar(wrf_file, "va", units="kt")
# wspd = getvar(wrf_file, "wspd_wdir", units="kt")[0,:]
#
# # Interpolate geopotential height, u, and v winds to 850 hPa
# ht_850 = interplevel(z, p, 850)
# u_850 = interplevel(ua, p, 850)
# v_850 = interplevel(va, p, 850)
# wspd_850 = interplevel(wspd, p, 850)
#
# # Get the lat/lon coordinates
# lats, lons = latlon_coords(ht_850)
#
# # Get the map projection information
# cart_proj = get_cartopy(ht_850)
#
# # Create the figure
# fig = plt.figure(figsize=(10,7.5))
# ax = plt.axes(projection=cart_proj)
#
#
# # Add the 850 hPa geopotential height contours
# levels = np.arange(0., 170., 10.)
# contours = plt.contour(to_np(lons), to_np(lats), to_np(ht_850), levels=levels,  colors="black",
#                        transform=crs.PlateCarree())
#
# plt.clabel(contours, inline=1, fontsize=10, fmt="%i")
#
# # Add the wind speed contours
# levels = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120]
# wspd_contours = plt.contourf(to_np(lons), to_np(lats), to_np(wspd_850), levels=levels,
#                                 cmap=get_cmap("rainbow"), transform=crs.PlateCarree())
#
# plt.colorbar(wspd_contours, ax=ax, orientation="horizontal",
#              pad=.05, shrink=.75)
#
# # Add the 850 hPa wind barbs, only plotting 10 barbs in each direction
# # Also, skipping the border barbs.
# thin = [int(x/10.) for x in lons.shape]
# plt.barbs(to_np(lons[::thin[0], ::thin[1]]),
#              to_np(lats[::thin[0], ::thin[1]]),
#              to_np(u_850[::thin[0], ::thin[1]]),
#              to_np(v_850[::thin[0], ::thin[1]]),
#              length=6,
#              transform=crs.PlateCarree())
#
# plt.title("850 MB Height (dm; black line), Wind Speed (kt; colorbar), Barbs (kt)")
# # Download and add the stateborders and coastlines
# provinces_50m = cfeature.NaturalEarthFeature('cultural', 'admin_1_states_provinces_lines', '50m', facecolor='none')
# ax.gridlines()
# ax.add_feature(cfeature.LAKES)
# ax.add_feature(cfeature.BORDERS)
# ax.add_feature(cfeature.COASTLINE)
# ax.add_feature(cfeature.RIVERS)
# ax.add_feature(provinces_50m)
# plt.show()
# plt.savefig('Wind.png',dpi=300)


# ## Cross section panel plot (doesn't work yet)
# # First define the latitude and longitude of your cross section.
# cross_start = CoordPair(lat=61.97, lon=-123.036)
# cross_end = CoordPair(lat=61.97, lon=-109.463)
#
# # Get the WRF variables
# slp = getvar(wrf_file, "slp", timeidx=-1)
# z = getvar(wrf_file, "z", timeidx=-1)
# dbz = getvar(wrf_file, "dbz", timeidx=-1)
# Z = 10**(dbz/10) # Use linear Z for interpolation
#
# # Compute the vertical cross-section interpolation.  Also, include the lat/lon
# # points along the cross-section in the metadata by setting latlon to True.
# z_cross = vertcross(Z, z, wrfin=wrf_file, start_point=cross_start, end_point=cross_end, latlon=True, meta=True)
#
# # Convert back to dBz after interpolation
# dbz_cross = 10.0 * np.log10(z_cross)
#
# # Get the lat/lon points
# lats, lons = latlon_coords(slp)
#
# # Get the cartopy projection object
# cart_proj = get_cartopy(slp)
#
# # Create a figure that will have 2 subplots (1 row, 2 columns)
# fig = plt.figure(figsize=(15,5))
# ax_slp = fig.add_subplot(1,2,1,projection=cart_proj)
# ax_dbz = fig.add_subplot(1,2,2)
#
# # Make the pressure contours
# slp_levels = np.arange(950.,1030.,5)
# slp_contours = ax_slp.contour(to_np(lons), to_np(lats), to_np(slp), levels=slp_levels, colors="black", zorder=3,
#                               transform=crs.PlateCarree())
# # Add contour labels for pressure
# ax_slp.clabel(slp_contours, fmt="%i")
#
# # Draw the cross section line
# ax_slp.plot([cross_start.lon, cross_end.lon],[cross_start.lat, cross_end.lat],color="yellow",marker="o",zorder=3,
#             transform=crs.PlateCarree())
# # Draw the oceans, land, and states
# provinces_50m = cfeature.NaturalEarthFeature('cultural', 'admin_1_states_provinces_lines', '50m', facecolor='none')
# ax_slp.gridlines()
# ax_slp.add_feature(cfeature.LAKES)
# ax_slp.add_feature(cfeature.BORDERS)
# ax_slp.add_feature(cfeature.COASTLINE)
# ax_slp.add_feature(cfeature.RIVERS)
# ax_slp.add_feature(provinces_50m)
#
# # Make the contour plot for dbz
# dbz_levels = np.arange(0.,20.,5.)
# dbz_contours = ax_dbz.contourf(to_np(dbz_cross), levels=dbz_levels, cmap=get_cmap("jet"))
# cb_dbz = fig.colorbar(dbz_contours, ax=ax_dbz)
# cb_dbz.ax.tick_params(labelsize=8)
#
# # Set the x-ticks to use latitude and longitude labels
# coord_pairs = to_np(dbz_cross.coords["xy_loc"])
# x_ticks = np.arange(coord_pairs.shape[0])
# x_labels = [pair.latlon_str() for pair in to_np(coord_pairs)]
# # Only keeping ~5 xticks
# thin = [int(x/5.) for x in x_ticks.shape]
# ax_dbz.set_xticks(x_ticks[1::thin[0]])
# ax_dbz.set_xticklabels(x_labels[::thin[0]], rotation=45, fontsize=8)
#
# # Set the y-ticks to be height
# vert_vals = to_np(dbz_cross.coords["vertical"])
# v_ticks = np.arange(vert_vals.shape[0])
# # Only keeping ~8 vertical ticks
# thin = [int(x/8.) for x in v_ticks.shape]
# ax_dbz.set_yticks(v_ticks[::thin[0]])
# ax_dbz.set_yticklabels(vert_vals[::thin[0]], fontsize=8)
#
# # Set the x-axis and  y-axis labels
# ax_dbz.set_xlabel("Latitude, Longitude", fontsize=12)
# ax_dbz.set_ylabel("Height (m)", fontsize=12)
#
# # Add a title
# ax_slp.set_title("Sea Level Pressure (hPa)", {"fontsize" : 14})
# ax_dbz.set_title("Cross-Section of Reflectivity (dBZ)", {"fontsize" : 14})
#
# plt.show()


## Get Surface Temperature
# t = getvar(wrf_file, "tc", timeidx=0)
# p = getvar(wrf_file, "pressure")
# tc = interplevel(t, p, 850)
# # Get the cartopy object and the lat,lon coords
# cart_proj = get_cartopy(tc)
# lats, lons = latlon_coords(tc)
# # Create a figure and get the GetAxes object
# fig = plt.figure(figsize=(10, 7.5))
# geo_axes = plt.axes(projection=cart_proj)
#
# provinces_50m = cfeature.NaturalEarthFeature('cultural', 'admin_1_states_provinces_lines', '50m', facecolor='none')
# geo_axes.gridlines()
# geo_axes.add_feature(cfeature.LAKES)
# geo_axes.add_feature(cfeature.BORDERS)
# geo_axes.add_feature(cfeature.COASTLINE)
# geo_axes.add_feature(cfeature.RIVERS)
# geo_axes.add_feature(provinces_50m)
#
# # Set the contour levels so that all plots match
# levels = np.arange(0., 24., 2.)
# # Make the contour lines and fill them.
# plt.contour(to_np(lons), to_np(lats),
#                to_np(tc), levels=levels, colors="black",
#                transform=crs.PlateCarree())
# plt.contourf(to_np(lons), to_np(lats),
#                 to_np(tc), levels=levels,
#                 transform=crs.PlateCarree(),
#                 cmap=get_cmap("jet"))
# # Add a color bar. The shrink often needs to be set by trial and error.
# plt.colorbar(ax=geo_axes, shrink=.86)
# plt.title('Tc')
# plt.show()