# load modules
import numpy as np
import pandas as pd
import os
from netCDF4 import Dataset
from wrf import getvar, latlon_coords
import datetime as dt
from pyproj import Proj, transform
import sys
from matplotlib import pyplot as plt
import netCDF4 as nc
import xarray as xr
from scipy import ndimage
from scipy.spatial import KDTree
from osgeo import gdal
from math import sin, cos, sqrt, atan2

# load data
df = pd.read_csv('/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/screen-level-obs/all_data.csv')
ds = Dataset('/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/wrfout_nc_files/wrfout_d01_2015.nc', 'r')
outfile = '/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/screen-level-obs/obs_at_d01.nc'

# make a new netCDF file based on the same grid as the d01 file.
def create_file_from_source(src_file, trg_file):
    src = ds
    trg = nc.Dataset(trg_file, mode='w')

    # Create the dimensions of the file
    for name, dim in src.dimensions.items():
        trg.createDimension(name, len(dim) if not dim.isunlimited() else None)

    # Copy the global attributes
    trg.setncatts({a:src.getncattr(a) for a in src.ncattrs()})

    # Create the variables in the file
    for name, var in src.variables.items():
        trg.createVariable(name, var.dtype, var.dimensions)

        # Copy the variable attributes
        trg.variables[name].setncatts({a:var.getncattr(a) for a in var.ncattrs()})

        # Copy the variables values (as 'f4' eventually)
        trg.variables[name][:] = src.variables[name][:]

    # Save the file
    trg.close()

# destination grid, create output file from input ds01 file
create_file_from_source(ds,outfile)

ds_out = Dataset(outfile, 'a')
# set all values to nan
ds_out['RH2'][:] = np.nan
ds_out['T2'][:] = np.nan

# drop the unnecessary columns to keep overview
final_table_columns = ['Date/Time (LST)', 'Latitude (y)', 'Longitude (x)', 'Temp (°C)', 'Rel Hum (%)', 'Station Name']
df.drop(columns=[col for col in df if col not in final_table_columns], inplace=True)

# Make sure the times of the csv file match those of the model
df = df[pd.DatetimeIndex(pd.to_datetime(df['Date/Time (LST)'])).year == 2015].reset_index(drop=True)
df = df[pd.DatetimeIndex(pd.to_datetime(df['Date/Time (LST)'])).month != 5].reset_index(drop=True)
df = df.sort_values(['Latitude (y)','Date/Time (LST)']).reset_index(drop=True)

# Create function to calculate the distance between two points on the world
radius_earth = 6471 # in km

def haversine_distance(lat1, lon1, lat2, lon2):
   r = radius_earth
   phi1 = np.radians(lat1)
   phi2 = np.radians(lat2)
   delta_phi = np.radians(lat2 - lat1)
   delta_lambda = np.radians(lon2 - lon1)
   a = np.sin(delta_phi / 2)**2 + np.cos(phi1) * np.cos(phi2) *   np.sin(delta_lambda / 2)**2
   res = r * (2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a)))
   return res

distance = np.zeros((26,140,200))
distance[distance[:]==0] = np.nan
time = pd.to_timedelta(ds['time'][:-1], unit='h') + pd.to_datetime('2000010100',format='%Y%m%d%H')

# Loop over all stations and calculate the nearest neighbor on the d01 grid.
for i in range(0,len(np.unique(df['Station Name']))):
        print('loop ' + str(i) + ': ' + str(np.unique(df['Station Name'])[i]))
        inds_all = np.asarray(np.where((df['Station Name'] == np.unique(df['Station Name'])[i])))
        print(df['Latitude (y)'][inds_all[0,0]])
        print(df['Longitude (x)'][inds_all[0,0]])
        for l in range(0,140):
            for k in range(0,200):
                distance[i, l, k] = haversine_distance(df['Latitude (y)'][inds_all[0,0]], df['Longitude (x)'][inds_all[0,0]], ds['lat'][l,k], ds['lon'][l,k])
                inds = np.asarray(np.where(distance[i,:,:] == np.nanmin(distance[i,:,:])))
        print(np.asarray(df['Temp (°C)'][inds_all[0,:]]))
        ds_out['T2'][:, inds[0,0], inds[1,0]] = np.asarray(df['Temp (°C)'][inds_all[0,:]])
        print(ds_out['T2'][:, inds[0,0], inds[1,0]])
        ds_out['RH2'][:, inds[0,0], inds[1,0]] = np.asarray(df['Rel Hum (%)'][inds_all[0,:]])

ds_out.close()
