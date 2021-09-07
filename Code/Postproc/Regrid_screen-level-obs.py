## Load Modules

## This script regrids d02 resolution data to d01 resolution data

infile1 = '/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/screen-level-obs/all_data.nc'
infile2 = '/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/wrfout_nc_files/wrfout_d01_2015.nc'
outfile = '/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/screen-level-obs/obs_at_d01.nc'
# lightning index to regrid
LIs = ['T', 'RH']

import numpy as np
import pandas as pd
import os
from netCDF4 import Dataset
from wrf import getvar, latlon_coords
import datetime as dt
from pyproj import Proj, transform
import sys
from pyresample.geometry import SwathDefinition
from pyresample.kd_tree import resample_nearest
from pyresample import bilinear, geometry
from matplotlib import pyplot as plt
import netCDF4 as nc
import xarray as xr
from scipy import ndimage

ds01 = Dataset(infile2, 'r')

# Create outfile with same resolution as d01
nf = Dataset(outfile, "w", format = 'NETCDF4')

# Create dimensions
time = nf.createDimension('time', None)
lat = nf.createDimension('lat', None)
lon = nf.createDimension('lon', None)

# Create variables
nf.createVariable('time', 'int', dimensions=('time'),zlib=True)
nf.variables['time'][:] = ds01['time'][:]
nf.createVariable('lat','f', dimensions=('lat'))
nf.variables['lat'][:] = np.unique(ds01['lat'][:])
nf.createVariable('lon', 'f4', dimensions=('lon'))
nf.variables['lon'][:] = np.unique(ds01['lon'][:])
nf.createVariable('T', 'f4', dimensions= ('time', 'lat', 'lon'),zlib=True)
nf.createVariable('RH', 'f4', dimensions= ('time', 'lat', 'lon'),zlib=True)

# Set metadata
nf.variables['time'].setncatts({'long_name': 'Time', 'units': ''})
nf.variables['lat'].setncatts({'long_name': 'Latitude', 'units': 'degrees_north'})
nf.variables['lon'].setncatts({'long_name': 'Longitude', 'units': 'degrees_east'})
nf.variables['T'].setncatts({'long_name': 'Dry Bulb Temperature', 'units': '°C'})
nf.variables['RH'].setncatts({'long_name': 'Relative Humidity', 'units': '%'})

# read source dataset at original resolution
ds02 = Dataset(infile1, 'r')
lat_out     = ds01.variables['lat'][:].data
lon_out     = ds01.variables['lon'][:].data
lat_in     = ds02.variables['lat'][:].data
lon_in     = ds02.variables['lon'][:].data

x, y = np.meshgrid(lon_in, lat_in)
def_a = SwathDefinition(lons=lon_out, lats=lat_out)
def_b = SwathDefinition(lons=x, lats=y)

for LI in LIs:
    for t in range(0,np.shape(ds02[LI])[0]):
        print(LI + " " + str(t) + " / " + str(np.shape(ds02[LI])[0]))
        nf[LI][t,:,:] = resample_nearest(def_b,ds02[LI][t,:,:],def_a,radius_of_influence = 70000,fill_value = np.nan)

nf.close()
