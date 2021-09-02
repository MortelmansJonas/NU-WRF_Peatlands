## Load Modules

## This script regrids d02 resolution data to d01 resolution data

infile_d01 = '/scratch/leuven/336/vsc33651/nu-wrf-dev/wrfout_nc_files/domain1_all.nc'
infile_d02 = '/scratch/leuven/336/vsc33651/nu-wrf-dev/wrfout_nc_files/domain2_all.nc'
outfile = '/scratch/leuven/317/vsc31786/nu-wrf-dev/wrfout_nc_files/domain2_at_domain1_LTG3.nc'
lightning index to regrid
LIs = ['LTG3']

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

def create_file_from_source(src_file, trg_file):
    src = nc.Dataset(src_file)
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
ds01 = Dataset(infile_d01, 'r')
create_file_from_source(infile,outfile)

# read source dataset at d02 resolution
ds02 = Dataset(infile_d02, 'r')
lat_out     = ds01.variables['lat'][:].data
lon_out     = ds01.variables['lon'][:].data
lat_in     = ds02.variables['lat'][:].data
lon_in     = ds02.variables['lon'][:].data

def_a = SwathDefinition(lons=lon_out, lats=lat_out)
def_b = SwathDefinition(lons=lon_in, lats=lat_in)

ds_out = nc.Dataset(outfile,"a")
for LI in LIs:
    for t in range(0,np.shape(ds02[LI])[0]):
        print(LI + " " + str(t) + " / " + str(np.shape(ds02[LI])[0]))
        # first do a 3x3 averaging, i.e. 9 km to 3 km resolution
        result = ndimage.generic_filter(ds02[LI][t,:,:].data, np.nanmean, size=3, mode='constant', cval=np.NaN)
        # then apply nearest neighbor search, ... worth checking why this is actually needed, ...
        # from d02 to d01 the 3x3 aggregation is probably already resulting in a grid that matches perfectly the d01 grid
        ds_out[LI][t,:,:] = resample_nearest(def_b,result,def_a,radius_of_influence = 70000,fill_value = np.nan)

ds_out.close()


