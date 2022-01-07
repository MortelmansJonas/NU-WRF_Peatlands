#!/usr/bin/env python
# ---------------------------------------------------------------------------------------------
# MODULES
# ---------------------------------------------------------------------------------------------
## This script regrids d02 resolution data to d01 resolution data

infile_d01 = '/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/wrfout_nc_files/d01_all.nc'
infile_d02 = '/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/wrfout_nc_files/d02_all.nc'
outfile = '//scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/wrfout_nc_files/domain2_at_domain1_CAPExP_CSI.nc'

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

# ---------------------------------------------------------------------------------------------
# CREATE .nc FILE FROM TARGET FILE (FILE WITH DESIRED RESOLUTION)
# ---------------------------------------------------------------------------------------------
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

#  destination grid, create output file from input ds01 file
ds01 = Dataset(infile_d01, 'r')
#create_file_from_source(infile_d01,outfile)

# read source dataset at d02 resolution
ds02 = Dataset(infile_d02, 'r')
lat_out     = ds01.variables['lat'][:].data
lon_out     = ds01.variables['lon'][:].data
lat_in     = ds02.variables['lat'][:].data
lon_in     = ds02.variables['lon'][:].data

def_a = SwathDefinition(lons=lon_out, lats=lat_out)
def_b = SwathDefinition(lons=lon_in, lats=lat_in)

# Set all other variables in the file to nan to avoid issues later on
ds_out = nc.Dataset(outfile,"a")
ds_out['LTG3'][:] = np.nan
ds_out['LPI'][:] = np.nan
ds_out['PR92W'][:] = np.nan
ds_out['CAPExP_R'][:] = np.nan
ds_out['T2'][:] = np.nan
ds_out['RH2'][:] = np.nan
ds_out['LH'][:] = np.nan
ds_out['SH'][:] = np.nan
ds_out['RAINC'][:] = np.nan
ds_out['RAINNC'][:] = np.nan

for t in range(0,np.shape(ds02['LTG3'])[0]):
    print(str(t) + " / " + str(np.shape(ds02['LTG3'])[0]))
    # first do a 3x3 averaging, i.e. 9 km to 3 km resolution
    result_CAPExP_CSI = ndimage.generic_filter(ds02['CAPExP_CSI'][t, :, :].data, np.nanmean, size=3, mode='constant',cval=np.NaN)
    # then apply nearest neighbor search, ... worth checking why this is actually needed, ...
    # from d02 to d01 the 3x3 aggregation is probably already resulting in a grid that matches perfectly the d01 grid
    ds_out['CAPExP_CSI'][t, :, :] = resample_nearest(def_b, result_CAPExP_CSI, def_a, radius_of_influence=70000, fill_value=np.nan)

ds_out.close()
