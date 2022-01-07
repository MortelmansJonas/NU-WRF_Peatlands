#!/usr/bin/env python
# ---------------------------------------------------------------------------------------------
# MODULES
# ---------------------------------------------------------------------------------------------
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

# ---------------------------------------------------------------------------------------------
# PUT ALL VARIABLES IN THE FILE
# ---------------------------------------------------------------------------------------------
infile = '/scratch/leuven/336/vsc33651/nu-wrf-dev/wrfout_nc_files/domain1_all.nc'
ds01 = Dataset(infile, 'r')
outfile = '/scratch/leuven/317/vsc31786/nu-wrf-dev/wrfout_nc_files/domain2_at_domain1_all_v4.nc'
create_file_from_source(infile,outfile)
ds_out = nc.Dataset(outfile,"a")

LI = 'LTG3'
infile = '/scratch/leuven/317/vsc31786/nu-wrf-dev/wrfout_nc_files/domain2_at_domain1_'+LI+'.nc'
ds_in = Dataset(infile, 'r')
ds_out[LI][:,:,:] = ds_in[LI][:,:,:]

LI = 'PR92H'
infile = '/scratch/leuven/317/vsc31786/nu-wrf-dev/wrfout_nc_files/domain2_at_domain1_'+LI+'.nc'
ds_in = Dataset(infile, 'r')
ds_out[LI][:,:,:] = ds_in[LI][:,:,:]

LI = 'PR92W'
infile = '/scratch/leuven/317/vsc31786/nu-wrf-dev/wrfout_nc_files/domain2_at_domain1_'+LI+'.nc'
ds_in = Dataset(infile, 'r')
ds_out[LI][:,:,:] = ds_in[LI][:,:,:]

LI = 'LPI'
infile = '/scratch/leuven/317/vsc31786/nu-wrf-dev/wrfout_nc_files/domain2_at_domain1_'+LI+'.nc'
ds_in = Dataset(infile, 'r')
ds_out[LI][:,:,:] = ds_in[LI][:,:,:]

LI = 'T2'
infile = '/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/wrfout_nc_files/domain2_at_domain1_'+LI+'.nc'
ds_in = Dataset(infile, 'r')
ds_out[LI][:,:,:] = ds_in[LI][:,:,:]

LI = 'RH2'
infile = '/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/wrfout_nc_files/domain2_at_domain1_'+LI+'.nc'
ds_in = Dataset(infile, 'r')
ds_out[LI][:,:,:] = ds_in[LI][:,:,:]

ds_out.close()


