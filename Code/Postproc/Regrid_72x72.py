## Load Modules

## This script regrids d02 resolution data to d01 resolution data

infile_d01 = '/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/wrfout_nc_files/data_calibrated_ax_d01.nc'
outfile = '/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/wrfout_nc_files/regrid_72x72n.nc'

LIs = ['LTG3', 'LPI', 'PR92W', 'CAPExP_R', 'CAPExP_CSI']

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

## Create new nc file
ds01 = Dataset(infile_d01, 'r')
trg = Dataset(outfile, mode='w')

lat_d01 = ds01['lat'][:]
lon_d01 = ds01['lon'][:]
lats = np.linspace(np.nanmin(lat_d01), np.nanmax(lat_d01),8)
lons = np.linspace(np.nanmin(lon_d01), np.nanmax(lon_d01),14)

lat,lon = np.meshgrid(lons,lats)

trg.createDimension('time', 13248)
trg.createDimension('lat', 8)
trg.createDimension('lon', 14)
trg.createVariable('time','int', dimensions=('time',),zlib=True)
trg.createVariable('lat',lat.dtype, dimensions=('lat','lon',),zlib=True)
trg.variables['lat'][:] = lat[:]
trg.createVariable('lon',lon.dtype, dimensions=('lat','lon',),zlib=True)
trg.variables['lon'][:] = lon[:]
trg.createVariable('LPI_d01', 'f4', dimensions=('time','lat','lon',), zlib=True)
trg.createVariable('LTG3_d01', 'f4', dimensions=('time','lat','lon',), zlib=True)
trg.createVariable('PR92W_d01', 'f4', dimensions=('time','lat','lon',), zlib=True)
trg.createVariable('CAPExP_R_d01', 'f4', dimensions=('time','lat','lon',), zlib=True)
trg.createVariable('CAPExP_CSI_d01', 'f4', dimensions=('time','lat','lon',), zlib=True)
trg.createVariable('LPI_d02', 'f4', dimensions=('time','lat','lon',), zlib=True)
trg.createVariable('LTG3_d02', 'f4', dimensions=('time','lat','lon',), zlib=True)
trg.createVariable('PR92W_d02', 'f4', dimensions=('time','lat','lon',), zlib=True)
trg.createVariable('CAPExP_R_d02', 'f4', dimensions=('time','lat','lon',), zlib=True)
trg.createVariable('CAPExP_CSI_d02', 'f4', dimensions=('time','lat','lon',), zlib=True)
trg.createVariable('Obs', 'f4', dimensions=('time','lat','lon',), zlib=True)

for i in range(0,8):
    for j in range(0,14):
        k=i*8
        l = (i+1)*8
        m = j*8
        n = (j+1)*8
        trg['CAPExP_R_d01'][:,i,j] = np.nanmean(ds01['CAPExP_R_d01'][:,k:l,m:n], axis=(1,2))
        trg['CAPExP_CSI_d01'][:, i, j] = np.nanmean(ds01['CAPExP_CSI_d01'][:, k:l, m:n], axis=(1, 2))
        trg['LTG3_d01'][:, i, j] = np.nanmean(ds01['LTG3_d01'][:, k:l, m:n], axis=(1, 2))
        trg['LPI_d01'][:, i, j] = np.nanmean(ds01['LPI_d01'][:, k:l, m:n], axis=(1, 2))
        trg['PR92W_d01'][:, i, j] = np.nanmean(ds01['PR92W_d01'][:, k:l, m:n], axis=(1, 2))
        trg['Obs'][:, i, j] = np.nanmean(ds01['Obs'][:, k:l, m:n], axis=(1, 2))
        trg['CAPExP_R_d02'][:,i,j] = np.nanmean(ds01['CAPExP_R_d02'][:,k:l,m:n], axis=(1,2))
        trg['CAPExP_CSI_d02'][:, i, j] = np.nanmean(ds01['CAPExP_CSI_d02'][:, k:l, m:n], axis=(1, 2))
        trg['LTG3_d02'][:, i, j] = np.nanmean(ds01['LTG3_d02'][:, k:l, m:n], axis=(1, 2))
        trg['LPI_d02'][:, i, j] = np.nanmean(ds01['LPI_d02'][:, k:l, m:n], axis=(1, 2))
        trg['PR92W_d02'][:, i, j] = np.nanmean(ds01['PR92W_d02'][:, k:l, m:n], axis=(1, 2))
trg.close()
