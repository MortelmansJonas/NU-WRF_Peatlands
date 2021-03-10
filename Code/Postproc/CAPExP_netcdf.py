## Load the necessary modules
import numpy as np
import pandas as pd
import netCDF4 as nc
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.axis as axis
from matplotlib.cm import get_cmap
from matplotlib.colors import from_levels_and_colors, ListedColormap, LinearSegmentedColormap
from numpy import meshgrid
from matplotlib import cm
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
from datetime import date, timedelta
from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()
from netCDF4 import Dataset
from wrf import (getvar, latlon_coords)
import os

f2015_d01 = '/staging/leuven/stg_00024/OUTPUT/michelb/nu-wrf-dev/Great_Slave_Lake/2015/WRFPRS_d01_2015.nc'
ds2015_d01 = Dataset(f2015_d01, 'r')

ds_wrf = Dataset('/scratch/leuven/336/vsc33651/nu-wrf-dev/wrfout_nc_files/wrfout_d01_2015.nc', 'r')
lats_2015_d01 = ds_wrf.variables['lat'][:]
lons_2015_d01 = ds_wrf.variables['lon'][:]
time_2015_d01 = ds_wrf.variables['time'][:]

ds_CP = Dataset('/scratch/leuven/336/vsc33651/nu-wrf-dev/wrfout_nc_files/2015_d01_CP.nc', mode='w', format='NETCDF4')
ds_CP.createDimension('time', 2208)
ds_CP.createDimension('lat', 140)
ds_CP.createDimension('lon', 200)
ds_CP.createVariable('time','int', dimensions=('time',),zlib=True)
ds_CP.variables['time'][:] = time_2015_d01[:]
ds_CP.createVariable('lat',lats_2015_d01.dtype, dimensions=('lat','lon'),zlib=True)
ds_CP.variables['lat'][:] = lats_2015_d01
ds_CP.createVariable('lon',lons_2015_d01.dtype, dimensions=('lat','lon'),zlib=True)
ds_CP.variables['lon'][:] = lons_2015_d01
ds_CP.createVariable('CAPExP','f4', dimensions=('time','lat','lon'),zlib=True)
ds_CP.variables['time'].setncatts({'long_name': 'time', 'units': 'hours since 2000-01-01 00:00'})
ds_CP.variables['lat'].setncatts({'long_name': 'latitude',  'units': 'degrees_north'})
ds_CP.variables['lon'].setncatts({'long_name': 'longitude', 'units': 'degrees_east'})
ds_CP.variables['CAPExP'].setncatts({'long_name': 'CAPE x P proxy',  'units': 'W m^{-2}'})
ds_CP['CAPExP'][:] = ds2015_d01['CAPE'][:]* ds_wrf['RAINC'][:]
ds_CP.close()