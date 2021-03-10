# ---------------------------------------------------------------------------------------------
# MODULES
# ---------------------------------------------------------------------------------------------
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import netCDF4 as nc
from matplotlib.colors import LogNorm, from_levels_and_colors, ListedColormap, LinearSegmentedColormap
from matplotlib import cm
from mpl_toolkits.basemap import Basemap
from datetime import timedelta
from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()
from netCDF4 import Dataset
import os
from wrf import getvar, latlon_coords

# ---------------------------------------------------------------------------------------------
# LIS OUTPUT
# ---------------------------------------------------------------------------------------------
date_from = '2005-09-30'
date_to = '2015-06-01'

path = '/staging/leuven/stg_00024/OUTPUT/michelb/nu-wrf-dev/Great_Slave_Lake/2015'
filename_root = 'LIS_RST_NOAHMP36'
domain = 'd01'
output_filename = os.path.join('/scratch/leuven/336/vsc33651/nu-wrf-dev/lis_spinup',filename_root+domain+'.nc')

startdate = pd.to_datetime(date_from)
enddate = pd.to_datetime(date_to)
nmonths = 12* (enddate.year - startdate.year) + (enddate.month - startdate.month)
day = str(startdate).split(' ')
month = np.linspace(1,nmonths,nmonths)
filename = os.path.join(path, filename_root+ '_'+
                        str(startdate.year)+str(startdate.strftime('%m'))+
                        str(startdate.strftime('%d'))+str(2330)+'.'+domain+'.nc')
ds_st = Dataset(filename, 'r')
# layers = np.arange(0,60,1, dtype = 'int')
ds = Dataset(output_filename, mode='w', format='NETCDF4')
print(ds_st)
timeunit = 'Months since 01-2005'
ds.createDimension('time', nmonths)
ds.createDimension('ntiles', np.shape(ds_st['lat'][:])[0])
ds.createDimension('layer_SSTC', 7)
ds.createDimension('layer_SH20',4)
ds.createVariable('time',month.dtype, dimensions=('time',),zlib=True)
ds.variables['time'][:] = month
ds.createVariable('ntiles',ds_st['lat'].dtype, dimensions=('ntiles'),zlib=True)
ds.variables['ntiles'][:] = np.arange(0,26018)
ds.createVariable('ZWT','f4', dimensions=('time','ntiles'),zlib=True)
ds.createVariable('WT','f4', dimensions=('time','ntiles'),zlib=True)
ds.createVariable('WA','f4', dimensions=('time','ntiles'),zlib=True)
ds.createVariable('SMCWTD','f4', dimensions=('time','ntiles'),zlib=True)
ds.createVariable('RECH','f4', dimensions=('time','ntiles'),zlib=True)
ds.createVariable('SH2O','f4', dimensions=('time','ntiles','layer_SH20'),zlib=True)
ds.createVariable('DEEPRECH','f4', dimensions=('time','ntiles'),zlib=True)
ds.createVariable('SSTC','f4', dimensions=('time','ntiles','layer_SSTC'),zlib=True)
ds.createVariable('TG','f4', dimensions=('time','ntiles'),zlib=True)

ds.variables['time'].setncatts({'long_name': 'time', 'units': timeunit})
ds.variables['ntiles'].setncatts({'long_name': 'number of tiles',  'units': '-'})
ds.variables['ZWT'].setncatts({'long_name': 'depth to water table',  'units': 'm'})
ds.variables['WT'].setncatts({'long_name': 'water in aquifer and saturated soil',  'units': 'mm'})
ds.variables['WA'].setncatts({'long_name': 'water storage in aquifer',  'units': 'mm'})
ds.variables['SMCWTD'].setncatts({'long_name': 'soil water content between bottom of the soil and water table',  'units': 'm³ m$^{-3}$'})
ds.variables['RECH'].setncatts({'long_name': 'recharge to or from the water table when shallow',  'units': 'm'})
ds.variables['SH2O'].setncatts({'long_name': 'Volumetric liquid soil moisture',  'units': 'm³ m$^{-3}$'})
ds.variables['DEEPRECH'].setncatts({'long_name': 'recharge to or from the water table when deep',  'units': 'm'})
ds.variables['SSTC'].setncatts({'long_name': 'snow/soil temperature',  'units': 'K'})
ds.variables['TG'].setncatts({'long_name': 'ground temperature (skin temperature)',  'units': 'K'})

for h in range(1,nmonths+1):
    months_passed = pd.to_timedelta(h-1, unit='M')
    a = startdate + months_passed - pd.to_timedelta(1.5, unit='d')
    year = str(a.year)
    month = str(a.strftime('%m'))
    if month == '01' or month == '03' or month == '05' or month == '07' or month == '08' or month == '10' or\
            month == '12':
        day = str(31)
    elif month == '02':
        if year == '2008' or year == '2012':
            day = str(29)
        else:
            day = str(28)
    else:
        day = str(30)
    hour = str(2330)
    filename = os.path.join(path, filename_root + '_' + year + month + day + hour + '.' + domain + '.nc')
    print('processing '+filename)
    ds_in = Dataset(filename, mode='r')

    ds.variables['ZWT'][h-1,:] = ds_in.variables['ZWT'][:]
    ds.variables['WT'][h-1,:] = ds_in.variables['WT'][:]
    ds.variables['WA'][h-1,:] = ds_in.variables['WA'][:]
    ds.variables['SMCWTD'][h-1,:] = ds_in.variables['SMCWTD'][:]
    ds.variables['RECH'][h-1,:] = ds_in.variables['RECH'][:]
    ds.variables['SH2O'][h-1,:,:] = ds_in.variables['SH2O'][:]
    ds.variables['DEEPRECH'][h-1,:] = ds_in.variables['DEEPRECH'][:]
    ds.variables['SSTC'][h-1,:,:] = ds_in.variables['SSTC'][:]
    ds.variables['TG'][h-1, :] = ds_in.variables['TG'][:]
ds.close()
