#!/usr/bin/env python
# ---------------------------------------------------------------------------------------------
# MODULES
# ---------------------------------------------------------------------------------------------
import numpy as np
import pandas as pd
from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()
from netCDF4 import Dataset
import os

# ---------------------------------------------------------------------------------------------
# LIS OUTPUT
# ---------------------------------------------------------------------------------------------
date_from = '2015-06-01 00'
date_to = '2015-07-16'

# Load data
path = '/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc31786/Great_Slave_Lake/2015_thompson'
filename_root = 'LIS_HIST_'
domain = 'd02'
output_filename = os.path.join('/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/wrfout_nc_files',filename_root+domain+'_2015_Thompson.nc')

startdate = pd.to_datetime(date_from, format = '%Y-%m-%d %H')
enddate = pd.to_datetime(date_to, format = '%Y-%m-%d')
time_passed = (enddate - startdate)
nhours = 24*time_passed.days + time_passed.seconds//3600
hour = np.linspace(0,nhours,nhours//3).astype(int)
filename = os.path.join(path, filename_root+
                        str(startdate.year)+str(startdate.strftime('%m'))+
                        str(startdate.strftime('%d'))+str(startdate.strftime('%H'))+ '00'+ '.'+domain+'.nc')
ds_st = Dataset(filename, 'r')

# Create .nc file
ds = Dataset(output_filename, mode='w', format='NETCDF4')
print(ds_st)
timeunit = 'hours since 01-06-2015 00:00:00'
ds.createDimension('time', nhours//3)
ds.createDimension('lat', np.shape(ds_st['lat'][:])[0])
ds.createDimension('lon', np.shape(ds_st['lon'][:])[1])
ds.createVariable('time','int', dimensions=('time',),zlib=True)
ds.variables['time'][:] = hour[:]
ds.createVariable('lat',ds_st['lat'].dtype, dimensions=('lat','lon',),zlib=True)
ds.variables['lat'][:] = ds_st['lat'][:]
ds.createVariable('lon',ds_st['lon'].dtype, dimensions=('lat','lon',),zlib=True)
ds.variables['lon'][:] = ds_st['lon'][:]
ds.createVariable('Albedo','f4', dimensions=('time','lat','lon'),zlib=True)
ds.createVariable('GF','f4', dimensions=('time','lat','lon'),zlib=True)

ds.variables['time'].setncatts({'long_name': 'time', 'units': timeunit})
ds.variables['lat'].setncatts({'long_name': 'latitude',  'units': 'degrees_north'})
ds.variables['lon'].setncatts({'long_name': 'longitude', 'units': 'degrees_east'})
ds.variables['Albedo'].setncatts({'long_name': 'Surface Albedo',  'units': '-'})
ds.variables['GF'].setncatts({'long_name': 'Green Vegetation Fraction',  'units': '-'})

# Fill .nc file
for h in range(0,nhours,3): # in steps of 3 because 3-hourly output
    hours_passed = pd.to_timedelta(h, unit='h')
    a = startdate + hours_passed

    filename = os.path.join(path, filename_root + str(2015) + str(a.strftime('%m')) +
                            str(a.strftime('%d')) + str(a.strftime('%H')) + '00' + '.' + domain + '.nc')
    print('processing '+filename)
    ds_in = Dataset(filename, mode='r')

    ds.variables['Albedo'][h/3,:,:] = ds_in.variables['Albedo_inst'][:]
    ds.variables['GF'][h/3,:,:] = ds_in.variables['Greenness_inst'][:]
ds.close()
