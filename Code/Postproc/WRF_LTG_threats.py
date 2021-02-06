## Load Modules
import numpy as np
import pandas as pd
import os
from netCDF4 import Dataset
from wrf import getvar, latlon_coords

## Obtain the WRF Output files
date_from = '2015-06-01'
date_to = '2015-09-01'

path = '/scratch/leuven/336/vsc33651/nu-wrf-dev/Great_Slave_Lake/2015/'
filename_root = 'wrfout_'
domain = 'd02'
output_filename = os.path.join('/scratch/leuven/336/vsc33651/nu-wrf-dev/wrfout_nc_files/2015/',filename_root+domain+'.nc')

# determine length of time vector
startdate = pd.to_datetime(date_from)
enddate = pd.to_datetime(date_to)
ndays = (pd.to_datetime(date_to) - pd.to_datetime(date_from)).days
nhours = ndays * 24
hours = np.linspace(1,nhours,nhours)
day = str(startdate).split(' ')
filename = os.path.join(path, filename_root+domain + '_'+day[0]+'_'+day[1])
ds_st = Dataset(filename, 'r')
ter = getvar(ds_st, "ter")
lat, lon = latlon_coords(ter)
layers = np.arange(0,60,1, dtype = 'int')
ds = Dataset(output_filename, mode='w', format='NETCDF4')
timeunit = 'hours since 2000-01-01 00:00'
ds.createDimension('time', nhours)
ds.createDimension('lat', np.shape(lat)[0])
ds.createDimension('lon', np.shape(lon)[1])
ds.createDimension('layer', np.shape(layers)[0])
ds.createVariable('time',hours.dtype, dimensions=('time',),zlib=True)
ds.variables['time'][:] = hours + (pd.to_datetime(date_from) - pd.to_datetime('2000-01-01')).days*24 - 1
ds.createVariable('lat',lat.dtype, dimensions=('lat','lon'),zlib=True)
ds.variables['lat'][:] = lat
ds.createVariable('lon',lon.dtype, dimensions=('lat','lon'),zlib=True)
ds.variables['lon'][:] = lon
ds.createVariable('layer',layers.dtype, dimensions=('layer'),zlib=True)
ds.variables['layer'][:] = layers
ds.createVariable('QGRAUP','f4', dimensions=('time','layer','lat','lon'),zlib=True)
ds.createVariable('QSNOW','f4', dimensions=('time','layer','lat','lon'),zlib=True)
ds.createVariable('QICE','f4', dimensions=('time','layer','lat','lon'),zlib=True)
ds.createVariable('LTG1_MAX','f4', dimensions=('time','lat','lon'),zlib=True)
ds.createVariable('LTG2_MAX','f4', dimensions=('time','lat','lon'),zlib=True)
ds.createVariable('LTG3_MAX','f4', dimensions=('time','lat','lon'),zlib=True)

ds.variables['time'].setncatts({'long_name': 'time', 'units': timeunit})
ds.variables['lat'].setncatts({'long_name': 'latitude',  'units': 'degrees_north'})
ds.variables['lon'].setncatts({'long_name': 'longitude', 'units': 'degrees_east'})
ds.variables['layer'].setncatts({'long_name': 'vertical layer', 'units': '-'})
ds.variables['QGRAUP'].setncatts({'long_name': 'Graupel mixing ratio',  'units': 'kg (kg)$^{1}$'})
ds.variables['QICE'].setncatts({'long_name': 'Ice mixing ratio',  'units': 'kg (kg)$^{1}$'})
ds.variables['QSNOW'].setncatts({'long_name': 'Snow mixing ratio',  'units': 'kg (kg)$^{1}$'})
ds.variables['LTG1_MAX'].setncatts({'long_name': 'Max lightning threat 1',  'units': 'flashes km$^{-2}$ (5min)$^{-1}$'})
ds.variables['LTG2_MAX'].setncatts({'long_name': 'Max lightning threat 2',  'units': 'flashes km$^{-2}$ (5min)$^{-1}$'})
ds.variables['LTG3_MAX'].setncatts({'long_name': 'Max lightning threat 3',  'units': 'flashes km$^{-2}$ (5min)$^{-1}$'})

for h in range(1,nhours+1):
    hours_passed = pd.to_timedelta(h-1, unit='h')
    a = startdate + hours_passed
    i = str(a)
    day = i.split(' ')
    filename = os.path.join(path, filename_root + domain + '_' + day[0] + '_' + day[1])
    print('processing '+filename)
    ds_in = Dataset(filename, mode='r')
    ## Get LTG threat 1 and 2 from the WRF files
    qgraup = ds_in.variables['QGRAUP'][:]
    qsnow = ds_in.variables['QSNOW'][:]
    qice = ds_in.variables['QICE'][:]
    ltg1_max = ds_in.variables['LTG1_MAX'][:]
    ltg2_max = ds_in.variables['LTG2_MAX'][:]
    weight_ltg1 = np.multiply(ltg1_max, 0.95)
    weight_ltg2 = np.multiply(ltg2_max, 0.05)
    ltg3_max = np.add(ltg1_max, ltg2_max)

    ds.variables['QGRAUP'][h-1,:,:,:] = qgraup
    ds.variables['QICE'][h-1,:,:,:] = qice
    ds.variables['QSNOW'][h-1,:,:,:] = qsnow
    ds.variables['LTG1_MAX'][h-1,:,:] = ltg1_max
    ds.variables['LTG2_MAX'][h-1,:,:] = ltg2_max
    ds.variables['LTG3_MAX'][h-1,:,:] = ltg3_max
ds.close()