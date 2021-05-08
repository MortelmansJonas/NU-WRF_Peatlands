# ---------------------------------------------------------------------------------------------
# MODULES
# ---------------------------------------------------------------------------------------------
import numpy as np
import pandas as pd
import os
from netCDF4 import Dataset
from wrf import getvar, latlon_coords

# ---------------------------------------------------------------------------------------------
# LOAD DATA
# ---------------------------------------------------------------------------------------------
date_from = '2015-06-01'
date_to = '2015-09-01'

path = '/scratch/leuven/336/vsc33651/nu-wrf-dev/Great_Slave_Lake/2015'
filename_root = 'wrfout_'
domain = 'd01'
output_filename = os.path.join('/scratch/leuven/336/vsc33651/nu-wrf-dev/wrfout_nc_files',filename_root+domain+'_2015_MERRA2.nc')

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
ds = Dataset(output_filename, mode='w', format='NETCDF4')
print(ds_st)
timeunit = 'hours since 2000-01-01 00:00'
ds.createDimension('time', None)
ds.createDimension('lat', np.shape(lat)[0])
ds.createDimension('lon', np.shape(lon)[1])
ds.createVariable('time',hours.dtype, dimensions=('time'),zlib=True)
ds.variables['time'][:] = hours + (pd.to_datetime(date_from) - pd.to_datetime('2000-01-01')).days*24 - 1
ds.createVariable('lat',lat.dtype, dimensions=('lat','lon'),zlib=True)
ds.variables['lat'][:] = lat
ds.createVariable('lon',lon.dtype, dimensions=('lat','lon'),zlib=True)
ds.variables['lon'][:] = lon
ds.createVariable('RAINC','f4', dimensions=('time','lat','lon'),zlib=True)
ds.createVariable('RAINNC','f4', dimensions=('time','lat','lon'),zlib=True)
ds.createVariable('RAINSH','f4', dimensions=('time','lat','lon'),zlib=True)
ds.createVariable('LH','f4', dimensions=('time','lat','lon'),zlib=True)
ds.createVariable('HFX','f4', dimensions=('time','lat','lon'),zlib=True)

ds.variables['time'].setncatts({'long_name': 'time',  'units': timeunit})
ds.variables['lat'].setncatts({'long_name': 'latitude',  'units': 'degrees_north'})
ds.variables['lon'].setncatts({'long_name': 'longitude', 'units': 'degrees_east'})
ds.variables['RAINC'].setncatts({'long_name': 'Cumulus precipitation',  'units': 'mm'})
ds.variables['RAINNC'].setncatts({'long_name': 'precipitation',  'units': 'mm'})
ds.variables['RAINSH'].setncatts({'long_name': 'Shallow Cumulus precipitation',  'units': 'mm'})
ds.variables['LH'].setncatts({'long_name': 'Latent heat flux',  'units': 'W/m²'})
ds.variables['HFX'].setncatts({'long_name': 'Sensible heat flux',  'units': 'W/m²'})

hr = pd.to_timedelta(startdate - pd.to_datetime('2015-06-01')).days*24 + \
     pd.to_timedelta(startdate - pd.to_datetime('2015-06-01')).seconds//3600
for h in range(hr,hr+nhours+1):
    hours_passed = pd.to_timedelta(h, unit='h')
    a = pd.to_datetime('2015-06-01') + hours_passed
    i = str(a)
    day = i.split(' ')
    filename = os.path.join(path, filename_root + domain + '_' + day[0] + '_' + day[1])
    print('processing '+filename)
    ds_in = Dataset(filename, mode='r')

    ds.variables['LH'][h,:,:] = ds_in.variables['LH'][:]
    ds.variables['HFX'][h, :, :] = ds_in.variables['HFX'][:]
    ds.variables['RAINC'][h, :, :] = ds_in.variables['RAINC'][:]
    ds.variables['RAINNC'][h, :, :] = ds_in.variables['RAINNC'][:]
    ds.variables['RAINSH'][h, :, :] = ds_in.variables['W_UP_MAX'][:]
ds.close()
