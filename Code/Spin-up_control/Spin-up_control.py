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

path = '/scratch/leuven/336/vsc33651/LIS/v03'
filename_root = 'LIS_HIST_'
domain = 'd02'
output_filename = os.path.join('/scratch/leuven/336/vsc33651/nu-wrf-dev/lis_spinup',filename_root+domain+'_v03.nc')

startdate = pd.to_datetime(date_from, format = '%Y-%m-%d %H')
enddate = pd.to_datetime(date_to, format = '%Y-%m-%d')
time_passed = (enddate - startdate)
nhours = 24*time_passed.days + time_passed.seconds//3600
hour = np.linspace(0,nhours,nhours//3).astype(int)
filename = os.path.join(path, filename_root+
                        str(startdate.year)+str(startdate.strftime('%m'))+
                        str(startdate.strftime('%d'))+str(startdate.strftime('%H'))+ '00'+ '.'+domain+'.nc')
ds_st = Dataset(filename, 'r')
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
ds.createVariable('Qle','f4', dimensions=('time','lat','lon'),zlib=True)
ds.createVariable('Qh','f4', dimensions=('time','lat','lon'),zlib=True)

ds.variables['time'].setncatts({'long_name': 'time', 'units': timeunit})
ds.variables['lat'].setncatts({'long_name': 'latitude',  'units': 'degrees_north'})
ds.variables['lon'].setncatts({'long_name': 'longitude', 'units': 'degrees_east'})
ds.variables['Qle'].setncatts({'long_name': 'Latent heat flux',  'units': 'W m$^{-2}$'})
ds.variables['Qh'].setncatts({'long_name': 'sensible heat flux',  'units': 'W m$^{-2}$'})

for h in range(0,nhours,3):
    hours_passed = pd.to_timedelta(h, unit='h')
    a = startdate + hours_passed

    filename = os.path.join(path, filename_root + str(2015) + str(a.strftime('%m')) +
                            str(a.strftime('%d')) + str(a.strftime('%H')) + '00' + '.' + domain + '.nc')
    print('processing '+filename)
    ds_in = Dataset(filename, mode='r')

    ds.variables['Qle'][h/3,:,:] = ds_in.variables['Qle_tavg'][:]
    ds.variables['Qh'][h/3,:,:] = ds_in.variables['Qh_tavg'][:]
ds.close()
