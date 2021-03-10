## Load Modules
import numpy as np
import pandas as pd
import os
from netCDF4 import Dataset
from wrf import getvar, latlon_coords
import datetime as dt

f2015_d02 = '/scratch/leuven/336/vsc33651/nu-wrf-dev/wrfout_nc_files/wrfout_d01_2015.nc'
ds2015_d02 = Dataset(f2015_d02, 'r')
lats_2015_d02 = ds2015_d02.variables['lat'][:]
lons_2015_d02 = ds2015_d02.variables['lon'][:]
time = ds2015_d02.variables['time'][:]
time = pd.to_datetime('2000010100',format='%Y%m%d%H') + pd.to_timedelta(time, unit='h')
print(time)
times = time - pd.to_timedelta(7, unit='h')
print(times)
year = times.year
time_obs = np.where([(year == 2015)])
time_obs1 = times[time_obs[1]]
months = time_obs1.month
summer_2015 = np.where([(months > 5) & (months < 9)])
time_2015 = times[summer_2015[1]]
print(time_2015)
LTG3_MAX = ds2015_d02['LTG3_MAX'][summer_2015[1],:,:]
RAINNC = ds2015_d02.variables['RAINNC'][summer_2015[1],:,:]
RAINC = ds2015_d02.variables['RAINC'][summer_2015[1],:,:]
TOTAL_PREC = RAINNC + RAINC
TOTAL_PREC_2 = np.zeros((len(time_2015), 140,200))
TOTAL_PREC_2[0,:,:] = TOTAL_PREC[0,:,:]
for i in range(1,len(TOTAL_PREC)):
    TOTAL_PREC_2[i,:,:] = np.subtract(TOTAL_PREC[i,:,:],TOTAL_PREC[i-1,:,:])

# CREATE NEW NETCDF FILE
ds_diurnal = Dataset('/scratch/leuven/336/vsc33651/nu-wrf-dev/wrfout_nc_files/2015_d01_diurnal.nc', mode='w', format='NETCDF4')
ds_diurnal.createDimension('time', 24)
ds_diurnal.createDimension('lat', 140)
ds_diurnal.createDimension('lon', 200)
ds_diurnal.createVariable('time','int', dimensions=('time',),zlib=True)
ds_diurnal.variables['time'][:] = np.linspace(0,24, num=24, endpoint=False)
ds_diurnal.createVariable('lat',lats_2015_d02.dtype, dimensions=('lat','lon'),zlib=True)
ds_diurnal.variables['lat'][:] = lats_2015_d02[:]
ds_diurnal.createVariable('lon',lons_2015_d02.dtype, dimensions=('lat','lon'),zlib=True)
ds_diurnal.variables['lon'][:] = lons_2015_d02[:]
ds_diurnal.createVariable('LTG1_MAX','f4', dimensions=('time','lat','lon'),zlib=True)
ds_diurnal.createVariable('LTG2_MAX','f4', dimensions=('time','lat','lon'),zlib=True)
ds_diurnal.createVariable('LTG3_MAX','f4', dimensions=('time','lat','lon'),zlib=True)
ds_diurnal.createVariable('CTOP2D','f4', dimensions=('time','lat','lon'),zlib=True)
ds_diurnal.createVariable('COD2D','f4', dimensions=('time','lat','lon'),zlib=True)
ds_diurnal.createVariable('RAINNC','f4', dimensions=('time','lat','lon'),zlib=True)
ds_diurnal.createVariable('RAINC','f4', dimensions=('time','lat','lon'),zlib=True)
ds_diurnal.createVariable('TOTAL_PREC','f4', dimensions=('time','lat','lon'),zlib=True)
ds_diurnal.variables['time'].setncatts({'long_name': 'time', 'units': 'hour of the day'})
ds_diurnal.variables['lat'].setncatts({'long_name': 'latitude',  'units': 'degrees_north'})
ds_diurnal.variables['lon'].setncatts({'long_name': 'longitude', 'units': 'degrees_east'})
ds_diurnal.variables['LTG1_MAX'].setncatts({'long_name': 'Max lightning threat 1',  'units': 'flashes km$^{-2}$ (5min)$^{-1}$'})
ds_diurnal.variables['LTG2_MAX'].setncatts({'long_name': 'Max lightning threat 2',  'units': 'flashes km$^{-2}$ (5min)$^{-1}$'})
ds_diurnal.variables['LTG3_MAX'].setncatts({'long_name': 'Max lightning threat 3',  'units': 'flashes km$^{-2}$ (5min)$^{-1}$'})
ds_diurnal.variables['CTOP2D'].setncatts({'long_name': 'Cloud top pressure',  'units': 'mbar'})
ds_diurnal.variables['COD2D'].setncatts({'long_name': 'Column cloud optical depth',  'units': '-'})
ds_diurnal.variables['RAINNC'].setncatts({'long_name': 'precipitation',  'units': 'mm'})
ds_diurnal.variables['RAINC'].setncatts({'long_name': 'cumulus precipitation',  'units': 'mm'})
ds_diurnal.variables['TOTAL_PREC'].setncatts({'long_name': 'hourly precipitation',  'units': 'mm'})

for i in range(0,24):
    time_select_vector = np.where(time_2015.hour == i)
    ds_diurnal['LTG3_MAX'][i, :, :] = np.mean(LTG3_MAX[time_select_vector[0], :,:], axis=0)
    ds_diurnal['TOTAL_PREC'][i, :, :] = np.mean(TOTAL_PREC_2[time_select_vector[0], :,:], axis=0)
print(np.amax(TOTAL_PREC_2))
ds_diurnal.close()
