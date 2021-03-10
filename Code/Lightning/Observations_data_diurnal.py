## Load Modules
import numpy as np
import pandas as pd
import os
from netCDF4 import Dataset
from wrf import getvar, latlon_coords
from datetime import timedelta
import datetime as dt
import pytz

nf= '/data/leuven/336/vsc33651/projects/nu-wrf-dev/Lightning_data/Slave_lake_hourly.nc'
ds_hourly = Dataset(nf, 'r')
lats = ds_hourly.variables['lat'][:]
lons = ds_hourly.variables['lon'][:]

time = pd.to_datetime(ds_hourly.variables['time'][:]) + timedelta(minutes=2) # Add 2 minutes again to get correct times
time_UTC = time.tz_localize('UTC')
times = time_UTC.tz_convert('MST')
year = times.year
time_obs = np.where([(year == 2015)])
time_obs1 = times[time_obs[1]]
Flashdensity_CC_2015 = ds_hourly.variables['Flashdensity_CC'][time_obs[1],:,:]
Flashdensity_CG_2015 = ds_hourly.variables['Flashdensity_CG'][time_obs[1],:,:]
months = time_obs1.month
summer_2015 = np.where([(months > 5) & (months < 9)])
time_2015 = times[summer_2015[1]]
time_2015_inds = np.where([(time_2015.month == 8) & (time_2015.day == 31) & (time_2015.hour >= 17)])
time_2015 = time_2015[0:1238]
print(time_2015)
F_CC_summer_2015 = Flashdensity_CC_2015[summer_2015[1],:,:]
F_CC_summer_2015 = np.where(np.isnan(F_CC_summer_2015),0, F_CC_summer_2015)
F_CG_summer_2015 = Flashdensity_CG_2015[summer_2015[1],:,:]
F_CG_summer_2015 = np.where(np.isnan(F_CG_summer_2015),0, F_CG_summer_2015)
total_flashes_2015 = np.add(F_CG_summer_2015, F_CG_summer_2015)
print(np.amax(total_flashes_2015[np.where(time_2015.hour ==8),:,:]))

## Create netcdf file to store diurnal data (length of time vector = 24)
ds_obs = Dataset('/scratch/leuven/336/vsc33651/nu-wrf-dev/Observations_diurnal.nc', mode='w', format='NETCDF4')
ds_obs.createDimension('time', 24)
ds_obs.createDimension('lat', len(lats))
ds_obs.createDimension('lon', len(lons))
ds_obs.createVariable('time','int', dimensions=('time',),zlib=True)
ds_obs.variables['time'][:] = np.linspace(0,24, num=24, endpoint=False)
ds_obs.createVariable('lat',lats.dtype, dimensions=('lat',),zlib=True)
ds_obs.variables['lat'][:] = lats[:]
ds_obs.createVariable('lon',lons.dtype, dimensions=('lon',),zlib=True)
ds_obs.variables['lon'][:] = lons[:]
ds_obs.createVariable('LF_diurnal','f4', dimensions=('time','lat','lon'),zlib=True)
ds_obs.variables['time'].setncatts({'long_name': 'time', 'units': 'hour of the day'})
ds_obs.variables['lat'].setncatts({'long_name': 'latitude',  'units': 'degrees_north'})
ds_obs.variables['lon'].setncatts({'long_name': 'longitude', 'units': 'degrees_east'})
ds_obs.variables['LF_diurnal'].setncatts({'long_name': 'average number of lightning flashes per hour of the day',
                                            'units': 'flashes (100(km)$^{2}$ hr)$^{-1}$'})

for i in range(0,24):
    time_select_vector = np.where(time_2015.hour == i)
    ds_obs['LF_diurnal'][i, :, :] = np.nanmean(total_flashes_2015[time_select_vector, :, :], axis=0)
    print(np.amax(ds_obs['LF_diurnal'][i,:,:]))
ds_obs.close()