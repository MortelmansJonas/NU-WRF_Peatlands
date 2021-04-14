# ---------------------------------------------------------------------------------------------
# MODULES
# ---------------------------------------------------------------------------------------------
import netCDF4 as nc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from netCDF4 import Dataset
import os
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import LogNorm, from_levels_and_colors, ListedColormap, LinearSegmentedColormap
from matplotlib import cm
import matplotlib.dates as mdates
from datetime import timedelta
# ---------------------------------------------------------------------------------------------
# LOAD DATA
# ---------------------------------------------------------------------------------------------
ds_d01_in = Dataset('/scratch/leuven/336/vsc33651/nu-wrf-dev/domain1_weekly.nc', 'r')
ds_d02_in = Dataset('/scratch/leuven/336/vsc33651/nu-wrf-dev/domain2_weekly.nc', 'r')

# ---------------------------------------------------------------------------------------------
# CREATE NC FILE
# ---------------------------------------------------------------------------------------------
ds_d01 = Dataset('/scratch/leuven/336/vsc33651/nu-wrf-dev/domain1_seasonality.nc', mode='w', format='NETCDF4')
ds_d02 = Dataset('/scratch/leuven/336/vsc33651/nu-wrf-dev/domain2_seasonality.nc', mode='w', format='NETCDF4')

ds_d01.createDimension('time', 13) # Because 13 weeks
ds_d01.createVariable('time', 'int', dimensions=('time',),zlib=True)
ds_d01.createVariable('LPI', 'f4', dimensions=('time',), zlib=True)
ds_d01.createVariable('LTG3', 'f4', dimensions=('time',), zlib=True)
ds_d01.createVariable('PR92H', 'f4', dimensions=('time',), zlib=True)
ds_d01.createVariable('PR92W', 'f4', dimensions=('time',), zlib=True)
ds_d01.createVariable('CAPExP', 'f4', dimensions=('time',), zlib=True)
ds_d01.createVariable('Obs', 'f4', dimensions=('time',), zlib=True)

ds_d02.createDimension('time', 13) #} Because 13 weeks
ds_d02.createVariable('time', 'int', dimensions=('time',),zlib=True)
ds_d02.createVariable('LPI', 'f4', dimensions=('time',), zlib=True)
ds_d02.createVariable('LTG3', 'f4', dimensions=('time',), zlib=True)
ds_d02.createVariable('PR92H', 'f4', dimensions=('time',), zlib=True)
ds_d02.createVariable('PR92W', 'f4', dimensions=('time',), zlib=True)
ds_d02.createVariable('Obs', 'f4', dimensions=('time',), zlib=True)

ds_obs = Dataset('/scratch/leuven/336/vsc33651/nu-wrf-dev/wrfout_nc_files/Slave_lake_daily_0.1d.nc')

timeseries = pd.date_range('2015-05-15', '2020-09-15', freq='D') # Because this is the range in *weekly for the timeseries.
year = timeseries.year
month = timeseries.month

i2015 = np.where((year == 2015) & (month > 5) & (month < 9)) # To get all indices in timeseries that are in the summer months for each year (to cut the gaps)
i2016 = np.where((year == 2016) & (month > 5) & (month < 9))
i2017 = np.where((year == 2017) & (month > 5) & (month < 9))
i2018 = np.where((year == 2018) & (month > 5) & (month < 9))
i2019 = np.where((year == 2019) & (month > 5) & (month < 9))
i2020 = np.where((year == 2020) & (month > 5) & (month < 9))

# ---------------------------------------------------------------------------------------------
# OBSERVATIONS
# ---------------------------------------------------------------------------------------------
# Calculate weekly average
data_obs = np.zeros((len(timeseries), 2))
time_obs = pd.to_datetime(ds_obs['time'][:], format='%Y-%m-%d') + timedelta(minutes=2)
year_obs = time_obs.year
month_obs = time_obs.month
i2015_obs = np.where((year_obs == 2015) & (month_obs > 5) & (month_obs < 9)) # The same as before, but this time for the observations.
i2016_obs = np.where((year_obs == 2016) & (month_obs > 5) & (month_obs < 9))
i2017_obs = np.where((year_obs == 2017) & (month_obs > 5) & (month_obs < 9))
i2018_obs = np.where((year_obs == 2018) & (month_obs > 5) & (month_obs < 9))
i2019_obs = np.where((year_obs == 2019) & (month_obs > 5) & (month_obs < 9))
i2020_obs = np.where((year_obs == 2020) & (month_obs > 5) & (month_obs < 9))

summer = np.where([(month_obs > 5) & (month_obs < 9)]) # To only get the summer months across all years.
inds_obs_lat = np.unique(np.where((ds_obs['lat'][:] > 58) & (ds_obs['lat'][:] <62))[0])
inds_obs_lon = np.unique(np.where((ds_obs['lon'][:] > -123) & (ds_obs['lon'][:] <-108))[0])
CC_d01 = np.divide(np.nanmean(ds_obs['Flashdensity_CC'][:,inds_obs_lat,inds_obs_lon] + ds_obs['Flashdensity_CG'][:,inds_obs_lat,inds_obs_lon], axis=(1,2)), 61.94)

summer_int = np.where([(month > 5) & (month < 9)])
date = pd.to_datetime(time_obs.strftime('%Y-%m-%d'))
data_obs = np.zeros((len(timeseries),2))
data_obs[:,0] = pd.to_datetime(timeseries)
duplicates = np.in1d(date[:], timeseries[summer_int[1]]).astype(int) # To only get times from observations that match timeseries dates
duplicates_inv = np.in1d(timeseries[:], date[summer[1]]).astype(int) # and viceversa
data_obs[np.where(duplicates_inv == 1),1] = CC_d01[np.where(duplicates == 1)]

weekly_obs = np.zeros((546))
weekly_obs[0:91] = data_obs[17:108,1] # Values on the right (17:108) are the values of i2015
weekly_obs[91:182] = data_obs[383:474,1]
weekly_obs[182:273] = data_obs[748:839,1]
weekly_obs[273:364] = data_obs[1113:1204,1]
weekly_obs[364:455] = data_obs[1478:1569,1]
weekly_obs[455:546] = data_obs[1844:1935,1]

w_obs = np.zeros((78))
for i in range(0,78):
    j = (i+1)*7
    k=i*7
    w_obs[i] = np.nanmean(weekly_obs[k:j])
# ---------------------------------------------------------------------------------------------
# SEASONALITY 2.0
# ---------------------------------------------------------------------------------------------
ds_d01['LPI'][:] = 0
ds_d01['LTG3'][:] = 0
ds_d01['PR92H'][:] = 0
ds_d01['PR92W'][:] = 0

ds_d02['LPI'][:] = 0
ds_d02['LTG3'][:] = 0
ds_d02['PR92H'][:] = 0
ds_d02['PR92W'][:] = 0

for i in range(0,13):
    indices = np.arange(i,78,13)
    ds_d01['LPI'][i] = np.nanmean(ds_d01_in['LPI'][indices])
    ds_d02['LPI'][i] = np.nanmean(ds_d02_in['LPI'][indices])
    ds_d01['LTG3'][i] = np.nanmean(ds_d01_in['LTG3'][indices])
    ds_d02['LTG3'][i] = np.nanmean(ds_d02_in['LTG3'][indices])
    ds_d01['PR92H'][i] = np.nanmean(ds_d01_in['PR92H'][indices])
    ds_d02['PR92H'][i] = np.nanmean(ds_d02_in['PR92H'][indices])
    ds_d01['PR92W'][i] = np.nanmean(ds_d01_in['PR92W'][indices])
    ds_d02['PR92W'][i] = np.nanmean(ds_d02_in['PR92W'][indices])
    ds_d01['CAPExP'][i] = np.nanmean(ds_d01_in['CAPExP'][indices])
    ds_d01['Obs'][i] = np.nanmean(w_obs[indices])
    ds_d02['Obs'][i] = np.nanmean(w_obs[indices])

ds_d01.close()
ds_d02.close()
