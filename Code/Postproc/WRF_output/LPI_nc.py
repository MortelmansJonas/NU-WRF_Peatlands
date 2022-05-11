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
date_from = '2020-06-01'
date_to = '2020-09-01'

path = '/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc31786/Great_Slave_Lake/2020_thompson'
filename_root = 'lpi_'
domain = 'd02'
output_filename = os.path.join('/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/wrfout_nc_files',filename_root+domain+'_2020_Thompson.nc')

# determine length of time vector
startdate = pd.to_datetime(date_from)
enddate = pd.to_datetime(date_to)
ndays = (pd.to_datetime(date_to) - pd.to_datetime(date_from)).days
nhours = ndays * 24
hours = np.linspace(1,nhours,nhours)
day = str(startdate).split(' ')
filename = os.path.join(path, 'wrfout_d02_2020-08-10_05:00:00')
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
ds.createVariable('LPI','f4', dimensions=('time','lat','lon'),zlib=True)

ds.variables['time'].setncatts({'long_name': 'time',  'units': timeunit})
ds.variables['lat'].setncatts({'long_name': 'latitude',  'units': 'degrees_north'})
ds.variables['lon'].setncatts({'long_name': 'longitude', 'units': 'degrees_east'})
ds.variables['LPI'].setncatts({'long_name': 'Lightning Potential Index',  'units': 'J kg$^{-1}$'})

hr = pd.to_timedelta(startdate - pd.to_datetime('2020-06-01')).days*24 + \
     pd.to_timedelta(startdate - pd.to_datetime('2020-06-01')).seconds//3600
for h in range(hr,hr+nhours+1):
    filename = os.path.join(path, filename_root + 'time_' + str(h+408) + '_' + domain + '.nc')
    print('processing '+filename)
    ds_in = Dataset(filename, mode='r')

    ds.variables['LPI'][h,:,:] = ds_in.variables['lpi'][:]

print(np.nanmax(ds['LPI'][:]))
ds.close()
