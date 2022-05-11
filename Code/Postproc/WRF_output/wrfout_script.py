# ---------------------------------------------------------------------------------------------
# MODULES
# ---------------------------------------------------------------------------------------------
import numpy as np
import pandas as pd
import os
from netCDF4 import Dataset
import xarray as xr
from wrf import (getvar, interplevel, vertcross,
                 CoordPair, ALL_TIMES, to_np,
                 get_cartopy, latlon_coords,
                 cartopy_xlim, cartopy_ylim)

# ---------------------------------------------------------------------------------------------
# LOAD DATA
# ---------------------------------------------------------------------------------------------
date_from = '2015-06-01'
date_to = '2015-09-01'

path = '/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc31786/Great_Slave_Lake/2015'
# path = '/scratch/leuven/336/vsc33651/nu-wrf-dev/Great_Slave_Lake/2015_lake'
filename_root = 'wrfout_'
domain = 'd01'
output_filename = os.path.join('/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/wrfout_nc_files',filename_root+domain+'_2015.nc')
# output_filename = os.path.join('/scratch/leuven/336/vsc33651/nu-wrf-dev/wrfout_nc_files',filename_root+domain+'_2015_lake.nc')

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
ds.createVariable('LTG1_MAX','f4', dimensions=('time','lat','lon'),zlib=True)
ds.createVariable('LTG2_MAX','f4', dimensions=('time','lat','lon'),zlib=True)
ds.createVariable('LTG3_MAX','f4', dimensions=('time','lat','lon'),zlib=True)
ds.createVariable('LTG3_MAX_new','f4', dimensions=('time','lat','lon'),zlib=True)
ds.createVariable('CTOP2D','f4', dimensions=('time','lat','lon'),zlib=True)
ds.createVariable('COD2D','f4', dimensions=('time','lat','lon'),zlib=True)
ds.createVariable('RAINC','f4', dimensions=('time','lat','lon'),zlib=True)
ds.createVariable('RAINNC','f4', dimensions=('time','lat','lon'),zlib=True)
ds.createVariable('W_UP_MAX','f4', dimensions=('time','lat','lon'),zlib=True)
ds.createVariable('HGT','f4', dimensions=('time','lat','lon'),zlib=True)
ds.createVariable('CAPE','f4', dimensions=('time','lat','lon'),zlib=True)
ds.createVariable('CSI','f4', dimensions=('time','lat','lon'),zlib=True)
ds.createVariable('LH','f4', dimensions=('time','lat','lon'),zlib=True)
ds.createVariable('SH','f4', dimensions=('time','lat','lon'),zlib=True)
ds.createVariable('RH2','f4', dimensions=('time','lat','lon'),zlib=True)
ds.createVariable('T2','f4', dimensions=('time','lat','lon'),zlib=True)
ds.createVariable('DBZ','f4', dimensions=('time','lat','lon'),zlib=True)

ds.variables['time'].setncatts({'long_name': 'time',  'units': timeunit})
ds.variables['lat'].setncatts({'long_name': 'latitude',  'units': 'degrees_north'})
ds.variables['lon'].setncatts({'long_name': 'longitude', 'units': 'degrees_east'})
ds.variables['LTG1_MAX'].setncatts({'long_name': 'Max lightning threat 1',  'units': 'flashes km$^{-2}$ (5min)$^{-1}$'})
ds.variables['LTG2_MAX'].setncatts({'long_name': 'Max lightning threat 2',  'units': 'flashes km$^{-2}$ (5min)$^{-1}$'})
ds.variables['LTG3_MAX'].setncatts({'long_name': 'Max lightning threat 3',  'units': 'flashes km$^{-2}$ (5min)$^{-1}$'})
ds.variables['LTG3_MAX_new'].setncatts({'long_name': 'Max lightning threat 3',  'units': 'flashes km$^{-2}$ (5min)$^{-1}$'})
ds.variables['CTOP2D'].setncatts({'long_name': 'Cloud top pressure',  'units': 'mbar'})
ds.variables['COD2D'].setncatts({'long_name': 'Column cloud optical depth',  'units': '-'})
ds.variables['RAINC'].setncatts({'long_name': 'Cumulus precipitation',  'units': 'mm'})
ds.variables['RAINNC'].setncatts({'long_name': 'precipitation',  'units': 'mm'})
ds.variables['W_UP_MAX'].setncatts({'long_name': 'maximal updraft vertical velocity',  'units': 'm/s'})
ds.variables['HGT'].setncatts({'long_name': 'elevation',  'units': 'm'})
ds.variables['CAPE'].setncatts({'long_name': 'Convective Available Potential Energy',  'units': 'J/kg'})
ds.variables['LH'].setncatts({'long_name': 'Latent heat flux',  'units': 'W/m²'})
ds.variables['SH'].setncatts({'long_name': 'Sensible heat flux',  'units': 'W/m²'})
ds.variables['RH2'].setncatts({'long_name': '2-meter relative humidity',  'units': '%'})
ds.variables['T2'].setncatts({'long_name': '2-meter temperature',  'units': 'Degrees Celsius'})
ds.variables['DBZ'].setncatts({'long_name': 'Reflectivity',  'units': 'dBZ'})

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

    ds.variables['LTG1_MAX'][h,:,:] = ds_in.variables['LTG1_MAX'][:]
    ds.variables['LTG2_MAX'][h,:,:] = ds_in.variables['LTG2_MAX'][:]
    ds.variables['LTG3_MAX'][h,:,:] = 0.95 * ds_in.variables['LTG1_MAX'][:] + 0.05 * ds_in.variables['LTG2_MAX'][:]
    ds.variables['LTG3_MAX_new'][h, :, :] = 0.80 * ds_in.variables['LTG1_MAX'][:] + 0.20 * ds_in.variables['LTG2_MAX'][:]
    ds.variables['CTOP2D'][h, :, :] = ds_in.variables['CTOP2D_OUT'][:]
    ds.variables['COD2D'][h, :, :] = ds_in.variables['COD2D_OUT'][:]
    ds.variables['RAINC'][h, :, :] = ds_in.variables['RAINC'][:]
    ds.variables['RAINNC'][h, :, :] = ds_in.variables['RAINNC'][:]
    ds.variables['W_UP_MAX'][h, :, :] = ds_in.variables['W_UP_MAX'][:]
    ds.variables['HGT'][h, :, :] = ds_in.variables['HGT'][:]
    ds.variables['CAPE'][h, :, :] = getvar(ds_in, 'cape_2d')[0,:,:]
    ds.variables['DBZ'][h, :, :] = getvar(ds_in, 'dbz')[:, :]
    ds.variables['CSI'][h, :, :] = ds_in.variables['CSI'][:]
    ds.variables['LH'][h, :, :] = ds_in.variables['LH'][:]
    ds.variables['SH'][h, :, :] = ds_in.variables['HFX'][:]
    ds.variables['RH2'][h, :, :] = getvar(ds_in, 'rh2')[:, :]
    ds.variables['T2'][h, :, :] = getvar(ds_in, 'T2')[:, :]

ds.close()