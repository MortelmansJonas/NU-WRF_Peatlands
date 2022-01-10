#!/usr/bin/env python
# ---------------------------------------------------------------------------------------------
# MODULES
# ---------------------------------------------------------------------------------------------
import numpy as np
import pandas as pd
from netCDF4 import Dataset

# ---------------------------------------------------------------------------------------------
# LOAD DATA
# ---------------------------------------------------------------------------------------------
print('load data')
ds_d01_in = Dataset('/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/wrfout_nc_files/data_calibrated_ax_d01.nc', mode='r')

lat_d01 = ds_d01_in['lat'][:]
lon_d01 = ds_d01_in['lon'][:]

# ---------------------------------------------------------------------------------------------
# CREATE NC FILE
# ---------------------------------------------------------------------------------------------
print('creating nc file')
# D01
ds_d01 = Dataset('/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/wrfout_nc_files/domain1_daily_regridded.nc', mode='w', format='NETCDF4')
ds_d01.createDimension('time', 552)
ds_d01.createDimension('lat', 63)
ds_d01.createDimension('lon', 109)
ds_d01.createVariable('time','int', dimensions=('time',),zlib=True)
ds_d01.createVariable('lat',lat_d01.dtype, dimensions=('lat','lon',),zlib=True)
ds_d01.variables['lat'][:] = lat_d01
ds_d01.createVariable('lon',lon_d01.dtype, dimensions=('lat','lon',),zlib=True)
ds_d01.variables['lon'][:] = lon_d01
ds_d01.createVariable('LPI', 'f4', dimensions=('time','lat','lon',), zlib=True)
ds_d01.createVariable('LTG3', 'f4', dimensions=('time','lat','lon',), zlib=True)
ds_d01.createVariable('PR92W', 'f4', dimensions=('time','lat','lon',), zlib=True)
ds_d01.createVariable('CAPExP_R', 'f4', dimensions=('time','lat','lon',), zlib=True)
ds_d01.createVariable('CAPExP_CSI', 'f4', dimensions=('time','lat','lon',), zlib=True)
ds_d01.createVariable('Obs', 'f4', dimensions=('time','lat','lon',), zlib=True)
# D02
ds_d02 = Dataset('/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/wrfout_nc_files/domain2_daily_regridded.nc', mode='w', format='NETCDF4')
ds_d02.createDimension('time', 552)
ds_d02.createDimension('lat', 63)
ds_d02.createDimension('lon', 109)
ds_d02.createVariable('time','int', dimensions=('time',),zlib=True)
ds_d02.createVariable('lat',lat_d01.dtype, dimensions=('lat','lon',),zlib=True)
ds_d02.variables['lat'][:] = lat_d01
ds_d02.createVariable('lon',lon_d01.dtype, dimensions=('lat','lon',),zlib=True)
ds_d02.variables['lon'][:] = lon_d01
ds_d02.createVariable('LPI', 'f4', dimensions=('time','lat','lon',), zlib=True)
ds_d02.createVariable('LTG3', 'f4', dimensions=('time','lat','lon',), zlib=True)
ds_d02.createVariable('PR92W', 'f4', dimensions=('time','lat','lon',), zlib=True)
ds_d02.createVariable('CAPExP_R', 'f4', dimensions=('time','lat','lon',), zlib=True)
ds_d02.createVariable('CAPExP_CSI', 'f4', dimensions=('time','lat','lon',), zlib=True)
ds_d02.createVariable('Obs', 'f4', dimensions=('time','lat','lon',), zlib=True)

# ---------------------------------------------------------------------------------------------
# FILL NC FILE
# ---------------------------------------------------------------------------------------------
ds_d02['time'][0:92] = pd.date_range('2015-06-01','2015-08-31', freq='D')
ds_d02['time'][92:184] = pd.date_range('2016-06-01','2016-08-31', freq='D')
ds_d02['time'][184:276] = pd.date_range('2017-06-01','2017-08-31', freq='D')
ds_d02['time'][276:368] = pd.date_range('2018-06-01','2018-08-31', freq='D')
ds_d02['time'][368:460] = pd.date_range('2019-06-01','2019-08-31', freq='D')
ds_d02['time'][460:552] = pd.date_range('2020-06-01','2020-08-31', freq='D')
ds_d01['time'][0:92] = pd.date_range('2015-06-01','2015-08-31', freq='D')
ds_d01['time'][92:184] = pd.date_range('2016-06-01','2016-08-31', freq='D')
ds_d01['time'][184:276] = pd.date_range('2017-06-01','2017-08-31', freq='D')
ds_d01['time'][276:368] = pd.date_range('2018-06-01','2018-08-31', freq='D')
ds_d01['time'][368:460] = pd.date_range('2019-06-01','2019-08-31', freq='D')
ds_d01['time'][460:552] = pd.date_range('2020-06-01','2020-08-31', freq='D')

print('for loop')
days = np.divide(13248,24).astype(int)
for i in range(0,days):
    print('loop' +str(i))
    j = (i+1)*24
    k = i*24
    ds_d02['LPI'][i,:,:] = np.nansum(ds_d01_in['LPI_d02'][k:j,:,:], axis = 0)
    ds_d01['LPI'][i, :, :] = np.nansum(ds_d01_in['LPI_d01'][k:j, :, :], axis=0)
    ds_d02['LTG3'][i, :, :] = np.nansum(ds_d01_in['LTG3_d02'][k:j,:,:], axis = 0)
    ds_d01['LTG3'][i, :, :] = np.nansum(ds_d01_in['LTG3_d01'][k:j,:,:], axis = 0)
    ds_d02['PR92W'][i, :, :] = np.nansum(ds_d01_in['PR92W_d02'][k:j,:,:], axis = 0)
    ds_d01['PR92W'][i, :, :] = np.nansum(ds_d01_in['PR92W_d01'][k:j,:,:], axis = 0)
    ds_d01['CAPExP_R'][i, :, :] = np.nansum(ds_d01_in['CAPExP_R_d01'][k:j, :, :], axis=0)
    ds_d01['CAPExP_CSI'][i, :, :] = np.nansum(ds_d01_in['CAPExP_CSI_d01'][k:j, :, :], axis=0)
    ds_d02['CAPExP_R'][i, :, :] = np.nansum(ds_d01_in['CAPExP_R_d02'][k:j, :, :], axis=0)
    ds_d02['CAPExP_CSI'][i, :, :] = np.nansum(ds_d01_in['CAPExP_CSI_d02'][k:j, :, :], axis=0)
    ds_d01['Obs'][i, :, :] = np.nansum(ds_d01['Obs'][k:j, :, :], axis=0)
    ds_d02['Obs'][i, :, :] = np.nansum(ds_d01['Obs'][k:j, :, :], axis=0)

ds_d01.close()
ds_d02.close()
