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
ds_d01_in = Dataset('/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/NU-WRF/wrfout_nc_files/regrid_72x72_Thompson.nc', 'r')
lat_d01 = ds_d01_in['lat'][:]
lon_d01 = ds_d01_in['lon'][:]

# ---------------------------------------------------------------------------------------------
# CREATE NC FILE
# ---------------------------------------------------------------------------------------------
print('creating nc file')
# D01
ds_d01 = Dataset('/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/NU-WRF/wrfout_nc_files/domain1_6h_72x72_Thompson.nc', mode='w', format='NETCDF4')
ds_d01.createDimension('time', 2208)
ds_d01.createDimension('lat', 8)
ds_d01.createDimension('lon', 14)
ds_d01.createVariable('time','int', dimensions=('time',),zlib=True)
ds_d01.createVariable('lat',lat_d01.dtype, dimensions=('lat','lon',),zlib=True)
ds_d01.variables['lat'][:] = lat_d01
ds_d01.createVariable('lon',lon_d01.dtype, dimensions=('lat','lon',),zlib=True)
ds_d01.variables['lon'][:] = lon_d01
ds_d01.createVariable('LPI_d01', 'f4', dimensions=('time','lat','lon',), zlib=True)
ds_d01.createVariable('LTG3_d01', 'f4', dimensions=('time','lat','lon',), zlib=True)
ds_d01.createVariable('PR92W_d01', 'f4', dimensions=('time','lat','lon',), zlib=True)
ds_d01.createVariable('CAPExP_R_d01', 'f4', dimensions=('time','lat','lon',), zlib=True)
ds_d01.createVariable('Obs', 'f4', dimensions=('time','lat','lon',), zlib=True)
ds_d01.createVariable('LPI_d02', 'f4', dimensions=('time','lat','lon',), zlib=True)
ds_d01.createVariable('LTG3_d02', 'f4', dimensions=('time','lat','lon',), zlib=True)
ds_d01.createVariable('PR92W_d02', 'f4', dimensions=('time','lat','lon',), zlib=True)
ds_d01.createVariable('CAPExP_R_d02', 'f4', dimensions=('time','lat','lon',), zlib=True)

# ---------------------------------------------------------------------------------------------
# FILL NC FILE
# ---------------------------------------------------------------------------------------------
ds_d01['time'][0:368] = pd.date_range('2015-06-01','2015-08-31 18:00:00', freq='6H')
ds_d01['time'][368:736] = pd.date_range('2016-06-01','2016-08-31 18:00:00', freq='6H')
ds_d01['time'][736:1104] = pd.date_range('2017-06-01','2017-08-31 18:00:00', freq='6H')
ds_d01['time'][1104:1472] = pd.date_range('2018-06-01','2018-08-31 18:00:00', freq='6H')
ds_d01['time'][1472:1840] = pd.date_range('2019-06-01','2019-08-31 18:00:00', freq='6H')
ds_d01['time'][1840:2208] = pd.date_range('2020-06-01','2020-08-31 18:00:00', freq='6H')

print('for loop')
days = np.divide(13248,6).astype(int)
for i in range(0,days):
    print('loop' +str(i))
    j = (i+1)*6
    k = i*6
    ds_d01['LPI_d02'][i,:,:] = np.nansum(ds_d01_in['LPI_d02'][k:j,:,:], axis = 0)
    ds_d01['LPI_d01'][i, :, :] = np.nansum(ds_d01_in['LPI_d01'][k:j, :, :], axis=0)
    ds_d01['LTG3_d02'][i, :, :] = np.nansum(ds_d01_in['LTG3_d02'][k:j,:,:], axis = 0)
    ds_d01['LTG3_d01'][i, :, :] = np.nansum(ds_d01_in['LTG3_d01'][k:j,:,:], axis = 0)
    ds_d01['PR92W_d02'][i, :, :] = np.nansum(ds_d01_in['PR92W_d02'][k:j,:,:], axis = 0)
    ds_d01['PR92W_d01'][i, :, :] = np.nansum(ds_d01_in['PR92W_d01'][k:j,:,:], axis = 0)
    ds_d01['CAPExP_R_d01'][i, :, :] = np.nansum(ds_d01_in['CAPExP_R_d01'][k:j, :, :], axis=0)
    ds_d01['CAPExP_R_d02'][i, :, :] = np.nansum(ds_d01_in['CAPExP_R_d02'][k:j, :, :], axis=0)
    ds_d01['Obs'][i, :, :] = np.nansum(ds_d01_in['Obs'][k:j, :, :], axis=0)

ds_d01.close()
