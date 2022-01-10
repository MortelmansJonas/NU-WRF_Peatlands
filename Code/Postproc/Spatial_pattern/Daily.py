# ---------------------------------------------------------------------------------------------
# MODULES
# ---------------------------------------------------------------------------------------------
import numpy as np
import pandas as pd
from netCDF4 import Dataset

# ---------------------------------------------------------------------------------------------
# LOAD DATA OF D02
# ---------------------------------------------------------------------------------------------
print('load data')
ds_d01_in = Dataset('/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/wrfout_nc_files/data_calibrated_ax_d01.nc', mode='r')
ds_d02_in = Dataset('/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/wrfout_nc_files/data_calibrated_ax_d02.nc', mode='r')
ds_obs_d01 = Dataset('//scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/wrfout_nc_files/CLDN_at_domain1_all_v4.nc','r')
ds_obs_d02 = Dataset('/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/wrfout_nc_files/CLDN_at_domain2_all.nc','r')

lat_d01 = ds_d01_in['lat'][:]
lat_d02 = ds_d02_in['lat'][:]
lon_d01 = ds_d01_in['lon'][:]
lon_d02 = ds_d02_in['lon'][:]
L_d01 = ds_obs_d01['Flashdensity_CC'][:].data + ds_obs_d01['Flashdensity_CG'][:].data
L_d02 = ds_obs_d02['Flashdensity_CC'][:].data + ds_obs_d02['Flashdensity_CG'][:].data

# ---------------------------------------------------------------------------------------------
# CREATE NC FILE
# ---------------------------------------------------------------------------------------------
print('creating nc file')
# D01
ds_d01 = Dataset('/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/wrfout_nc_files/domain1_daily.nc', mode='w', format='NETCDF4')
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
ds_d02 = Dataset('/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/wrfout_nc_files/domain2_daily.nc', mode='w', format='NETCDF4')
ds_d02.createDimension('time', 552)
ds_d02.createDimension('lat', 170)
ds_d02.createDimension('lon', 309)
ds_d02.createVariable('time','int', dimensions=('time',),zlib=True)
ds_d02.createVariable('lat',lat_d02.dtype, dimensions=('lat','lon',),zlib=True)
ds_d02.variables['lat'][:] = lat_d02
ds_d02.createVariable('lon',lon_d02.dtype, dimensions=('lat','lon',),zlib=True)
ds_d02.variables['lon'][:] = lon_d02
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
    ds_d02['LPI'][i,:,:] = np.nansum(ds_d02_in['LPI'][k:j,:,:], axis = 0)
    ds_d01['LPI'][i, :, :] = np.nansum(ds_d01_in['LPI_d01'][k:j, :, :], axis=0)
    ds_d02['LTG3'][i, :, :] = np.nansum(ds_d02_in['LTG3'][k:j,:,:], axis = 0)
    ds_d01['LTG3'][i, :, :] = np.nansum(ds_d01_in['LTG3_d01'][k:j,:,:], axis = 0)
    ds_d02['PR92W'][i, :, :] = np.nansum(ds_d02_in['PR92W'][k:j,:,:], axis = 0)
    ds_d01['PR92W'][i, :, :] = np.nansum(ds_d01_in['PR92W_d01'][k:j,:,:], axis = 0)
    ds_d01['CAPExP_R'][i, :, :] = np.nansum(ds_d01_in['CAPExP_R_d01'][k:j, :, :], axis=0)
    ds_d01['CAPExP_CSI'][i, :, :] = np.nansum(ds_d01_in['CAPExP_CSI_d01'][k:j, :, :], axis=0)
    ds_d02['CAPExP_R'][i, :, :] = np.nansum(ds_d02_in['CAPExP_R'][k:j, :, :], axis=0)
    ds_d02['CAPExP_CSI'][i, :, :] = np.nansum(ds_d02_in['CAPExP_CSI'][k:j, :, :], axis=0)
    ds_d01['Obs'][i, :, :] = np.nansum(L_d01[k:j, :, :], axis=0)
    ds_d02['Obs'][i, :, :] = np.nansum(L_d02[k:j, :, :], axis=0)

ds_d01.close()
ds_d02.close()
