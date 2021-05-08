# ---------------------------------------------------------------------------------------------
# MODULES
# ---------------------------------------------------------------------------------------------
import numpy as np
from netCDF4 import Dataset

# ---------------------------------------------------------------------------------------------
# CREATE NC FILE
# ---------------------------------------------------------------------------------------------
ds_mccaul_2015_d01 = Dataset('/scratch/leuven/336/vsc33651/nu-wrf-dev/wrfout_nc_files/wrfout_d02_2015_v03.nc','r')
inds_lat = np.unique(np.where((ds_mccaul_2015_d01['lat'][:] > 58) & (ds_mccaul_2015_d01['lat'][:] <62))[0])
inds_lon = np.unique(np.where((ds_mccaul_2015_d01['lon'][:] > -123) & (ds_mccaul_2015_d01['lon'][:] <-108))[1])
lat = ds_mccaul_2015_d01['lat'][inds_lat,inds_lon]
lon = ds_mccaul_2015_d01['lon'][inds_lat,inds_lon]
ds = Dataset('/scratch/leuven/336/vsc33651/nu-wrf-dev/domain2_all_MERRA2.nc', mode='w', format='NETCDF4')
ds.createDimension('time', 13248)
ds.createDimension('lat', 170)
ds.createDimension('lon', 309)
ds.createVariable('time','int', dimensions=('time',),zlib=True)
ds.createVariable('lat',lat.dtype, dimensions=('lat','lon',),zlib=True)
ds.variables['lat'][:] = lat
ds.createVariable('lon',lon.dtype, dimensions=('lat','lon',),zlib=True)
ds.variables['lon'][:] = lon
ds.createVariable('LH', 'f4', dimensions=('time','lat','lon',), zlib=True)
ds.createVariable('HFX', 'f4', dimensions=('time','lat','lon',), zlib=True)
ds.createVariable('RAIN', 'f4', dimensions=('time','lat','lon',), zlib=True)

# ---------------------------------------------------------------------------------------------
# LOAD DATA
# ---------------------------------------------------------------------------------------------
print('2015')
ds_2015 = Dataset('/scratch/leuven/336/vsc33651/nu-wrf-dev/wrfout_nc_files/wrfout_d02_2015_MERRA2.nc', 'r')
ds['LH'][0:2208,:,:] = ds_2015['LH'][0:-1,inds_lat, inds_lon]
ds['HFX'][0:2208,:,:] = ds_2015['HFX'][0:-1,inds_lat, inds_lon]
TOTAL_PREC_2015 = ds_2015['RAINC'][0:-1, inds_lat, inds_lon] + ds_2015['RAINNC'][0:-1, inds_lat, inds_lon] + ds_2015['RAINSH'][0:-1, inds_lat, inds_lon]
ds['RAIN'][0,:,:] = TOTAL_PREC_2015[0,:,:]
for i in range(1,2208):
    print('loop_' + str(i))
    ds['RAIN'][i,:,:] = np.subtract(TOTAL_PREC_2015[i,:,:],TOTAL_PREC_2015[i-1,:,:])

print('2016')
ds_2016 = Dataset('/scratch/leuven/336/vsc33651/nu-wrf-dev/wrfout_nc_files/wrfout_d02_2016_MERRA2.nc', 'r')
TOTAL_PREC_2016 = ds_2016['RAINC'][0:-1, inds_lat, inds_lon] + ds_2016['RAINNC'][0:-1, inds_lat, inds_lon] + ds_2016['RAINSH'][0:-1, inds_lat, inds_lon]
ds['RAIN'][2208,:,:] = TOTAL_PREC_2016[0,:,:]
for i in range(1,2208):
    ds['RAIN'][i+2208,:,:] = np.subtract(TOTAL_PREC_2016[i,:,:],TOTAL_PREC_2016[i-1,:,:])

ds['LH'][2208:4416,:,:] = ds_2016['LH'][0:-1,inds_lat, inds_lon]
ds['HFX'][2208:4416,:,:] = ds_2016['HFX'][0:-1,inds_lat, inds_lon]

print('2017')
ds_2017 = Dataset('/scratch/leuven/336/vsc33651/nu-wrf-dev/wrfout_nc_files/wrfout_d02_2017_MERRA2.nc', 'r')
TOTAL_PREC_2017 = ds_2017['RAINC'][0:-1, inds_lat, inds_lon] + ds_2017['RAINNC'][0:-1, inds_lat, inds_lon] + ds_2017['RAINSH'][0:-1, inds_lat, inds_lon]
ds['RAIN'][4416,:,:] = TOTAL_PREC_2017[0,:,:]
for i in range(1,2208):
    ds['RAIN'][i+4416,:,:] = np.subtract(TOTAL_PREC_2017[i,:,:],TOTAL_PREC_2017[i-1,:,:])

ds['LH'][4416:6624,:,:] = ds_2017['LH'][0:-1,inds_lat, inds_lon]
ds['HFX'][4416,:6624,:] = ds_2017['HFX'][0:-1,inds_lat, inds_lon]

print('2018')
ds_2018 = Dataset('/scratch/leuven/336/vsc33651/nu-wrf-dev/wrfout_nc_files/wrfout_d02_2018_MERRA2.nc', 'r')
TOTAL_PREC_2018 = ds_2018['RAINC'][0:-1, inds_lat, inds_lon] + ds_2018['RAINNC'][0:-1, inds_lat, inds_lon] + ds_2018['RAINSH'][0:-1, inds_lat, inds_lon]
ds['RAIN'][6624,:,:] = TOTAL_PREC_2018[0,:,:]
for i in range(1,2208):
    ds['RAIN'][i+6624,:,:] = np.subtract(TOTAL_PREC_2018[i,:,:],TOTAL_PREC_2018[i-1,:,:])

ds['LH'][6624:8832,:,:] = ds_2018['LH'][0:-1,inds_lat, inds_lon]
ds['HFX'][6624:8832,:,:] = ds_2018['HFX'][0:-1,inds_lat, inds_lon]

print('2019')
ds_2019 = Dataset('/scratch/leuven/336/vsc33651/nu-wrf-dev/wrfout_nc_files/wrfout_d02_2019_MERRA2.nc', 'r')
TOTAL_PREC_2019 = ds_2019['RAINC'][0:-1, inds_lat, inds_lon] + ds_2019['RAINNC'][0:-1, inds_lat, inds_lon] + ds_2019['RAINSH'][0:-1, inds_lat, inds_lon]
ds['RAIN'][8832,:,:] = TOTAL_PREC_2019[0,:,:]
for i in range(1,2208):
    ds['RAIN'][i+8832,:,:] = np.subtract(TOTAL_PREC_2019[i,:,:],TOTAL_PREC_2019[i-1,:,:])

ds['LH'][8832:11040,:,:] = ds_2019['LH'][0:-1,inds_lat, inds_lon]
ds['HFX'][8832:11040,:,:] = ds_2019['HFX'][0:-1,inds_lat, inds_lon]

print('2020')
ds_2020 = Dataset('/scratch/leuven/336/vsc33651/nu-wrf-dev/wrfout_nc_files/wrfout_d02_2020_MERRA2.nc', 'r')
TOTAL_PREC_2020 = ds_2020['RAINC'][0:-1, inds_lat, inds_lon] + ds_2020['RAINNC'][0:-1, inds_lat, inds_lon] + ds_2020['RAINSH'][0:-1, inds_lat, inds_lon]
ds['RAIN'][11040,:,:] = TOTAL_PREC_2020[0,:,:]
for i in range(1,2208):
    ds['RAIN'][i+11040,:,:] = np.subtract(TOTAL_PREC_2020[i,:,:],TOTAL_PREC_2020[i-1,:,:])

ds['LH'][11040:13248,:,:] = ds_2020['LH'][0:-1,inds_lat, inds_lon]
ds['HFX'][11040:13248,:,:] = ds_2020['HFX'][0:-1,inds_lat, inds_lon]

ds.close()
