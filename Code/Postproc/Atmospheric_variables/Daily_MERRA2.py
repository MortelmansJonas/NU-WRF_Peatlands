# ---------------------------------------------------------------------------------------------
# MODULES
# ---------------------------------------------------------------------------------------------
import numpy as np
from netCDF4 import Dataset

# ---------------------------------------------------------------------------------------------
# LOAD DATA
# ---------------------------------------------------------------------------------------------
ds_d01 = Dataset('/scratch/leuven/336/vsc33651/nu-wrf-dev/domain1_all_MERRA2.nc', mode='r')
ds_d02 = Dataset('/scratch/leuven/336/vsc33651/nu-wrf-dev/domain2_all_MERRA2.nc', mode='r')
ds_M2 = Dataset('/scratch/leuven/336/vsc33651/nu-wrf-dev/MERRA2_at_domain1_all.nc', mode='r')
lat_d01 = ds_d01['lat'][:]
lon_d01 = ds_d01['lon'][:]
lat_d02 = ds_d02['lat'][:]
lon_d02 = ds_d02['lon'][:]
d01_precipitation = ds_d01['RAIN'][:]
d01_precipitation[d01_precipitation<0]=0
d02_precipitation = ds_d02['RAIN'][:]
d02_precipitation[d02_precipitation<0]=0
M2_precipitation = ds_M2['PRECTOT'][:] * 3600
# ---------------------------------------------------------------------------------------------
# CALCULATE DAILY INSTEAD OF HOURLY
# ---------------------------------------------------------------------------------------------
print('for loop')
days = np.divide(13248,24).astype(int)
M2_prec = np.zeros((days,63,109))
M2_LH = np.zeros((days,63,109))
M2_SH = np.zeros((days,63,109))
d01_prec = np.zeros((days,63,109))
d01_LH = np.zeros((days,63,109))
d01_SH = np.zeros((days,63,109))
d02_prec = np.zeros((days,170,309))
d02_LH = np.zeros((days,170,309))
d02_SH = np.zeros((days,170,309))

for i in range(0,days):
    print('loop' +str(i))
    j = (i+1)*24
    k = i*24
    M2_prec[i,:,:] = np.sum(M2_precipitation[k:j,:,:], axis = 0)
    M2_LH[i,:,:] = np.mean(ds_M2['LHLAND'][k:j,:,:], axis=0)
    M2_SH[i,:,:] = np.mean(ds_M2['SHLAND'][k:j,:,:], axis=0)
    d01_prec[i,:,:] = np.sum(d01_precipitation[k:j,:,:], axis = 0)
    d01_LH[i,:,:] = np.mean(ds_d01['LH'][k:j,:,:], axis=0)
    d01_SH[i,:,:] = np.mean(ds_d01['HFX'][k:j,:,:], axis=0)
    d02_prec[i,:,:] = np.sum(d02_precipitation[k:j,:,:], axis = 0)
    d02_LH[i,:,:] = np.mean(ds_d02['LH'][k:j,:,:], axis=0)
    d02_SH[i,:,:] = np.mean(ds_d02['HFX'][k:j,:,:], axis=0)

avg_M2_prec = np.nanmean(M2_prec, axis=0)
avg_M2_LH = np.nanmean(M2_LH, axis=0)
avg_M2_SH = np.nanmean(M2_SH, axis=0)
avg_d01_prec = np.nanmean(d01_prec, axis=0)
avg_d01_LH = np.nanmean(d01_LH, axis=0)
avg_d01_SH = np.nanmean(d01_SH, axis=0)
avg_d02_prec = np.nanmean(d02_prec, axis=0)
avg_d02_LH = np.nanmean(d02_LH, axis=0)
avg_d02_SH = np.nanmean(d02_SH, axis=0)

# ---------------------------------------------------------------------------------------------
# CREATE NC FILE
# ---------------------------------------------------------------------------------------------
print('creating nc file')
# D01
ds_d01 = Dataset('/scratch/leuven/336/vsc33651/nu-wrf-dev/domain1_M2.nc', mode='w', format='NETCDF4')
ds_d01.createDimension('time', 1)
ds_d01.createDimension('lat', 63)
ds_d01.createDimension('lon', 109)
ds_d01.createVariable('time','int', dimensions=('time',),zlib=True)
ds_d01.createVariable('lat',lat_d01.dtype, dimensions=('lat','lon',),zlib=True)
ds_d01.variables['lat'][:] = lat_d01
ds_d01.createVariable('lon',lon_d01.dtype, dimensions=('lat','lon',),zlib=True)
ds_d01.variables['lon'][:] = lon_d01
ds_d01.createVariable('M2_PREC', 'f4', dimensions=('time','lat','lon',), zlib=True)
ds_d01.createVariable('M2_LH', 'f4', dimensions=('time','lat','lon',), zlib=True)
ds_d01.createVariable('M2_SH', 'f4', dimensions=('time','lat','lon',), zlib=True)
ds_d01.createVariable('d01_PREC', 'f4', dimensions=('time','lat','lon',), zlib=True)
ds_d01.createVariable('d01_LH', 'f4', dimensions=('time','lat','lon',), zlib=True)
ds_d01.createVariable('d01_SH', 'f4', dimensions=('time','lat','lon',), zlib=True)
# D02
ds_d02 = Dataset('/scratch/leuven/336/vsc33651/nu-wrf-dev/domain2_M2.nc', mode='w', format='NETCDF4')
ds_d02.createDimension('time', 1)
ds_d02.createDimension('lat', 170)
ds_d02.createDimension('lon', 309)
ds_d02.createVariable('time','int', dimensions=('time',),zlib=True)
ds_d02.createVariable('lat',lat_d02.dtype, dimensions=('lat','lon',),zlib=True)
ds_d02.variables['lat'][:] = lat_d02
ds_d02.createVariable('lon',lon_d02.dtype, dimensions=('lat','lon',),zlib=True)
ds_d02.variables['lon'][:] = lon_d02
ds_d02.createVariable('d02_PREC', 'f4', dimensions=('time','lat','lon',), zlib=True)
ds_d02.createVariable('d02_LH', 'f4', dimensions=('time','lat','lon',), zlib=True)
ds_d02.createVariable('d02_SH', 'f4', dimensions=('time','lat','lon',), zlib=True)

ds_d01['M2_PREC'][:] = avg_M2_prec
ds_d01['M2_LH'][:] = avg_M2_LH
ds_d01['M2_SH'][:] = avg_M2_SH
ds_d01['d01_PREC'][:] = avg_d01_prec
ds_d01['d01_LH'][:] = avg_d01_LH
ds_d01['d01_SH'][:] = avg_d01_SH
ds_d02['d02_PREC'][:] = avg_d02_prec
ds_d02['d02_LH'][:] = avg_d02_LH
ds_d02['d02_SH'][:] = avg_d02_SH

ds_d01.close()
ds_d02.close()
