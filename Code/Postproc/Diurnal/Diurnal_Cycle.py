# ---------------------------------------------------------------------------------------------
# MODULES
# ---------------------------------------------------------------------------------------------
import numpy as np
import pandas as pd
from netCDF4 import Dataset

# ---------------------------------------------------------------------------------------------
# CREATE NC FILE
# ---------------------------------------------------------------------------------------------
print('creating nc file')
# D01
ds_d01 = Dataset('/scratch/leuven/336/vsc33651/nu-wrf-dev/domain1_diurnal.nc', mode='w', format='NETCDF4')
ds_d01.createDimension('time', 24)
ds_d01.createVariable('time','int', dimensions=('time',),zlib=True)
ds_d01.variables['time'][:] = np.linspace(0,24, num=24, endpoint=False)
ds_d01.createVariable('LPI', 'f4', dimensions=('time',), zlib=True)
ds_d01.createVariable('LTG3', 'f4', dimensions=('time',), zlib=True)
ds_d01.createVariable('PR92H', 'f4', dimensions=('time',), zlib=True)
ds_d01.createVariable('PR92W', 'f4', dimensions=('time',), zlib=True)
ds_d01.createVariable('CAPExP', 'f4', dimensions=('time',), zlib=True)
ds_d01.createVariable('Obs', 'f4', dimensions=('time',), zlib=True)
# D02
ds_d02 = Dataset('/scratch/leuven/336/vsc33651/nu-wrf-dev/domain2_diurnal.nc', mode='w', format='NETCDF4')
ds_d02.createDimension('time', 24)
ds_d02.createVariable('time','int', dimensions=('time',),zlib=True)
ds_d02.variables['time'][:] = np.linspace(0,24, num=24, endpoint=False)
ds_d02.createVariable('LPI', 'f4', dimensions=('time',), zlib=True)
ds_d02.createVariable('LTG3', 'f4', dimensions=('time',), zlib=True)
ds_d02.createVariable('PR92H', 'f4', dimensions=('time',), zlib=True)
ds_d02.createVariable('PR92W', 'f4', dimensions=('time',), zlib=True)
ds_d02.createVariable('Obs', 'f4', dimensions=('time',), zlib=True)

# ---------------------------------------------------------------------------------------------
# LOAD DATA
# ---------------------------------------------------------------------------------------------
ds_d02_in = Dataset('/scratch/leuven/336/vsc33651/nu-wrf-dev/domain2_at_domain1_all_v4.nc', mode='r')
ds_d01_in = Dataset('/scratch/leuven/336/vsc33651/nu-wrf-dev/domain1_all.nc', mode='r')
ds_obs_in = Dataset('/scratch/leuven/336/vsc33651/nu-wrf-dev/CLDN_at_domain1_all_v4.nc','r')
L = ds_obs_in['Flashdensity_CC'][:].data + ds_obs_in['Flashdensity_CG'][:].data

time = pd.to_datetime('2000010100',format='%Y%m%d%H') + pd.to_timedelta(ds_d02_in['time'][:], unit='h')

# ---------------------------------------------------------------------------------------------
# DIURNAL CYCLE
# ---------------------------------------------------------------------------------------------
print('diurnal cycle')
# GET THE TIMES RIGHT
times = pd.to_datetime(time) - pd.to_timedelta(7, unit='h')
months = times.month
summer = np.where([(months > 5) & (months < 9)])
time_summer = times[summer[1]]

# EXTRACT THE RIGHT DATA
print('lpi_d02')
lpi_d02 = np.nanmean(ds_d02_in['LPI'][:], axis = (1,2))
lpi_summer_d02 = lpi_d02[summer[1]]
print('lpi_d01')
lpi_d01 = np.nanmean(ds_d01_in['LPI'][:], axis = (1,2))
lpi_summer_d01 = lpi_d01[summer[1]]
print('ltg_d01')
LTG3_d01 = np.nanmean(ds_d01_in['LTG3'][:], axis = (1,2))
LTG3_summer_d01 = LTG3_d01[summer[1]]
print('ltg_d02')
LTG3_d02 = np.nanmean(ds_d02_in['LTG3'][:], axis = (1,2))
LTG3_summer_d02 = LTG3_d02[summer[1]]
print('PR92H_d01')
PR92_H_d01 = np.nanmean(ds_d01_in['PR92H'][:], axis = (1,2))
PR92_H_summer_d01 = PR92_H_d01[summer[1]]
print('PR92H_d02')
PR92_H_d02 = np.nanmean(ds_d02_in['PR92H'][:], axis = (1,2))
PR92_H_summer_d02 = PR92_H_d02[summer[1]]
print('PR92W_d01')
PR92_W_d01 = np.nanmean(ds_d01_in['PR92W'][:], axis = (1,2))
PR92_W_summer_d01 = PR92_W_d01[summer[1]]
print('PR92W_d02')
PR92_W_d02 = np.nanmean(ds_d02_in['PR92W'][:], axis = (1,2))
PR92_W_summer_d02 = PR92_W_d02[summer[1]]
print('CAPExP_d01')
CAPExP_d01 = np.nanmean(ds_d01_in['CAPExP'][:], axis = (1,2))
CAPExP_summer_d01 = CAPExP_d01[summer[1]]
print('Obs_d01')
Obs_d01 = np.nanmean(L[:], axis = (1,2))
Obs_summer_d01 = Obs_d01[summer[1]]
print('Obs_d02')
Obs_d02 = np.nanmean(L[:], axis = (1,2))
Obs_summer_d02 = Obs_d02[summer[1]]

for i in range(0,24):
    print('loop'+str(i))
    time_select_vector = np.where(time_summer.hour == i)
    ds_d02['LPI'][i] = np.nanmean(lpi_summer_d02[time_select_vector[0]])
    ds_d01['LPI'][i] = np.nanmean(lpi_summer_d01[time_select_vector[0]])
    ds_d02['LTG3'][i] = np.nanmean(LTG3_summer_d02[time_select_vector[0]])
    ds_d01['LTG3'][i] = np.nanmean(LTG3_summer_d01[time_select_vector[0]])
    ds_d02['PR92W'][i] = np.nanmean(PR92_W_summer_d02[time_select_vector[0]])
    ds_d01['PR92W'][i] = np.nanmean(PR92_W_summer_d01[time_select_vector[0]])
    ds_d02['PR92H'][i] = np.nanmean(PR92_H_summer_d02[time_select_vector[0]])
    ds_d01['PR92H'][i] = np.nanmean(PR92_H_summer_d01[time_select_vector[0]])
    ds_d01['CAPExP'][i] = np.nanmean(CAPExP_summer_d01[time_select_vector[0]])
    ds_d02['Obs'][i] = np.nanmean(Obs_summer_d02[time_select_vector[0]])
    ds_d01['Obs'][i] = np.nanmean(Obs_summer_d01[time_select_vector[0]])

ds_d01.close()
ds_d02.close()
