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
ds_d01 = Dataset('/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/NU-WRF/wrfout_nc_files/d01_diurnal_Thompson.nc', 'w')
ds_d01.createDimension('time', 24)
ds_d01.createVariable('time', 'int', dimensions=('time',),zlib=True)
ds_d01.variables['time'][:] = np.linspace(0,24,num=24, endpoint=False)
ds_d01.createVariable('LPI', 'f4', dimensions=('time',), zlib=True)
ds_d01.createVariable('LTG3', 'f4', dimensions=('time',), zlib=True)
ds_d01.createVariable('PR92W', 'f4', dimensions=('time',), zlib=True)
ds_d01.createVariable('CAPExP_R', 'f4', dimensions=('time',), zlib=True)
ds_d01.createVariable('Obs', 'f4', dimensions=('time',), zlib=True)
# D02
ds_d02 = Dataset('/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/NU-WRF/wrfout_nc_files/d02_diurnal_Thompson.nc', 'w')
ds_d02.createDimension('time', 24)
ds_d02.createVariable('time', 'int', dimensions=('time',),zlib=True)
ds_d02.variables['time'][:] = np.linspace(0,24,num=24, endpoint=False)
ds_d02.createVariable('LPI', 'f4', dimensions=('time',), zlib=True)
ds_d02.createVariable('LTG3', 'f4', dimensions=('time',), zlib=True)
ds_d02.createVariable('PR92W', 'f4', dimensions=('time',), zlib=True)
ds_d02.createVariable('CAPExP_R', 'f4', dimensions=('time',), zlib=True)
ds_d02.createVariable('Obs', 'f4', dimensions=('time',), zlib=True)

# ---------------------------------------------------------------------------------------------
# LOAD DATA
# ---------------------------------------------------------------------------------------------
print('load data')
ds_d01_in = Dataset('/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/NU-WRF/wrfout_nc_files/data_calibrated_ax_d01_Thompson.nc', 'r')

time = pd.to_datetime('2000010100', format='%Y%m%d%H') + pd.to_timedelta(ds_d01_in['time'][:], unit='h')

# ---------------------------------------------------------------------------------------------
# DIURNAL CYCLE
# ---------------------------------------------------------------------------------------------
print('diurnal cycle')
# Get the times right (convert to local time so that the peak is in the local afternoon)
times = pd.to_datetime(time) - pd.to_timedelta(7, unit='h')
months = times.month
summer = np.where([(months > 5) & (months <9)])
time_summer = times[summer[1]]

# Extract the right data
print('lpi_d01')
lpi_d01 = np.nanmean(ds_d01_in['LPI_d01'][:], axis=(1,2))
lpi_summer_d01 = lpi_d01[summer[1]]
print('lpi_d02')
lpi_d02 = np.nanmean(ds_d01_in['LPI_d02'][:], axis=(1,2))
lpi_summer_d02 = lpi_d02[summer[1]]
print('ltg3_d01')
ltg3_d01 = np.nanmean(ds_d01_in['LTG3_d01'][:], axis=(1,2))
ltg3_summer_d01 = ltg3_d01[summer[1]]
print('ltg3_d02')
ltg3_d02 = np.nanmean(ds_d01_in['LTG3_d02'][:], axis=(1,2))
ltg3_summer_d02 = ltg3_d02[summer[1]]
print('pr92w_d01')
pr92w_d01 = np.nanmean(ds_d01_in['PR92W_d01'][:], axis=(1,2))
pr92w_summer_d01 = pr92w_d01[summer[1]]
print('pr92w_d02')
pr92w_d02 = np.nanmean(ds_d01_in['PR92W_d02'][:], axis=(1,2))
pr92w_summer_d02 = pr92w_d02[summer[1]]
print('capexp_R_d01')
capexp_r_d01 = np.nanmean(ds_d01_in['CAPExP_R_d01'][:], axis=(1,2))
capexp_r_summer_d01 = capexp_r_d01[summer[1]]
print('capexp_R_d02')
capexp_r_d02 = np.nanmean(ds_d01_in['CAPExP_R_d02'][:], axis=(1,2))
capexp_r_summer_d02 = capexp_r_d02[summer[1]]
print('obs_d01')
obs_d01 = np.nanmean(ds_d01_in['Obs'][:], axis=(1,2))
obs_summer_d01 = obs_d01[summer[1]]
print('obs_d02')
obs_d02 = np.nanmean(ds_d01_in['Obs'][:], axis=(1,2))
obs_summer_d02 = obs_d02[summer[1]]

# ---------------------------------------------------------------------------------------------
# FILL .nc FILE
# ---------------------------------------------------------------------------------------------
for i in range(0,24):
    print('loop ' + str(i))
    time_select_vector = np.where(time_summer.hour == i)
    ds_d01['LPI'][i] = np.nanmean(lpi_summer_d01[time_select_vector[0]])
    ds_d02['LPI'][i] = np.nanmean(lpi_summer_d02[time_select_vector[0]])
    ds_d01['LTG3'][i] = np.nanmean(ltg3_summer_d01[time_select_vector[0]])
    ds_d02['LTG3'][i] = np.nanmean(ltg3_summer_d02[time_select_vector[0]])
    ds_d01['PR92W'][i] = np.nanmean(pr92w_summer_d01[time_select_vector[0]])
    ds_d02['PR92W'][i] = np.nanmean(pr92w_summer_d02[time_select_vector[0]])
    ds_d01['CAPExP_R'][i] = np.nanmean(capexp_r_summer_d01[time_select_vector[0]])
    ds_d02['CAPExP_R'][i] = np.nanmean(capexp_r_summer_d02[time_select_vector[0]])
    ds_d01['Obs'][i] = np.nanmean(obs_summer_d01[time_select_vector[0]])
    ds_d02['Obs'][i] = np.nanmean(obs_summer_d02[time_select_vector[0]])

ds_d01.close()
ds_d02.close()
