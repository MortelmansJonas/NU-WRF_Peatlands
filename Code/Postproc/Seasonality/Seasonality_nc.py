# ---------------------------------------------------------------------------------------------
# MODULES
# ---------------------------------------------------------------------------------------------
import numpy as np
import pandas as pd
from netCDF4 import Dataset

# ---------------------------------------------------------------------------------------------
# LOAD DATA
# ---------------------------------------------------------------------------------------------
ds_d01_in = Dataset('/scratch/leuven/336/vsc33651/nu-wrf-dev/domain1_weekly_regridded.nc', 'r')
ds_d02_in = Dataset('/scratch/leuven/336/vsc33651/nu-wrf-dev/domain2_weekly_regridded.nc', 'r')

# ---------------------------------------------------------------------------------------------
# CREATE NC FILE
# ---------------------------------------------------------------------------------------------
ds_d01 = Dataset('/scratch/leuven/336/vsc33651/nu-wrf-dev/domain1_seasonality.nc', mode='w', format='NETCDF4')
ds_d02 = Dataset('/scratch/leuven/336/vsc33651/nu-wrf-dev/domain2_seasonality.nc', mode='w', format='NETCDF4')

ds_d01.createDimension('time', 13)
ds_d01.createVariable('time', 'int', dimensions=('time',),zlib=True)
ds_d01.createVariable('LPI', 'f4', dimensions=('time',), zlib=True)
ds_d01.createVariable('LTG3', 'f4', dimensions=('time',), zlib=True)
ds_d01.createVariable('PR92H', 'f4', dimensions=('time',), zlib=True)
ds_d01.createVariable('PR92W', 'f4', dimensions=('time',), zlib=True)
ds_d01.createVariable('CAPExP', 'f4', dimensions=('time',), zlib=True)
ds_d01.createVariable('Obs', 'f4', dimensions=('time',), zlib=True)

ds_d02.createDimension('time', 13)
ds_d02.createVariable('time', 'int', dimensions=('time',),zlib=True)
ds_d02.createVariable('LPI', 'f4', dimensions=('time',), zlib=True)
ds_d02.createVariable('LTG3', 'f4', dimensions=('time',), zlib=True)
ds_d02.createVariable('PR92H', 'f4', dimensions=('time',), zlib=True)
ds_d02.createVariable('PR92W', 'f4', dimensions=('time',), zlib=True)
ds_d02.createVariable('Obs', 'f4', dimensions=('time',), zlib=True)

timeseries = pd.date_range('2015-05-15', '2020-09-15', freq='D')
year = timeseries.year
month = timeseries.month

i2015 = np.where((year == 2015) & (month > 5) & (month < 9))
i2016 = np.where((year == 2016) & (month > 5) & (month < 9))
i2017 = np.where((year == 2017) & (month > 5) & (month < 9))
i2018 = np.where((year == 2018) & (month > 5) & (month < 9))
i2019 = np.where((year == 2019) & (month > 5) & (month < 9))
i2020 = np.where((year == 2020) & (month > 5) & (month < 9))

# ---------------------------------------------------------------------------------------------
# SEASONALITY 2.0
# ---------------------------------------------------------------------------------------------

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
    ds_d01['Obs'][i] = np.nanmean(ds_d01_in['Obs'][indices])
    ds_d02['Obs'][i] = np.nanmean(ds_d02_in['Obs'][indices])

ds_d01.close()
ds_d02.close()
