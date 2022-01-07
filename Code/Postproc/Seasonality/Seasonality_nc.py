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
ds_d02_in = Dataset('/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/wrfout_nc_files/domain2_weekly_regridded.nc', 'r')
ds_d01_in = Dataset('/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/wrfout_nc_files/domain1_weekly_regridded.nc', 'r')

# ---------------------------------------------------------------------------------------------
# CREATE NC FILE
# ---------------------------------------------------------------------------------------------
ds_d01 = Dataset('/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/wrfout_nc_files/d01_seasonality.nc', 'w')
ds_d01.createDimension('time', 13)
ds_d01.createVariable('time', 'int', dimensions=('time',),zlib=True)
ds_d01.createVariable('LPI', 'f4', dimensions=('time',), zlib=True)
ds_d01.createVariable('LTG3', 'f4', dimensions=('time',), zlib=True)
ds_d01.createVariable('PR92W', 'f4', dimensions=('time',), zlib=True)
ds_d01.createVariable('CAPExP_R', 'f4', dimensions=('time',), zlib=True)
ds_d01.createVariable('CAPExP_CSI', 'f4', dimensions=('time',), zlib=True)
ds_d01.createVariable('Obs', 'f4', dimensions=('time',), zlib=True)
# D02
ds_d02 = Dataset('/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/wrfout_nc_files/d02_seasonality.nc', 'w')
ds_d02.createDimension('time', 13)
ds_d02.createVariable('time', 'int', dimensions=('time',),zlib=True)
ds_d02.createVariable('LPI', 'f4', dimensions=('time',), zlib=True)
ds_d02.createVariable('LTG3', 'f4', dimensions=('time',), zlib=True)
ds_d02.createVariable('PR92W', 'f4', dimensions=('time',), zlib=True)
ds_d02.createVariable('CAPExP_R', 'f4', dimensions=('time',), zlib=True)
ds_d02.createVariable('CAPExP_CSI', 'f4', dimensions=('time',), zlib=True)
ds_d02.createVariable('Obs', 'f4', dimensions=('time',), zlib=True)

timeseries= pd.date_range('2015-05-15', '2020-09-15', freq = 'D')
year = timeseries.year
month=timeseries.month

# Only get the lightning season
i2015 = np.where((year == 2015) & (month > 5) & (month < 9))
i2016 = np.where((year == 2016) & (month > 5) & (month < 9))
i2017 = np.where((year == 2017) & (month > 5) & (month < 9))
i2018 = np.where((year == 2018) & (month > 5) & (month < 9))
i2019 = np.where((year == 2019) & (month > 5) & (month < 9))
i2020 = np.where((year == 2020) & (month > 5) & (month < 9))

# ---------------------------------------------------------------------------------------------
# SEASONALITY
# ---------------------------------------------------------------------------------------------
for i in range(0,13):
    indices = np.arange(i,78,13)
    ds_d01['LPI'][i] = np.nanmean(ds_d01_in['LPI'][indices])
    ds_d02['LPI'][i] = np.nanmean(ds_d02_in['LPI'][indices])
    ds_d01['LTG3'][i] = np.nanmean(ds_d01_in['LTG3'][indices])
    ds_d02['LTG3'][i] = np.nanmean(ds_d02_in['LTG3'][indices])
    ds_d01['PR92W'][i] = np.nanmean(ds_d01_in['PR92W'][indices])
    ds_d02['PR92W'][i] = np.nanmean(ds_d02_in['PR92W'][indices])
    ds_d01['CAPExP_R'][i] = np.nanmean(ds_d01_in['CAPExP_R'][indices])
    ds_d02['CAPExP_R'][i] = np.nanmean(ds_d02_in['CAPExP_R'][indices])
    ds_d01['CAPExP_CSI'][i] = np.nanmean(ds_d01_in['CAPExP_CSI'][indices])
    ds_d02['CAPExP_CSI'][i] = np.nanmean(ds_d02_in['CAPExP_CSI'][indices])
    ds_d01['Obs'][i] = np.nanmean(ds_d01_in['Obs'][indices])
    ds_d02['Obs'][i] = np.nanmean(ds_d02_in['Obs'][indices])

ds_d01.close()
ds_d02.close()
