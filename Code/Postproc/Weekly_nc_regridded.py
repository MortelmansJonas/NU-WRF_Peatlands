# ---------------------------------------------------------------------------------------------
# MODULES
# ---------------------------------------------------------------------------------------------
import numpy as np
import pandas as pd
from netCDF4 import Dataset

# ---------------------------------------------------------------------------------------------
# LOAD DATA OF D02
# ---------------------------------------------------------------------------------------------
ds_d01_daily = Dataset('/scratch/leuven/336/vsc33651/nu-wrf-dev/domain1_daily_regridded.nc', mode='r')
ds_d02_daily = Dataset('/scratch/leuven/336/vsc33651/nu-wrf-dev/domain2_daily_regridded.nc', mode='r')
# ---------------------------------------------------------------------------------------------
# CREATE NC FILE
# ---------------------------------------------------------------------------------------------
print('creating nc file')
# D01
timeseries = pd.date_range('2015-05-15','2020-09-15', freq='D')
ds_d01 = Dataset('/scratch/leuven/336/vsc33651/nu-wrf-dev/domain1_weekly_regridded.nc', mode='w', format='NETCDF4')
ds_d01.createDimension('weeks', 78)
ds_d01.createDimension('timeseries', 1951)
ds_d01.createVariable('weeks','int', dimensions=('weeks',),zlib=True)
ds_d01.createVariable('timeseries','int', dimensions=('timeseries',),zlib=True)
ds_d01.createVariable('LPI', 'f4', dimensions=('weeks',), zlib=True)
ds_d01.createVariable('LTG3', 'f4', dimensions=('weeks',), zlib=True)
ds_d01.createVariable('PR92H', 'f4', dimensions=('weeks',), zlib=True)
ds_d01.createVariable('PR92W', 'f4', dimensions=('weeks',), zlib=True)
ds_d01.createVariable('CAPExP', 'f4', dimensions=('weeks',), zlib=True)
ds_d01.createVariable('Obs', 'f4', dimensions=('weeks',), zlib=True)
ds_d01.createVariable('data_LPI', 'f4', dimensions=('timeseries',), zlib=True)
ds_d01.createVariable('data_LTG3', 'f4', dimensions=('timeseries',), zlib=True)
ds_d01.createVariable('data_PR92H', 'f4', dimensions=('timeseries',), zlib=True)
ds_d01.createVariable('data_PR92W', 'f4', dimensions=('timeseries',), zlib=True)
ds_d01.createVariable('data_CAPExP', 'f4', dimensions=('timeseries',), zlib=True)
ds_d01.createVariable('data_Obs', 'f4', dimensions=('timeseries',), zlib=True)
# D02
ds_d02 = Dataset('/scratch/leuven/336/vsc33651/nu-wrf-dev/domain2_weekly_regridded.nc', mode='w', format='NETCDF4')
ds_d02.createDimension('weeks', 78)
ds_d02.createDimension('timeseries', 1951 )
ds_d02.createVariable('weeks','int', dimensions=('weeks',),zlib=True)
ds_d02.createVariable('timeseries', 'int', dimensions=('timeseries',),zlib=True)
ds_d02.createVariable('LPI', 'f4', dimensions=('weeks',), zlib=True)
ds_d02.createVariable('LTG3', 'f4', dimensions=('weeks',), zlib=True)
ds_d02.createVariable('PR92H', 'f4', dimensions=('weeks',), zlib=True)
ds_d02.createVariable('PR92W', 'f4', dimensions=('weeks',), zlib=True)
ds_d02.createVariable('Obs', 'f4', dimensions=('weeks',), zlib=True)
ds_d02.createVariable('data_LPI', 'f4', dimensions=('timeseries',), zlib=True)
ds_d02.createVariable('data_LTG3', 'f4', dimensions=('timeseries',), zlib=True)
ds_d02.createVariable('data_PR92H', 'f4', dimensions=('timeseries',), zlib=True)
ds_d02.createVariable('data_PR92W', 'f4', dimensions=('timeseries',), zlib=True)
ds_d02.createVariable('data_Obs', 'f4', dimensions=('timeseries',), zlib=True)
# ---------------------------------------------------------------------------------------------
# FILL NC FILE
# ---------------------------------------------------------------------------------------------
timeseries_unit = 'days since 2000-01-01 00:00'
weeks_unit = 'weeks since 2000-01-01 00:00'
ds_d02['weeks'][0:13] = pd.date_range('2015-06-01 00:00:00','2015-08-31 00:00:00', freq='W')
ds_d02['weeks'][13:26] = pd.date_range('2016-06-01 00:00:00','2016-08-31 00:00:00', freq='W')
ds_d02['weeks'][26:39] = pd.date_range('2017-06-01 00:00:00','2017-08-31 00:00:00', freq='W')
ds_d02['weeks'][39:52] = pd.date_range('2018-06-01','2018-08-31 00:00:00', freq='W')
ds_d02['weeks'][52:65] = pd.date_range('2019-06-01 00:00:00','2019-08-31 00:00:00', freq='W')
ds_d02['weeks'][65:78] = pd.date_range('2020-06-01 00:00:00','2020-08-31 00:00:00', freq='W')
ds_d01['weeks'][0:13] = pd.date_range('2015-06-01 00:00:00','2015-08-31 00:00:00', freq='W')
ds_d01['weeks'][13:26] = pd.date_range('2016-06-01 00:00:00','2016-08-31 00:00:00', freq='W')
ds_d01['weeks'][26:39] = pd.date_range('2017-06-01 00:00:00','2017-08-31 00:00:00', freq='W')
ds_d01['weeks'][39:52] = pd.date_range('2018-06-01 00:00:00','2018-08-31 00:00:00', freq='W')
ds_d01['weeks'][52:65] = pd.date_range('2019-06-01 00:00:00','2019-08-31 00:00:00', freq='W')
ds_d01['weeks'][65:782] = pd.date_range('2020-06-01 00:00:00','2020-08-31 00:00:00', freq='W')
ds_d01['timeseries'][:] = timeseries
ds_d02['timeseries'][:] = timeseries
ds_d02.variables['timeseries'].setncatts({'long_name': 'timeseries', 'units': timeseries_unit})
ds_d01.variables['timeseries'].setncatts({'long_name': 'timeseries', 'units': timeseries_unit})
ds_d02.variables['weeks'].setncatts({'long_name': 'weeks', 'units':weeks_unit})
ds_d01.variables['weeks'].setncatts({'long_name': 'weeks', 'units': weeks_unit})

# SET EVERYTHING TO 0
ds_d02['LPI'][:] = 0
ds_d02['LTG3'][:] = 0
ds_d02['PR92H'][:] = 0
ds_d02['PR92W'][:] = 0
ds_d02['Obs'][:] = 0
ds_d01['LPI'][:] = 0
ds_d01['LTG3'][:] = 0
ds_d01['PR92H'][:] = 0
ds_d01['PR92W'][:] = 0
ds_d01['CAPExP'][:] = 0
ds_d01['Obs'][:] = 0
ds_d02['data_LPI'][:] = 0
ds_d02['data_LTG3'][:] = 0
ds_d02['data_PR92H'][:] = 0
ds_d02['data_PR92W'][:] = 0
ds_d02['data_Obs'][:] = 0
ds_d01['data_LPI'][:] = 0
ds_d01['data_LTG3'][:] = 0
ds_d01['data_PR92H'][:] = 0
ds_d01['data_PR92W'][:] = 0
ds_d01['data_CAPExP'][:] = 0
ds_d01['data_Obs'][:] = 0


print('time')
year = timeseries.year
month = timeseries.month
i2015 = np.where((year == 2015) & (month > 5) & (month < 9))
i2016 = np.where((year == 2016) & (month > 5) & (month < 9))
i2017 = np.where((year == 2017) & (month > 5) & (month < 9))
i2018 = np.where((year == 2018) & (month > 5) & (month < 9))
i2019 = np.where((year == 2019) & (month > 5) & (month < 9))
i2020 = np.where((year == 2020) & (month > 5) & (month < 9))
# LPI DOMAIN 1
print('LPI DATA')
ds_d01['data_LPI'][i2015] = np.nanmean(ds_d01_daily['LPI'][0:92], axis=(1,2))
ds_d01['data_LPI'][i2016] = np.nanmean(ds_d01_daily['LPI'][92:184], axis=(1,2))
ds_d01['data_LPI'][i2017] = np.nanmean(ds_d01_daily['LPI'][184:276], axis=(1,2))
ds_d01['data_LPI'][i2018] = np.nanmean(ds_d01_daily['LPI'][276:368], axis=(1,2))
ds_d01['data_LPI'][i2019] = np.nanmean(ds_d01_daily['LPI'][368:460], axis=(1,2))
ds_d01['data_LPI'][i2020] = np.nanmean(ds_d01_daily['LPI'][460:552], axis=(1,2))
print(ds_d01['data_LPI'][:])
print(np.amax(ds_d01['data_LPI'][:]))
# LPI DOMAIN 2
ds_d02['data_LPI'][i2015] = np.mean(ds_d02_daily['LPI'][0:92], axis=(1,2))
ds_d02['data_LPI'][i2016] = np.mean(ds_d02_daily['LPI'][92:184], axis=(1,2))
ds_d02['data_LPI'][i2017] = np.mean(ds_d02_daily['LPI'][184:276], axis=(1,2))
ds_d02['data_LPI'][i2018] = np.mean(ds_d02_daily['LPI'][276:368], axis=(1,2))
ds_d02['data_LPI'][i2019] = np.mean(ds_d02_daily['LPI'][368:460], axis=(1,2))
ds_d02['data_LPI'][i2020] = np.mean(ds_d02_daily['LPI'][460:552], axis=(1,2))
# LTG3 DOMAIN 1
print('LTG3 DATA')
ds_d01['data_LTG3'][i2015] = np.mean(ds_d01_daily['LTG3'][0:92], axis=(1,2))
ds_d01['data_LTG3'][i2016] = np.mean(ds_d01_daily['LTG3'][92:184], axis=(1,2))
ds_d01['data_LTG3'][i2017] = np.mean(ds_d01_daily['LTG3'][184:276], axis=(1,2))
ds_d01['data_LTG3'][i2018] = np.mean(ds_d01_daily['LTG3'][276:368], axis=(1,2))
ds_d01['data_LTG3'][i2019] = np.mean(ds_d01_daily['LTG3'][368:460], axis=(1,2))
ds_d01['data_LTG3'][i2020] = np.mean(ds_d01_daily['LTG3'][460:552], axis=(1,2))
# LTG3 DOMAIN 2
ds_d02['data_LTG3'][i2015] = np.mean(ds_d02_daily['LTG3'][0:92], axis=(1,2))
ds_d02['data_LTG3'][i2016] = np.mean(ds_d02_daily['LTG3'][92:184], axis=(1,2))
ds_d02['data_LTG3'][i2017] = np.mean(ds_d02_daily['LTG3'][184:276], axis=(1,2))
ds_d02['data_LTG3'][i2018] = np.mean(ds_d02_daily['LTG3'][276:368], axis=(1,2))
ds_d02['data_LTG3'][i2019] = np.mean(ds_d02_daily['LTG3'][368:460], axis=(1,2))
ds_d02['data_LTG3'][i2020] = np.mean(ds_d02_daily['LTG3'][460:552], axis=(1,2))
# PR92H DOMAIN 1
print('PR92H DATA')
ds_d01['data_PR92H'][i2015] = np.nanmean(ds_d01_daily['PR92H'][0:92], axis=(1,2))
ds_d01['data_PR92H'][i2016] = np.nanmean(ds_d01_daily['PR92H'][92:184], axis=(1,2))
ds_d01['data_PR92H'][i2017] = np.nanmean(ds_d01_daily['PR92H'][184:276], axis=(1,2))
ds_d01['data_PR92H'][i2018] = np.nanmean(ds_d01_daily['PR92H'][276:368], axis=(1,2))
ds_d01['data_PR92H'][i2019] = np.nanmean(ds_d01_daily['PR92H'][368:460], axis=(1,2))
ds_d01['data_PR92H'][i2020] = np.nanmean(ds_d01_daily['PR92H'][460:552], axis=(1,2))
# PR92H DOMAIN 2
ds_d02['data_PR92H'][i2015] = np.nanmean(ds_d02_daily['PR92H'][0:92], axis=(1,2))
ds_d02['data_PR92H'][i2016] = np.nanmean(ds_d02_daily['PR92H'][92:184], axis=(1,2))
ds_d02['data_PR92H'][i2017] = np.nanmean(ds_d02_daily['PR92H'][184:276], axis=(1,2))
ds_d02['data_PR92H'][i2018] = np.nanmean(ds_d02_daily['PR92H'][276:368], axis=(1,2))
ds_d02['data_PR92H'][i2019] = np.nanmean(ds_d02_daily['PR92H'][368:460], axis=(1,2))
ds_d02['data_PR92H'][i2020] = np.nanmean(ds_d02_daily['PR92H'][460:552], axis=(1,2))
# PR92W DOMAIN 1
print('PR92W DATA')
ds_d01['data_PR92W'][i2015] = np.mean(ds_d01_daily['PR92W'][0:92], axis=(1,2))
ds_d01['data_PR92W'][i2016] = np.mean(ds_d01_daily['PR92W'][92:184], axis=(1,2))
ds_d01['data_PR92W'][i2017] = np.mean(ds_d01_daily['PR92W'][184:276], axis=(1,2))
ds_d01['data_PR92W'][i2018] = np.mean(ds_d01_daily['PR92W'][276:368], axis=(1,2))
ds_d01['data_PR92W'][i2019] = np.mean(ds_d01_daily['PR92W'][368:460], axis=(1,2))
ds_d01['data_PR92W'][i2020] = np.mean(ds_d01_daily['PR92W'][460:552], axis=(1,2))
# PR92W DOMAIN 2
ds_d02['data_PR92W'][i2015] = np.mean(ds_d02_daily['PR92W'][0:92], axis=(1,2))
ds_d02['data_PR92W'][i2016] = np.mean(ds_d02_daily['PR92W'][92:184], axis=(1,2))
ds_d02['data_PR92W'][i2017] = np.mean(ds_d02_daily['PR92W'][184:276], axis=(1,2))
ds_d02['data_PR92W'][i2018] = np.mean(ds_d02_daily['PR92W'][276:368], axis=(1,2))
ds_d02['data_PR92W'][i2019] = np.mean(ds_d02_daily['PR92W'][368:460], axis=(1,2))
ds_d02['data_PR92W'][i2020] = np.mean(ds_d02_daily['PR92W'][460:552], axis=(1,2))
# CAPExP DOMAIN 1
print('CAPExP DATA')
ds_d01['data_CAPExP'][i2015] = np.mean(ds_d01_daily['CAPExP'][0:92], axis=(1,2))
ds_d01['data_CAPExP'][i2016] = np.mean(ds_d01_daily['CAPExP'][92:184], axis=(1,2))
ds_d01['data_CAPExP'][i2017] = np.mean(ds_d01_daily['CAPExP'][184:276], axis=(1,2))
ds_d01['data_CAPExP'][i2018] = np.mean(ds_d01_daily['CAPExP'][276:368], axis=(1,2))
ds_d01['data_CAPExP'][i2019] = np.mean(ds_d01_daily['CAPExP'][368:460], axis=(1,2))
ds_d01['data_CAPExP'][i2020] = np.mean(ds_d01_daily['CAPExP'][460:552], axis=(1,2))
# Obs DOMAIN 1
print('Obs DATA')
ds_d01['data_Obs'][i2015] = np.nanmean(ds_d01_daily['Obs'][0:92], axis=(1,2))
ds_d01['data_Obs'][i2016] = np.nanmean(ds_d01_daily['Obs'][92:184], axis=(1,2))
ds_d01['data_Obs'][i2017] = np.nanmean(ds_d01_daily['Obs'][184:276], axis=(1,2))
ds_d01['data_Obs'][i2018] = np.nanmean(ds_d01_daily['Obs'][276:368], axis=(1,2))
ds_d01['data_Obs'][i2019] = np.nanmean(ds_d01_daily['Obs'][368:460], axis=(1,2))
ds_d01['data_Obs'][i2020] = np.nanmean(ds_d01_daily['Obs'][460:552], axis=(1,2))
print(ds_d01['data_Obs'][:])
print(np.amax(ds_d01['data_Obs'][:]))
# Obs DOMAIN 2
print('Obs DATA')
ds_d02['data_Obs'][i2015] = np.nanmean(ds_d02_daily['Obs'][0:92], axis=(1,2))
ds_d02['data_Obs'][i2016] = np.nanmean(ds_d02_daily['Obs'][92:184], axis=(1,2))
ds_d02['data_Obs'][i2017] = np.nanmean(ds_d02_daily['Obs'][184:276], axis=(1,2))
ds_d02['data_Obs'][i2018] = np.nanmean(ds_d02_daily['Obs'][276:368], axis=(1,2))
ds_d02['data_Obs'][i2019] = np.nanmean(ds_d02_daily['Obs'][368:460], axis=(1,2))
ds_d02['data_Obs'][i2020] = np.nanmean(ds_d02_daily['Obs'][460:552], axis=(1,2))
print(ds_d02['data_Obs'][:])
print(np.amax(ds_d02['data_Obs'][:]))


print('np.where DATA')
ds_d01['data_LPI'][:] = np.where(month == 6, ds_d01['data_LPI'][:],
                np.where(month == 7, ds_d01['data_LPI'][:],
                         np.where(month == 8, ds_d01['data_LPI'][:], np.nan)))
ds_d02['data_LPI'][:] = np.where(month == 6, ds_d02['data_LPI'][:],
                np.where(month == 7, ds_d02['data_LPI'][:],
                         np.where(month == 8, ds_d02['data_LPI'][:], np.nan)))
ds_d01['data_LTG3'][:] = np.where(month == 6, ds_d01['data_LTG3'][:],
                np.where(month == 7, ds_d01['data_LTG3'][:],
                         np.where(month == 8, ds_d01['data_LTG3'][:], np.nan)))
ds_d02['data_LTG3'][:] = np.where(month == 6, ds_d02['data_LTG3'][:],
                np.where(month == 7, ds_d02['data_LTG3'][:],
                         np.where(month == 8, ds_d02['data_LTG3'][:], np.nan)))
ds_d01['data_PR92H'][:] = np.where(month == 6, ds_d01['data_PR92H'][:],
                np.where(month == 7, ds_d01['data_PR92H'][:],
                         np.where(month == 8, ds_d01['data_PR92H'][:], np.nan)))
ds_d02['data_PR92H'][:] = np.where(month == 6, ds_d02['data_PR92H'][:],
                np.where(month == 7, ds_d02['data_PR92H'][:],
                         np.where(month == 8, ds_d02['data_PR92H'][:], np.nan)))
ds_d01['data_PR92W'][:] = np.where(month == 6, ds_d01['data_PR92W'][:],
                np.where(month == 7, ds_d01['data_PR92W'][:],
                         np.where(month == 8, ds_d01['data_PR92W'][:], np.nan)))
ds_d02['data_PR92W'][:] = np.where(month == 6, ds_d02['data_PR92W'][:],
                np.where(month == 7, ds_d02['data_PR92W'][:],
                         np.where(month == 8, ds_d02['data_PR92W'][:], np.nan)))
ds_d01['data_CAPExP'][:] = np.where(month == 6, ds_d01['data_CAPExP'][:],
                np.where(month == 7, ds_d01['data_CAPExP'][:],
                         np.where(month == 8, ds_d01['data_CAPExP'][:], np.nan)))
ds_d01['data_Obs'][:] = np.where(month == 6, ds_d01['data_Obs'][:],
                np.where(month == 7, ds_d01['data_Obs'][:],
                         np.where(month == 8, ds_d01['data_Obs'][:], np.nan)))
ds_d02['data_Obs'][:] = np.where(month == 6, ds_d02['data_Obs'][:],
                np.where(month == 7, ds_d02['data_Obs'][:],
                         np.where(month == 8, ds_d02['data_Obs'][:], np.nan)))

# WEEKLY AVERAGE
print('weekly average')
weekly_LPI_d02 = np.zeros((546))
weekly_LPI_d01 = np.zeros((546))
weekly_LTG3_d02 = np.zeros((546))
weekly_LTG3_d01 = np.zeros((546))
weekly_PR92_H_d02 = np.zeros((546))
weekly_PR92_H_d01 = np.zeros((546))
weekly_PR92_W_d02 = np.zeros((546))
weekly_PR92_W_d01 = np.zeros((546))
weekly_CAPExP_d01 = np.zeros((546))
weekly_Obs_d02 = np.zeros((546))
weekly_Obs_d01 = np.zeros((546))
# LPI DOMAIN 2
print('weekly average LPI')
weekly_LPI_d02[0:91] = ds_d02['data_LPI'][17:108]
weekly_LPI_d02[91:182] = ds_d02['data_LPI'][383:474]
weekly_LPI_d02[182:273] = ds_d02['data_LPI'][748:839]
weekly_LPI_d02[273:364] = ds_d02['data_LPI'][1113:1204]
weekly_LPI_d02[364:455] = ds_d02['data_LPI'][1478:1569]
weekly_LPI_d02[455:546] = ds_d02['data_LPI'][1844:1935]
# LPI DOMAIN 1
weekly_LPI_d01[0:91] = ds_d01['data_LPI'][17:108]
weekly_LPI_d01[91:182] = ds_d01['data_LPI'][383:474]
weekly_LPI_d01[182:273] = ds_d01['data_LPI'][748:839]
weekly_LPI_d01[273:364] = ds_d01['data_LPI'][1113:1204]
weekly_LPI_d01[364:455] = ds_d01['data_LPI'][1478:1569]
weekly_LPI_d01[455:546] = ds_d01['data_LPI'][1844:1935]
# LTG3 DOMAIN 2
print('weekly average LTG')
weekly_LTG3_d02[0:91] = ds_d02['data_LTG3'][17:108]
weekly_LTG3_d02[91:182] = ds_d02['data_LTG3'][383:474]
weekly_LTG3_d02[182:273] = ds_d02['data_LTG3'][748:839]
weekly_LTG3_d02[273:364] = ds_d02['data_LTG3'][1113:1204]
weekly_LTG3_d02[364:455] = ds_d02['data_LTG3'][1478:1569]
weekly_LTG3_d02[455:546] = ds_d02['data_LTG3'][1844:1935]
# LTG3 DOMAIN 1
weekly_LTG3_d01[0:91] = ds_d01['data_LTG3'][17:108]
weekly_LTG3_d01[91:182] = ds_d01['data_LTG3'][383:474]
weekly_LTG3_d01[182:273] = ds_d01['data_LTG3'][748:839]
weekly_LTG3_d01[273:364] = ds_d01['data_LTG3'][1113:1204]
weekly_LTG3_d01[364:455] = ds_d01['data_LTG3'][1478:1569]
weekly_LTG3_d01[455:546] = ds_d01['data_LTG3'][1844:1935]
# # PR92_H DOMAIN 2
print('weekly average PR92H')
weekly_PR92_H_d02[0:91] = ds_d02['data_PR92H'][17:108]
weekly_PR92_H_d02[91:182] = ds_d02['data_PR92H'][383:474]
weekly_PR92_H_d02[182:273] = ds_d02['data_PR92H'][748:839]
weekly_PR92_H_d02[273:364] = ds_d02['data_PR92H'][1113:1204]
weekly_PR92_H_d02[364:455] = ds_d02['data_PR92H'][1478:1569]
weekly_PR92_H_d02[455:546] = ds_d02['data_PR92H'][1844:1935]
# PR92_H DOMAIN 1
weekly_PR92_H_d01[0:91] = ds_d01['data_PR92H'][17:108]
weekly_PR92_H_d01[91:182] = ds_d01['data_PR92H'][383:474]
weekly_PR92_H_d01[182:273] = ds_d01['data_PR92H'][748:839]
weekly_PR92_H_d01[273:364] = ds_d01['data_PR92H'][1113:1204]
weekly_PR92_H_d01[364:455] = ds_d01['data_PR92H'][1478:1569]
weekly_PR92_H_d01[455:546] = ds_d01['data_PR92H'][1844:1935]
# PR92_W DOMAIN 2
print('weekly average PR92W')
weekly_PR92_W_d02[0:91] = ds_d02['data_PR92W'][17:108]
weekly_PR92_W_d02[91:182] = ds_d02['data_PR92W'][383:474]
weekly_PR92_W_d02[182:273] = ds_d02['data_PR92W'][748:839]
weekly_PR92_W_d02[273:364] = ds_d02['data_PR92W'][1113:1204]
weekly_PR92_W_d02[364:455] = ds_d02['data_PR92W'][1478:1569]
weekly_PR92_W_d02[455:546] = ds_d02['data_PR92W'][1844:1935]
# PR92_W DOMAIN 1
weekly_PR92_W_d01[0:91] = ds_d01['data_PR92W'][17:108]
weekly_PR92_W_d01[91:182] = ds_d01['data_PR92W'][383:474]
weekly_PR92_W_d01[182:273] = ds_d01['data_PR92W'][748:839]
weekly_PR92_W_d01[273:364] = ds_d01['data_PR92W'][1113:1204]
weekly_PR92_W_d01[364:455] = ds_d01['data_PR92W'][1478:1569]
weekly_PR92_W_d01[455:546] = ds_d01['data_PR92W'][1844:1935]
# CAPExP DOMAIN 1
print('weekly average CAPExP')
weekly_CAPExP_d01[0:91] = ds_d01['data_CAPExP'][17:108]
weekly_CAPExP_d01[91:182] = ds_d01['data_CAPExP'][383:474]
weekly_CAPExP_d01[182:273] = ds_d01['data_CAPExP'][748:839]
weekly_CAPExP_d01[273:364] = ds_d01['data_CAPExP'][1113:1204]
weekly_CAPExP_d01[364:455] = ds_d01['data_CAPExP'][1478:1569]
weekly_CAPExP_d01[455:546] = ds_d01['data_CAPExP'][1844:1935]
# Obs DOMAIN 2
print('weekly average LTG')
weekly_Obs_d02[0:91] = ds_d02['data_Obs'][17:108]
weekly_Obs_d02[91:182] = ds_d02['data_Obs'][383:474]
weekly_Obs_d02[182:273] = ds_d02['data_Obs'][748:839]
weekly_Obs_d02[273:364] = ds_d02['data_Obs'][1113:1204]
weekly_Obs_d02[364:455] = ds_d02['data_Obs'][1478:1569]
weekly_Obs_d02[455:546] = ds_d02['data_Obs'][1844:1935]
# Obs DOMAIN 1
weekly_Obs_d01[0:91] = ds_d01['data_Obs'][17:108]
weekly_Obs_d01[91:182] = ds_d01['data_Obs'][383:474]
weekly_Obs_d01[182:273] = ds_d01['data_Obs'][748:839]
weekly_Obs_d01[273:364] = ds_d01['data_Obs'][1113:1204]
weekly_Obs_d01[364:455] = ds_d01['data_Obs'][1478:1569]
weekly_Obs_d01[455:546] = ds_d01['data_Obs'][1844:1935]

for i in range(0,78):
    print('for loop weekly'+str(i+1))
    j = (i+1)*7
    k=i*7
    ds_d02['LPI'][i] = np.nanmean(weekly_LPI_d02[k:j])
    ds_d01['LPI'][i] = np.nanmean(weekly_LPI_d01[k:j])
    ds_d01['LTG3'][i] = np.nanmean(weekly_LTG3_d01[k:j])
    ds_d02['LTG3'][i] = np.nanmean(weekly_LTG3_d02[k:j])
    ds_d01['PR92H'][i] = np.nanmean(weekly_PR92_H_d01[k:j])
    ds_d02['PR92H'][i] = np.nanmean(weekly_PR92_H_d02[k:j])
    ds_d01['PR92W'][i] = np.nanmean(weekly_PR92_W_d01[k:j])
    ds_d02['PR92W'][i] = np.nanmean(weekly_PR92_W_d02[k:j])
    ds_d01['CAPExP'][i] = np.nanmean(weekly_CAPExP_d01[k:j])
    ds_d01['Obs'][i] = np.nanmean(weekly_Obs_d01[k:j])
    ds_d02['Obs'][i] = np.nanmean(weekly_Obs_d02[k:j])
ds_d01.close()
ds_d02.close()
