# ---------------------------------------------------------------------------------------------
# MODULES
# ---------------------------------------------------------------------------------------------
import netCDF4 as nc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from netCDF4 import Dataset
import os
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import LogNorm, from_levels_and_colors, ListedColormap, LinearSegmentedColormap
from matplotlib import cm
import matplotlib.dates as mdates
from datetime import timedelta

# ---------------------------------------------------------------------------------------------
# CREATE NC FILE
# ---------------------------------------------------------------------------------------------
ds_mccaul_2015_d01 = Dataset('/scratch/leuven/336/vsc33651/nu-wrf-dev/wrfout_nc_files/wrfout_d01_2015_v03.nc','r')
inds_lat = np.unique(np.where((ds_mccaul_2015_d01['lat'][:] > 58) & (ds_mccaul_2015_d01['lat'][:] <62))[0])
inds_lon = np.unique(np.where((ds_mccaul_2015_d01['lon'][:] > -123) & (ds_mccaul_2015_d01['lon'][:] <-108))[1])
lat = ds_mccaul_2015_d01['lat'][inds_lat,inds_lon]
lon = ds_mccaul_2015_d01['lon'][inds_lat,inds_lon]
ds = Dataset('/scratch/leuven/317/vsc31786/nu-wrf-dev/domain1_all_mb.nc', mode='w', format='NETCDF4')
ds.createDimension('time', 13248)
ds.createDimension('lat', 63)
ds.createDimension('lon', 109)
ds.createVariable('time','int', dimensions=('time',),zlib=True)
ds.createVariable('lat',lat.dtype, dimensions=('lat','lon',),zlib=True)
ds.variables['lat'][:] = lat
ds.createVariable('lon',lon.dtype, dimensions=('lat','lon',),zlib=True)
ds.variables['lon'][:] = lon
ds.createVariable('LPI', 'f4', dimensions=('time','lat','lon',), zlib=True)
ds.createVariable('LTG3', 'f4', dimensions=('time','lat','lon',), zlib=True)
ds.createVariable('PR92H', 'f4', dimensions=('time','lat','lon',), zlib=True)
ds.createVariable('PR92W', 'f4', dimensions=('time','lat','lon',), zlib=True)
ds.createVariable('CAPExP', 'f4', dimensions=('time','lat','lon',), zlib=True)
ds.createVariable('W', 'f4', dimensions=('time','lat','lon',), zlib=True)
ds.createVariable('H', 'f4', dimensions=('time','lat','lon',), zlib=True)

# ---------------------------------------------------------------------------------------------
# LOAD DATA OF D01
# ---------------------------------------------------------------------------------------------
lis = Dataset('/scratch/leuven/336/vsc33651/nu-wrf-dev/Great_Slave_Lake/2015/lis_input.d01.nc')
height = lis['ELEVATION'][inds_lat,inds_lon]
## 2015
# LPI
ds_lpi_2015_d01 = Dataset('/scratch/leuven/336/vsc33651/nu-wrf-dev/wrfout_nc_files/lpi_time_d01_2015_v03.nc','r')
ds['LPI'][0:2208] = ds_lpi_2015_d01['LPI'][:,inds_lat, inds_lon]
ds['time'][0:2208] = ds_lpi_2015_d01['time'][:]
# McCaul
ds['LTG3'][0:2208] = np.multiply(ds_mccaul_2015_d01.variables['LTG3_MAX'][0:-1,inds_lat, inds_lon],12)
ds['W'][0:2208] = ds_mccaul_2015_d01.variables['W_UP_MAX'][0:-1,inds_lat, inds_lon]
# PR92
ds_upp_2015_d01 = Dataset('/scratch/leuven/336/vsc33651/nu-wrf-dev/Great_Slave_Lake/2015/WRFPRS_d01_2015.nc', mode='r')
c_d01 = np.multiply(0.97241, np.exp(0.048203*0.01548)) # TAKE CALIBRATION FACTOR FOR RESOLUTION INTO ACCOUNT
Ctop_2015 = np.subtract(ds_upp_2015_d01['Ctop'][408:,inds_lat, inds_lon].data, height) # Get top height rel. to ground
Ctop_2015_km = np.divide(Ctop_2015, 1000) # get top height in km
Ctop_2015_km[Ctop_2015_km < 0.0] = 0.0
ds['H'][0:2208,:,:] = Ctop_2015_km[:]
power_2015 = np.power(Ctop_2015_km, 4.9) # power for PR92 scheme
PR92H_2015 = np.multiply(0.0000344, power_2015) # Formula
PR92_H_2015_ph = np.multiply(PR92H_2015, 60) # To get flashes per hour instead of minute
ctopp_2015 = np.where((ds_mccaul_2015_d01.variables['CTOP2D'][0:-1,inds_lat, inds_lon] <440)
                      & (ds_mccaul_2015_d01.variables['CTOP2D'][0:-1,inds_lat, inds_lon] >0),1,0) # conv cloud
COD_2015 = np.where(ds_mccaul_2015_d01.variables['COD2D'][0:-1,inds_lat, inds_lon]>23,1,0) # conv cloud
conv_2015 = np.multiply(ctopp_2015, COD_2015) # to locate convective clouds
PR92_H_2015 = np.multiply(PR92_H_2015_ph, conv_2015) # to only get conv clouds
PR92H_2015_spat = np.multiply(PR92_H_2015, c_d01) # apply spatial calibration factor (PR94)
PR92H_2015_deg = np.multiply(PR92H_2015_spat,0.0625) # apply custom calibration factor
print(np.nanmax(PR92H_2015_deg))
ds['PR92H'][0:2208] = PR92H_2015_deg
print(np.nanmax(ds['PR92H'][:]))

power_w_2015 = np.power(ds_mccaul_2015_d01.variables['W_UP_MAX'][0:-1,inds_lat, inds_lon], 4.54)
PR92_W_2015 = np.multiply(0.000005, power_w_2015)
PR92_w_2015_ph = np.multiply(PR92_W_2015, 60)
PR92_w_2015_c = np.multiply(PR92_w_2015_ph, c_d01)
PR92W_2015_deg = np.multiply(PR92_w_2015_c,0.0625)
ds['PR92W'][0:2208] = PR92W_2015_deg
print(np.nanmax(ds['PR92H'][:]))
# 2016
# LPI
ds_lpi_2016_d01 = Dataset('/scratch/leuven/336/vsc33651/nu-wrf-dev/wrfout_nc_files/lpi_time_d01_2016_v03.nc','r')
ds['LPI'][2208:4416] = ds_lpi_2016_d01['LPI'][:,inds_lat, inds_lon]
ds['time'][2208:4416] = ds_lpi_2016_d01['time'][:]

# McCaul
ds_mccaul_2016_d01 = Dataset('/scratch/leuven/336/vsc33651/nu-wrf-dev/wrfout_nc_files/wrfout_d01_2016_v03.nc','r')
ds['LTG3'][2208:4416] = np.multiply(ds_mccaul_2016_d01.variables['LTG3_MAX'][0:-1,inds_lat, inds_lon],12)
ds['W'][2208:4416] = ds_mccaul_2016_d01.variables['W_UP_MAX'][0:-1,inds_lat, inds_lon]
# PR92
ds_upp_2016_d01 = Dataset('/scratch/leuven/336/vsc33651/nu-wrf-dev/Great_Slave_Lake/2016/WRFPRS_d02_2016.nc', mode='r')
Ctop_2016 = np.subtract(ds_upp_2016_d01['Ctop'][408:,inds_lat, inds_lon].data, height) # Get top height rel. to ground
Ctop_2016_km = np.divide(Ctop_2016, 1000) # get top height in km
Ctop_2016_km[Ctop_2016_km < 0.0] = 0.0
ds['H'][2208:4416] = Ctop_2016_km[:]
power_2016 = np.power(Ctop_2016_km, 4.9) # power for PR92 scheme
PR92H_2016 = np.multiply(0.0000344, power_2016) # Formula
PR92_H_2016_ph = np.multiply(PR92H_2016, 60) # To get flashes per hour instead of minute
ctopp_2016 = np.where((ds_mccaul_2016_d01.variables['CTOP2D'][0:-1,inds_lat, inds_lon] <440)
                      & (ds_mccaul_2016_d01.variables['CTOP2D'][0:-1,inds_lat, inds_lon] >0),1,0) # conv cloud
COD_2016 = np.where(ds_mccaul_2016_d01.variables['COD2D'][0:-1,inds_lat, inds_lon]>23,1,0) # conv cloud
conv_2016 = np.multiply(ctopp_2016, COD_2016) # to locate convective clouds
PR92_H_2016 = np.multiply(PR92_H_2016_ph, conv_2016) # to only get conv clouds
PR92H_2016_spat = np.multiply(PR92_H_2016, c_d01) # apply spatial calibration factor (PR94)
PR92H_2016_deg = np.multiply(PR92H_2016_spat,0.0625) # apply custom calibration factor
print(np.nanmax(PR92H_2016_deg))
ds['PR92H'][2208:4416] = PR92H_2016_deg
print(np.nanmax(ds['PR92H'][:]))

power_w_2016 = np.power(ds_mccaul_2016_d01.variables['W_UP_MAX'][0:-1,inds_lat, inds_lon], 4.54)
PR92_W_2016 = np.multiply(0.000005, power_w_2016)
PR92_w_2016_ph = np.multiply(PR92_W_2016, 60)
PR92_w_2016_c = np.multiply(PR92_w_2016_ph, c_d01)
PR92W_2016_deg = np.multiply(PR92_w_2016_c,0.0625)
ds['PR92W'][2208:4416] = PR92W_2016_deg
print(np.nanmax(ds['PR92H'][:]))
# 2017
# LPI
ds_lpi_2017_d01 = Dataset('/scratch/leuven/336/vsc33651/nu-wrf-dev/wrfout_nc_files/lpi_time_d01_2017_v03.nc','r')
ds['time'][4416:6624] = ds_lpi_2017_d01['time'][:]
ds['LPI'][4416:6624] = ds_lpi_2017_d01['LPI'][:,inds_lat, inds_lon]
# McCaul
ds_mccaul_2017_d01 = Dataset('/scratch/leuven/336/vsc33651/nu-wrf-dev/wrfout_nc_files/wrfout_d01_2017_v03.nc','r')
ds['LTG3'][4416:6624] = np.multiply(ds_mccaul_2017_d01.variables['LTG3_MAX'][0:-1,inds_lat, inds_lon],12)
ds['W'][4416:6624] = ds_mccaul_2017_d01.variables['W_UP_MAX'][0:-1,inds_lat, inds_lon]
# PR92
ds_upp_2017_d01 = Dataset('/scratch/leuven/336/vsc33651/nu-wrf-dev/Great_Slave_Lake/2017/WRFPRS_d01_2017.nc', mode='r')
Ctop_2017 = np.subtract(ds_upp_2017_d01['Ctop'][408:,inds_lat, inds_lon].data, height) # Get top height rel. to ground
Ctop_2017_km = np.divide(Ctop_2017, 1000) # get top height in km
Ctop_2017_km[Ctop_2017_km < 0.0] = 0.0
ds['H'][4416:6624] = Ctop_2017_km[:]
power_2017 = np.power(Ctop_2017_km, 4.9) # power for PR92 scheme
PR92H_2017 = np.multiply(0.0000344, power_2017) # Formula
PR92_H_2017_ph = np.multiply(PR92H_2017, 60) # To get flashes per hour instead of minute
ctopp_2017 = np.where((ds_mccaul_2017_d01.variables['CTOP2D'][0:-1,inds_lat, inds_lon] <440)
                      & (ds_mccaul_2017_d01.variables['CTOP2D'][0:-1,inds_lat, inds_lon] >0),1,0) # conv cloud
COD_2017 = np.where(ds_mccaul_2017_d01.variables['COD2D'][0:-1,inds_lat, inds_lon]>23,1,0) # conv cloud
conv_2017 = np.multiply(ctopp_2017, COD_2017) # to locate convective clouds
PR92_H_2017 = np.multiply(PR92_H_2017_ph, conv_2017) # to only get conv clouds
PR92H_2017_spat = np.multiply(PR92_H_2017, c_d01) # apply spatial calibration factor (PR94)
PR92H_2017_deg = np.multiply(PR92H_2017_spat,0.0625) # apply custom calibration factor
print(np.nanmax(PR92H_2017_deg))
ds['PR92H'][4416:6624] = PR92H_2017_deg
print(np.nanmax(ds['PR92H'][:]))

power_w_2017 = np.power(ds_mccaul_2017_d01.variables['W_UP_MAX'][0:-1,inds_lat, inds_lon], 4.54)
PR92_W_2017 = np.multiply(0.000005, power_w_2017)
PR92_w_2017_ph = np.multiply(PR92_W_2017, 60)
PR92_w_2017_c = np.multiply(PR92_w_2017_ph, c_d01)
PR92W_2017_deg = np.multiply(PR92_w_2017_c,0.0625)
ds['PR92W'][4416:6624] = PR92W_2017_deg
print(np.nanmax(ds['PR92H'][:]))
# 2018
# LPI
ds_lpi_2018_d01 = Dataset('/scratch/leuven/336/vsc33651/nu-wrf-dev/wrfout_nc_files/lpi_time_d01_2018_v03.nc','r')
ds['time'][6624:8832] = ds_lpi_2018_d01['time'][:]
ds['LPI'][6624:8832] = ds_lpi_2018_d01['LPI'][:,inds_lat, inds_lon]
# McCaul
ds_mccaul_2018_d01 = Dataset('/scratch/leuven/336/vsc33651/nu-wrf-dev/wrfout_nc_files/wrfout_d01_2018_v03.nc','r')
ds['LTG3'][6624:8832] = np.multiply(ds_mccaul_2018_d01.variables['LTG3_MAX'][0:-1,inds_lat, inds_lon],12)
ds['W'][6624:8832] = ds_mccaul_2018_d01.variables['W_UP_MAX'][0:-1,inds_lat, inds_lon]
# PR92
ds_upp_2018_d01 = Dataset('/scratch/leuven/336/vsc33651/nu-wrf-dev/Great_Slave_Lake/2018/WRFPRS_d01_2018.nc', mode='r')
Ctop_2018 = np.subtract(ds_upp_2018_d01['Ctop'][408:,inds_lat, inds_lon].data, height) # Get top height rel. to ground
Ctop_2018_km = np.divide(Ctop_2018, 1000) # get top height in km
Ctop_2018_km[Ctop_2018_km < 0.0] = 0.0
ds['H'][6624:8832] = Ctop_2018_km[:]
power_2018 = np.power(Ctop_2018_km, 4.9) # power for PR92 scheme
PR92H_2018 = np.multiply(0.0000344, power_2018) # Formula
PR92_H_2018_ph = np.multiply(PR92H_2018, 60) # To get flashes per hour instead of minute
ctopp_2018 = np.where((ds_mccaul_2018_d01.variables['CTOP2D'][0:-1,inds_lat, inds_lon] <440)
                      & (ds_mccaul_2018_d01.variables['CTOP2D'][0:-1,inds_lat, inds_lon] >0),1,0) # conv cloud
COD_2018 = np.where(ds_mccaul_2018_d01.variables['COD2D'][0:-1,inds_lat, inds_lon]>23,1,0) # conv cloud
conv_2018 = np.multiply(ctopp_2018, COD_2018) # to locate convective clouds
PR92_H_2018 = np.multiply(PR92_H_2018_ph, conv_2018) # to only get conv clouds
PR92H_2018_spat = np.multiply(PR92_H_2018, c_d01) # apply spatial calibration factor (PR94)
PR92H_2018_deg = np.multiply(PR92H_2018_spat,0.0625) # apply custom calibration factor
print(np.nanmax(PR92H_2018_deg))
ds['PR92H'][6624:8832] = PR92H_2018_deg
print(np.nanmax(ds['PR92H'][:]))

power_w_2018 = np.power(ds_mccaul_2018_d01.variables['W_UP_MAX'][0:-1,inds_lat, inds_lon], 4.54)
PR92_W_2018 = np.multiply(0.000005, power_w_2018)
PR92_w_2018_ph = np.multiply(PR92_W_2018, 60)
PR92_w_2018_c = np.multiply(PR92_w_2018_ph, c_d01)
PR92W_2018_deg = np.multiply(PR92_w_2018_c,0.0625)
ds['PR92W'][6624:8832] = PR92W_2018_deg
print(np.nanmax(ds['PR92H'][:]))
# 2019
# LPI
ds_lpi_2019_d01 = Dataset('/scratch/leuven/336/vsc33651/nu-wrf-dev/wrfout_nc_files/lpi_time_d01_2019_v03.nc','r')
ds['time'][8832:11040] = ds_lpi_2019_d01['time'][:]
ds['LPI'][8832:11040] = ds_lpi_2019_d01['LPI'][:,inds_lat, inds_lon]
# McCaul
ds_mccaul_2019_d01 = Dataset('/scratch/leuven/336/vsc33651/nu-wrf-dev/wrfout_nc_files/wrfout_d01_2019_v03.nc','r')
ds['LTG3'][8832:11040] = np.multiply(ds_mccaul_2019_d01.variables['LTG3_MAX'][0:-1,inds_lat, inds_lon],12)
ds['W'][8832:11040] = ds_mccaul_2019_d01.variables['W_UP_MAX'][0:-1,inds_lat, inds_lon]
# PR92
ds_upp_2019_d01 = Dataset('/scratch/leuven/336/vsc33651/nu-wrf-dev/Great_Slave_Lake/2019/WRFPRS_d01_2019.nc', mode='r')
Ctop_2019 = np.subtract(ds_upp_2019_d01['Ctop'][408:,inds_lat, inds_lon].data, height) # Get top height rel. to ground
Ctop_2019_km = np.divide(Ctop_2019, 1000) # get top height in km
Ctop_2019_km[Ctop_2019_km < 0.0] = 0.0
ds['H'][8832:11040] = Ctop_2019_km[:]
power_2019 = np.power(Ctop_2019_km, 4.9) # power for PR92 scheme
PR92H_2019 = np.multiply(0.0000344, power_2019) # Formula
PR92_H_2019_ph = np.multiply(PR92H_2019, 60) # To get flashes per hour instead of minute
ctopp_2019 = np.where((ds_mccaul_2019_d01.variables['CTOP2D'][0:-1,inds_lat, inds_lon] <440)
                      & (ds_mccaul_2019_d01.variables['CTOP2D'][0:-1,inds_lat, inds_lon] >0),1,0) # conv cloud
COD_2019 = np.where(ds_mccaul_2019_d01.variables['COD2D'][0:-1,inds_lat, inds_lon]>23,1,0) # conv cloud
conv_2019 = np.multiply(ctopp_2019, COD_2019) # to locate convective clouds
PR92_H_2019 = np.multiply(PR92_H_2019_ph, conv_2019) # to only get conv clouds
PR92H_2019_spat = np.multiply(PR92_H_2019, c_d01) # apply spatial calibration factor (PR94)
PR92H_2019_deg = np.multiply(PR92H_2019_spat,0.0625) # apply custom calibration factor
print(np.nanmax(PR92H_2019_deg))
ds['PR92H'][8832:11040] = PR92H_2019_deg
print(np.nanmax(ds['PR92H'][:]))

power_w_2019 = np.power(ds_mccaul_2019_d01.variables['W_UP_MAX'][0:-1,inds_lat, inds_lon], 4.54)
PR92_W_2019 = np.multiply(0.000005, power_w_2019)
PR92_w_2019_ph = np.multiply(PR92_W_2019, 60)
PR92_w_2019_c = np.multiply(PR92_w_2019_ph, c_d01)
PR92W_2019_deg = np.multiply(PR92_w_2019_c,0.0625)
ds['PR92W'][8832:11040] = PR92W_2019_deg
print(np.nanmax(ds['PR92H'][:]))
# 2020
# LPI
ds_lpi_2020_d01 = Dataset('/scratch/leuven/336/vsc33651/nu-wrf-dev/wrfout_nc_files/lpi_time_d01_2020_v03.nc','r')
ds['time'][11040:13248] = ds_lpi_2020_d01['time'][:]
ds['LPI'][11040:13248] = ds_lpi_2020_d01['LPI'][:,inds_lat, inds_lon]
# McCaul
ds_mccaul_2020_d01 = Dataset('/scratch/leuven/336/vsc33651/nu-wrf-dev/wrfout_nc_files/wrfout_d01_2020_v03.nc','r')
ds['LTG3'][11040:13248] = np.multiply(ds_mccaul_2020_d01.variables['LTG3_MAX'][0:-1,inds_lat, inds_lon],12)
ds['W'][11040:13248] = ds_mccaul_2020_d01.variables['W_UP_MAX'][0:-1,inds_lat, inds_lon]
# PR92
ds_upp_2020_d01 = Dataset('/scratch/leuven/336/vsc33651/nu-wrf-dev/Great_Slave_Lake/2020/WRFPRS_d01_2020.nc', mode='r')
Ctop_2020 = np.subtract(ds_upp_2020_d01['Ctop'][408:,inds_lat, inds_lon].data, height) # Get top height rel. to ground
Ctop_2020_km = np.divide(Ctop_2020, 1000) # get top height in km
Ctop_2020_km[Ctop_2020_km < 0.0] = 0.0
ds['H'][11040:13248] = Ctop_2020_km[:]
power_2020 = np.power(Ctop_2020_km, 4.9) # power for PR92 scheme
PR92H_2020 = np.multiply(0.0000344, power_2020) # Formula
PR92_H_2020_ph = np.multiply(PR92H_2020, 60) # To get flashes per hour instead of minute
ctopp_2020 = np.where((ds_mccaul_2020_d01.variables['CTOP2D'][0:-1,inds_lat, inds_lon] <440)
                      & (ds_mccaul_2020_d01.variables['CTOP2D'][0:-1,inds_lat, inds_lon] >0),1,0) # conv cloud
COD_2020 = np.where(ds_mccaul_2020_d01.variables['COD2D'][0:-1,inds_lat, inds_lon]>23,1,0) # conv cloud
conv_2020 = np.multiply(ctopp_2020, COD_2020) # to locate convective clouds
PR92_H_2020 = np.multiply(PR92_H_2020_ph, conv_2020) # to only get conv clouds
PR92H_2020_spat = np.multiply(PR92_H_2020, c_d01) # apply spatial calibration factor (PR94)
PR92H_2020_deg = np.multiply(PR92H_2020_spat, 0.0625) # apply custom calibration factor
print(np.nanmax(PR92H_2020_deg))
ds['PR92H'][11040:13248] = PR92H_2020_deg
print(np.nanmax(ds['PR92H'][:]))

power_w_2020 = np.power(ds_mccaul_2020_d01.variables['W_UP_MAX'][0:-1,inds_lat, inds_lon], 4.54)
PR92_W_2020 = np.multiply(0.000005, power_w_2020)
PR92_w_2020_ph = np.multiply(PR92_W_2020, 60)
PR92_w_2020_c = np.multiply(PR92_w_2020_ph, c_d01)
PR92W_2020_deg = np.multiply(PR92_w_2020_c,0.0625)
ds['PR92W'][11040:13248] = PR92W_2020_deg
print(np.nanmax(ds['PR92H'][:]))

# CAPExP
print('CAPExP')
RAINC = np.zeros((13248,63,109))
RAINC[0,:,:] = ds_mccaul_2015_d01['RAINC'][0,inds_lat,inds_lon]
RAINC[2208,:,:] = ds_mccaul_2016_d01['RAINC'][0,inds_lat,inds_lon]
RAINC[4416,:,:] = ds_mccaul_2017_d01['RAINC'][0,inds_lat,inds_lon]
RAINC[6624,:,:] = ds_mccaul_2018_d01['RAINC'][0,inds_lat,inds_lon]
RAINC[8832,:,:] = ds_mccaul_2019_d01['RAINC'][0,inds_lat,inds_lon]
RAINC[11040,:,:] = ds_mccaul_2020_d01['RAINC'][0,inds_lat,inds_lon]
CAPE = np.zeros((13248,63,109))
CAPE[0:2208,:,:] = ds_upp_2015_d01['CAPE'][408:,inds_lat,inds_lon]
CAPE[2208:4416,:,:] = ds_upp_2016_d01['CAPE'][408:,inds_lat,inds_lon]
CAPE[4416:6624,:,:] = ds_upp_2017_d01['CAPE'][408:,inds_lat,inds_lon]
CAPE[6624:8832,:,:] = ds_upp_2018_d01['CAPE'][408:,inds_lat,inds_lon]
CAPE[8832:11040,:,:] = ds_upp_2019_d01['CAPE'][408:,inds_lat,inds_lon]
CAPE[11040:13248,:,:] = ds_upp_2020_d01['CAPE'][408:,inds_lat,inds_lon]

for i in range(1, 2208):
    print('loop 1/6 ' + str(i))
    RAINC[i,:,:] = np.subtract(ds_mccaul_2015_d01['RAINC'][i,inds_lat,inds_lon], ds_mccaul_2015_d01['RAINC'][i-1,inds_lat, inds_lon])
for i in range(1,2208):
    print('loop 2/6 ' + str(i))
    RAINC[i+2208, :, :] = np.subtract(ds_mccaul_2016_d01['RAINC'][i,inds_lat, inds_lon], ds_mccaul_2016_d01['RAINC'][i-1,inds_lat, inds_lon])
for i in range(1,2208):
    print('loop 3/6 ' + str(i))
    RAINC[i+4416,:,:] = np.subtract(ds_mccaul_2017_d01['RAINC'][i,inds_lat, inds_lon], ds_mccaul_2017_d01['RAINC'][i-1,inds_lat, inds_lon])
for i in range(1,2208):
    print('loop 4/6 ' + str(i))
    RAINC[i+6624,:,:] = np.subtract(ds_mccaul_2018_d01['RAINC'][i,inds_lat, inds_lon], ds_mccaul_2018_d01['RAINC'][i-1,inds_lat, inds_lon])
for i in range(1,2208):
    print('loop 5/6 ' + str(i))
    RAINC[i+8832,:,:] = np.subtract(ds_mccaul_2019_d01['RAINC'][i,inds_lat, inds_lon], ds_mccaul_2019_d01['RAINC'][i-1,inds_lat, inds_lon])
for i in range(1,2208):
    print('loop 6/6 ' + str(i))
    RAINC[i+11040,:,:] = np.subtract(ds_mccaul_2020_d01['RAINC'][i,inds_lat, inds_lon], ds_mccaul_2020_d01['RAINC'][i-1,inds_lat, inds_lon])

print('CAPExP')
RAINC[:] = np.divide(RAINC[:],3600) # To convert units to km/m²/s
ds['CAPExP'][:] = np.multiply(CAPE[:], RAINC[:])
ds['CAPExP'][:] = np.multiply(ds['CAPExP'][:], 0.000000000013) # constant, see Romps et al. 2014 (eta/E)
ds['CAPExP'][:] = np.multiply(ds['CAPExP'][:], 3600000000) # To convert from (m2s)-1 to (km2hour)-1

#plt.plot(ds['H'][:,50,50],'.')

ds.close()
