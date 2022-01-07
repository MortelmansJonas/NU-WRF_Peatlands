#!/usr/bin/env python
# ---------------------------------------------------------------------------------------------
# MODULES
# ---------------------------------------------------------------------------------------------
import numpy as np
from netCDF4 import Dataset

# ---------------------------------------------------------------------------------------------
# CREATE NC FILE
# ---------------------------------------------------------------------------------------------
ds_mccaul_2015_d01 = Dataset('/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/wrfout_nc_files/wrfout_d01_2015.nc','r')
inds_lat = np.unique(np.where((ds_mccaul_2015_d01['lat'][:] > 58) & (ds_mccaul_2015_d01['lat'][:] <62))[0])
inds_lon = np.unique(np.where((ds_mccaul_2015_d01['lon'][:] > -123) & (ds_mccaul_2015_d01['lon'][:] <-108))[1])
lat = ds_mccaul_2015_d01['lat'][inds_lat,inds_lon]
lon = ds_mccaul_2015_d01['lon'][inds_lat,inds_lon]
ds = Dataset('/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/wrfout_nc_files/d01_all.nc', mode='w', format='NETCDF4')
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
ds.createVariable('PR92W', 'f4', dimensions=('time','lat','lon',), zlib=True)
ds.createVariable('CAPExP_R', 'f4', dimensions=('time','lat','lon',), zlib=True)
ds.createVariable('CAPExP_CSI', 'f4', dimensions=('time','lat','lon',), zlib=True)
ds.createVariable('W', 'f4', dimensions=('time','lat','lon',), zlib=True)
ds.createVariable('T2', 'f4', dimensions=('time','lat','lon',), zlib=True)
ds.createVariable('RH2', 'f4', dimensions=('time','lat','lon',), zlib=True)
ds.createVariable('LH', 'f4', dimensions=('time','lat','lon',), zlib=True)
ds.createVariable('SH', 'f4', dimensions=('time','lat','lon',), zlib=True)
ds.createVariable('RAINC', 'f4', dimensions=('time','lat','lon',), zlib=True)
ds.createVariable('RAINNC', 'f4', dimensions=('time','lat','lon',), zlib=True)

# ---------------------------------------------------------------------------------------------
# LOAD DATA OF D01
# ---------------------------------------------------------------------------------------------
## 2015
# LPI
ds_lpi_2015_d01 = Dataset('/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/wrfout_nc_files/lpi_d01_2015.nc','r')
ds['LPI'][0:2208] = ds_lpi_2015_d01['LPI'][0:-1,inds_lat, inds_lon]
ds['time'][0:2208] = ds_lpi_2015_d01['time'][0:-1]
# McCaul
ds['LTG3'][0:2208] = np.multiply(ds_mccaul_2015_d01.variables['LTG3_MAX'][0:-1,inds_lat, inds_lon],12)
ds['W'][0:2208] = ds_mccaul_2015_d01.variables['W_UP_MAX'][0:-1,inds_lat, inds_lon]
# PR92
c_d01 = np.multiply(0.97241, np.exp(0.048203*0.01548)) # TAKE CALIBRATION FACTOR FOR RESOLUTION INTO ACCOUNT
power_w_2015 = np.power(ds_mccaul_2015_d01.variables['W_UP_MAX'][0:-1,inds_lat, inds_lon], 4.54)
PR92_W_2015 = np.multiply(0.000005, power_w_2015)
PR92_w_2015_ph = np.multiply(PR92_W_2015, 60)
PR92_w_2015_c = np.multiply(PR92_w_2015_ph, c_d01)
PR92W_2015_deg = np.multiply(PR92_w_2015_c,0.0625)
ds['PR92W'][0:2208] = PR92W_2015_deg
# Atmospheric variables
ds['LH'][0:2208] = ds_mccaul_2015_d01.variables['LH'][:-1, inds_lat, inds_lon]
ds['SH'][0:2208] = ds_mccaul_2015_d01.variables['SH'][:-1, inds_lat, inds_lon]
ds['T2'][0:2208] = ds_mccaul_2015_d01.variables['T2'][:-1, inds_lat, inds_lon]
ds['RH2'][0:2208] = ds_mccaul_2015_d01.variables['RH2'][:-1, inds_lat, inds_lon]


## 2016
# LPI
ds_lpi_2016_d01 = Dataset('/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/wrfout_nc_files/lpi_d01_2016.nc','r')
ds['LPI'][2208:4416] = ds_lpi_2016_d01['LPI'][0:-1,inds_lat, inds_lon]
ds['time'][2208:4416] = ds_lpi_2016_d01['time'][0:-1]

# McCaul
ds_mccaul_2016_d01 = Dataset('/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/wrfout_nc_files/wrfout_d01_2016.nc','r')
ds['LTG3'][2208:4416] = np.multiply(ds_mccaul_2016_d01.variables['LTG3_MAX'][0:-1,inds_lat, inds_lon],12)
ds['W'][2208:4416] = ds_mccaul_2016_d01.variables['W_UP_MAX'][0:-1,inds_lat, inds_lon]
# PR92
power_w_2016 = np.power(ds_mccaul_2016_d01.variables['W_UP_MAX'][0:-1,inds_lat, inds_lon], 4.54)
PR92_W_2016 = np.multiply(0.000005, power_w_2016)
PR92_w_2016_ph = np.multiply(PR92_W_2016, 60)
PR92_w_2016_c = np.multiply(PR92_w_2016_ph, c_d01)
PR92W_2016_deg = np.multiply(PR92_w_2016_c,0.0625)
ds['PR92W'][2208:4416] = PR92W_2016_deg

# Atmospheric variables
ds['LH'][2208:4416] = ds_mccaul_2016_d01.variables['LH'][:-1, inds_lat, inds_lon]
ds['SH'][2208:4416] = ds_mccaul_2016_d01.variables['SH'][:-1, inds_lat, inds_lon]
ds['T2'][2208:4416] = ds_mccaul_2016_d01.variables['T2'][:-1, inds_lat, inds_lon]
ds['RH2'][2208:4416] = ds_mccaul_2016_d01.variables['RH2'][:-1, inds_lat, inds_lon]


# 2017
# LPI
ds_lpi_2017_d01 = Dataset('/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/wrfout_nc_files/lpi_d01_2017.nc','r')
ds['time'][4416:6624] = ds_lpi_2017_d01['time'][0:-1]
ds['LPI'][4416:6624] = ds_lpi_2017_d01['LPI'][0:-1,inds_lat, inds_lon]
# McCaul
ds_mccaul_2017_d01 = Dataset('/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/wrfout_nc_files/wrfout_d01_2017.nc','r')
ds['LTG3'][4416:6624] = np.multiply(ds_mccaul_2017_d01.variables['LTG3_MAX'][0:-1,inds_lat, inds_lon],12)
ds['W'][4416:6624] = ds_mccaul_2017_d01.variables['W_UP_MAX'][0:-1,inds_lat, inds_lon]
# PR92
power_w_2017 = np.power(ds_mccaul_2017_d01.variables['W_UP_MAX'][0:-1,inds_lat, inds_lon], 4.54)
PR92_W_2017 = np.multiply(0.000005, power_w_2017)
PR92_w_2017_ph = np.multiply(PR92_W_2017, 60)
PR92_w_2017_c = np.multiply(PR92_w_2017_ph, c_d01)
PR92W_2017_deg = np.multiply(PR92_w_2017_c,0.0625)
ds['PR92W'][4416:6624] = PR92W_2017_deg

# Atmospheric variables
ds['LH'][4416:6624] = ds_mccaul_2017_d01.variables['LH'][:-1, inds_lat, inds_lon]
ds['SH'][4416:6624] = ds_mccaul_2017_d01.variables['SH'][:-1, inds_lat, inds_lon]
ds['T2'][4416:6624] = ds_mccaul_2017_d01.variables['T2'][:-1, inds_lat, inds_lon]
ds['RH2'][4416:6624] = ds_mccaul_2017_d01.variables['RH2'][:-1, inds_lat, inds_lon]

# 2018
# LPI
ds_lpi_2018_d01 = Dataset('/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/wrfout_nc_files/lpi_d01_2018.nc','r')
ds['time'][6624:8832] = ds_lpi_2018_d01['time'][0:-1]
ds['LPI'][6624:8832] = ds_lpi_2018_d01['LPI'][0:-1,inds_lat, inds_lon]
# McCaul
ds_mccaul_2018_d01 = Dataset('/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/wrfout_nc_files/wrfout_d01_2018.nc','r')
ds['LTG3'][6624:8832] = np.multiply(ds_mccaul_2018_d01.variables['LTG3_MAX'][0:-1,inds_lat, inds_lon],12)
ds['W'][6624:8832] = ds_mccaul_2018_d01.variables['W_UP_MAX'][0:-1,inds_lat, inds_lon]
# PR92
power_w_2018 = np.power(ds_mccaul_2018_d01.variables['W_UP_MAX'][0:-1,inds_lat, inds_lon], 4.54)
PR92_W_2018 = np.multiply(0.000005, power_w_2018)
PR92_w_2018_ph = np.multiply(PR92_W_2018, 60)
PR92_w_2018_c = np.multiply(PR92_w_2018_ph, c_d01)
PR92W_2018_deg = np.multiply(PR92_w_2018_c,0.0625)
ds['PR92W'][6624:8832] = PR92W_2018_deg

# Atmospheric variables
ds['LH'][6624:8832] = ds_mccaul_2018_d01.variables['LH'][:-1, inds_lat, inds_lon]
ds['SH'][6624:8832] = ds_mccaul_2018_d01.variables['SH'][:-1, inds_lat, inds_lon]
ds['T2'][6624:8832] = ds_mccaul_2018_d01.variables['T2'][:-1, inds_lat, inds_lon]
ds['RH2'][6624:8832] = ds_mccaul_2018_d01.variables['RH2'][:-1, inds_lat, inds_lon]

# 2019
# LPI
ds_lpi_2019_d01 = Dataset('/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/wrfout_nc_files/lpi_d01_2019.nc','r')
ds['time'][8832:11040] = ds_lpi_2019_d01['time'][0:-1]
ds['LPI'][8832:11040] = ds_lpi_2019_d01['LPI'][0:-1,inds_lat, inds_lon]
# McCaul
ds_mccaul_2019_d01 = Dataset('/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/wrfout_nc_files/wrfout_d01_2019.nc','r')
ds['LTG3'][8832:11040] = np.multiply(ds_mccaul_2019_d01.variables['LTG3_MAX'][0:-1,inds_lat, inds_lon],12)
ds['W'][8832:11040] = ds_mccaul_2019_d01.variables['W_UP_MAX'][0:-1,inds_lat, inds_lon]
# PR92
power_w_2019 = np.power(ds_mccaul_2019_d01.variables['W_UP_MAX'][0:-1,inds_lat, inds_lon], 4.54)
PR92_W_2019 = np.multiply(0.000005, power_w_2019)
PR92_w_2019_ph = np.multiply(PR92_W_2019, 60)
PR92_w_2019_c = np.multiply(PR92_w_2019_ph, c_d01)
PR92W_2019_deg = np.multiply(PR92_w_2019_c,0.0625)
ds['PR92W'][8832:11040] = PR92W_2019_deg

# Atmospheric variables
ds['LH'][8832:11040] = ds_mccaul_2019_d01.variables['LH'][:-1, inds_lat, inds_lon]
ds['SH'][8832:11040] = ds_mccaul_2019_d01.variables['SH'][:-1, inds_lat, inds_lon]
ds['T2'][8832:11040] = ds_mccaul_2019_d01.variables['T2'][:-1, inds_lat, inds_lon]
ds['RH2'][8832:11040] = ds_mccaul_2019_d01.variables['RH2'][:-1, inds_lat, inds_lon]

# 2020
# LPI
ds_lpi_2020_d01 = Dataset('/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/wrfout_nc_files/lpi_d01_2020.nc','r')
ds['time'][11040:13248] = ds_lpi_2020_d01['time'][0:-1]
ds['LPI'][11040:13248] = ds_lpi_2020_d01['LPI'][0:-1,inds_lat, inds_lon]
# McCaul
ds_mccaul_2020_d01 = Dataset('/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/wrfout_nc_files/wrfout_d01_2020.nc','r')
ds['LTG3'][11040:13248] = np.multiply(ds_mccaul_2020_d01.variables['LTG3_MAX'][0:-1,inds_lat, inds_lon],12)
ds['W'][11040:13248] = ds_mccaul_2020_d01.variables['W_UP_MAX'][0:-1,inds_lat, inds_lon]
# PR92
power_w_2020 = np.power(ds_mccaul_2020_d01.variables['W_UP_MAX'][0:-1,inds_lat, inds_lon], 4.54)
PR92_W_2020 = np.multiply(0.000005, power_w_2020)
PR92_w_2020_ph = np.multiply(PR92_W_2020, 60)
PR92_w_2020_c = np.multiply(PR92_w_2020_ph, c_d01)
PR92W_2020_deg = np.multiply(PR92_w_2020_c,0.0625)
ds['PR92W'][11040:13248] = PR92W_2020_deg
# Atmospheric variables
ds['LH'][11040:13248] = ds_mccaul_2020_d01.variables['LH'][:-1, inds_lat, inds_lon]
ds['SH'][11040:13248] = ds_mccaul_2020_d01.variables['SH'][:-1, inds_lat, inds_lon]
ds['T2'][11040:13248] = ds_mccaul_2020_d01.variables['T2'][:-1, inds_lat, inds_lon]
ds['RH2'][11040:13248] = ds_mccaul_2020_d01.variables['RH2'][:-1, inds_lat, inds_lon]

# CAPExP
print('CAPExP')
RAINC = np.zeros((13248,63,109))
RAINC[0,:,:] = ds_mccaul_2015_d01['RAINC'][0,inds_lat,inds_lon]
RAINC[2208,:,:] = ds_mccaul_2016_d01['RAINC'][0,inds_lat,inds_lon]
RAINC[4416,:,:] = ds_mccaul_2017_d01['RAINC'][0,inds_lat,inds_lon]
RAINC[6624,:,:] = ds_mccaul_2018_d01['RAINC'][0,inds_lat,inds_lon]
RAINC[8832,:,:] = ds_mccaul_2019_d01['RAINC'][0,inds_lat,inds_lon]
RAINC[11040,:,:] = ds_mccaul_2020_d01['RAINC'][0,inds_lat,inds_lon]
RAINNC = np.zeros((13248,63,109))
RAINNC[0,:,:] = ds_mccaul_2015_d01['RAINNC'][0,inds_lat,inds_lon]
RAINNC[2208,:,:] = ds_mccaul_2016_d01['RAINNC'][0,inds_lat,inds_lon]
RAINNC[4416,:,:] = ds_mccaul_2017_d01['RAINNC'][0,inds_lat,inds_lon]
RAINNC[6624,:,:] = ds_mccaul_2018_d01['RAINNC'][0,inds_lat,inds_lon]
RAINNC[8832,:,:] = ds_mccaul_2019_d01['RAINNC'][0,inds_lat,inds_lon]
RAINNC[11040,:,:] = ds_mccaul_2020_d01['RAINNC'][0,inds_lat,inds_lon]
CAPE = np.zeros((13248,63,109))
CAPE[0:2208,:,:] = ds_mccaul_2015_d01['CAPE'][:-1,inds_lat,inds_lon]
CAPE[2208:4416,:,:] = ds_mccaul_2016_d01['CAPE'][:-1,inds_lat,inds_lon]
CAPE[4416:6624,:,:] = ds_mccaul_2017_d01['CAPE'][:-1,inds_lat,inds_lon]
CAPE[6624:8832,:,:] = ds_mccaul_2018_d01['CAPE'][:-1,inds_lat,inds_lon]
CAPE[8832:11040,:,:] = ds_mccaul_2019_d01['CAPE'][:-1,inds_lat,inds_lon]
CAPE[11040:13248,:,:] = ds_mccaul_2020_d01['CAPE'][:-1,inds_lat,inds_lon]
CAPE[np.isnan(CAPE)] = 0
CSI = np.zeros((13248,63,109))
CSI[0:2208,:,:] = ds_mccaul_2015_d01['CSI'][:-1,inds_lat,inds_lon]
CSI[2208:4416,:,:] = ds_mccaul_2016_d01['CSI'][:-1,inds_lat,inds_lon]
CSI[4416:6624,:,:] = ds_mccaul_2017_d01['CSI'][:-1,inds_lat,inds_lon]
CSI[6624:8832,:,:] = ds_mccaul_2018_d01['CSI'][:-1,inds_lat,inds_lon]
CSI[8832:11040,:,:] = ds_mccaul_2019_d01['CSI'][:-1,inds_lat,inds_lon]
CSI[11040:13248,:,:] = ds_mccaul_2020_d01['CSI'][:-1,inds_lat,inds_lon]

for i in range(1, 2208):
    print('loop 1/6 ' + str(i))
    RAINC[i,:,:] = np.subtract(ds_mccaul_2015_d01['RAINC'][i,inds_lat,inds_lon], ds_mccaul_2015_d01['RAINC'][i-1,inds_lat, inds_lon])
    RAINNC[i, :, :] = np.subtract(ds_mccaul_2015_d01['RAINNC'][i, inds_lat, inds_lon],
                                 ds_mccaul_2015_d01['RAINNC'][i - 1, inds_lat, inds_lon])
for i in range(1,2208):
    print('loop 2/6 ' + str(i))
    RAINC[i+2208, :, :] = np.subtract(ds_mccaul_2016_d01['RAINC'][i,inds_lat, inds_lon],
                                      ds_mccaul_2016_d01['RAINC'][i-1,inds_lat, inds_lon])
    RAINNC[i+2208, :, :] = np.subtract(ds_mccaul_2016_d01['RAINNC'][i,inds_lat, inds_lon],
                                      ds_mccaul_2016_d01['RAINNC'][i-1,inds_lat, inds_lon])
for i in range(1,2208):
    print('loop 3/6 ' + str(i))
    RAINC[i+4416,:,:] = np.subtract(ds_mccaul_2017_d01['RAINC'][i,inds_lat, inds_lon],
                                    ds_mccaul_2017_d01['RAINC'][i-1,inds_lat, inds_lon])
    RAINNC[i+4416,:,:] = np.subtract(ds_mccaul_2017_d01['RAINNC'][i,inds_lat, inds_lon],
                                    ds_mccaul_2017_d01['RAINNC'][i-1,inds_lat, inds_lon])
for i in range(1,2208):
    print('loop 4/6 ' + str(i))
    RAINC[i+6624,:,:] = np.subtract(ds_mccaul_2018_d01['RAINC'][i,inds_lat, inds_lon],
                                    ds_mccaul_2018_d01['RAINC'][i-1,inds_lat, inds_lon])
    RAINNC[i+6624,:,:] = np.subtract(ds_mccaul_2018_d01['RAINNC'][i,inds_lat, inds_lon],
                                    ds_mccaul_2018_d01['RAINNC'][i-1,inds_lat, inds_lon])
for i in range(1,2208):
    print('loop 5/6 ' + str(i))
    RAINC[i+8832,:,:] = np.subtract(ds_mccaul_2019_d01['RAINC'][i,inds_lat, inds_lon],
                                    ds_mccaul_2019_d01['RAINC'][i-1,inds_lat, inds_lon])
    RAINNC[i+8832,:,:] = np.subtract(ds_mccaul_2019_d01['RAINNC'][i,inds_lat, inds_lon],
                                    ds_mccaul_2019_d01['RAINNC'][i-1,inds_lat, inds_lon])
for i in range(1,2208):
    print('loop 6/6 ' + str(i))
    RAINC[i+11040,:,:] = np.subtract(ds_mccaul_2020_d01['RAINC'][i,inds_lat, inds_lon],
                                     ds_mccaul_2020_d01['RAINC'][i-1,inds_lat, inds_lon])
    RAINNC[i+11040,:,:] = np.subtract(ds_mccaul_2020_d01['RAINNC'][i,inds_lat, inds_lon],
                                     ds_mccaul_2020_d01['RAINNC'][i-1,inds_lat, inds_lon])

print('CAPExP')
ds['RAINC'][:] = RAINC
ds['RAINNC'][:] = RAINNC
RAIN = np.add(RAINC, RAINNC)
RAIN[:] = np.divide(RAIN[:],3600) # To convert units to km/m²/s

# CAPExP WITH PRECIPITATION THRESHOLDS
K = np.zeros((13248,63,109))
for i in range(0,13248):
    print('loop ' + str(i+1) + ' of 13248')
    for j in range(0,63):
        for k in range(0,109):
            if RAIN[i,j,k] > 20:
                K[i,j-1:j+2,k-1:k+2] = 1
            elif RAIN[i,j,k] >= 2*np.average(RAIN[i,j-2:j+3,k-2:k+3]):
                K[i,j-1:j+2,k-1:k+2] = 1

PREC = np.multiply(RAIN, K)
ds['CAPExP_R'][:] = np.multiply(CAPE[:], PREC[:])
ds['CAPExP_R'][:] = np.multiply(ds['CAPExP_R'][:], 0.000000000013) # constant, see Romps et al. 2014 (eta/E)
ds['CAPExP_R'][:] = np.multiply(ds['CAPExP_R'][:], 3600000000) # To convert from (m2s)-1 to (km2hour)-1

# CAPExP WITH CSI
K = np.zeros((13248,63,109))
K = np.where(CSI == 3, 1, 0)
PREC = np.multiply(RAIN, K)
ds['CAPExP_CSI'][:] = np.multiply(np.multiply(np.multiply(PREC[:], CAPE[:]),0.000000000013),3600000000)


ds.close()
