import pygrib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from netCDF4 import Dataset

date_from = '2015-06-01'
date_to = '2015-06-02'

path = '/scratch/leuven/336/vsc33651/nu-wrf-dev/Great_Slave_Lake/noahmp36_modis_merra2/postprd/'
filename_root = 'WRFPRS_'
domain = 'd02'
output_filename = os.path.join('/scratch/leuven/317/vsc31786/',filename_root+domain+'.nc')

# determine length of time vector
ndays = (pd.to_datetime(date_to) - pd.to_datetime(date_from)).days
nhours = ndays * 24
hours = np.linspace(1,nhours,nhours)
# with leading zeros
length_counter = len(str(nhours))
cstr = "%0"+str(length_counter)+"d"
timestep = cstr % 1
filename = os.path.join(path,filename_root+domain+'.'+timestep)
gr=pygrib.open(filename)
lat = np.array(gr[497].values) #To get the actual values
lon = np.array(gr[498].values) #To get the actual values
lon = np.subtract(lon, 360) #To go from 0 - 360 to -180 - 180

ds = Dataset(output_filename, mode='w', format='NETCDF4')
timeunit = 'hours since 2000-01-01 00:00'
ds.createDimension('time', nhours)
ds.createDimension('lat', np.shape(lat)[0])
ds.createDimension('lon', np.shape(lon)[1])
ds.createVariable('time',hours.dtype, dimensions=('time',),zlib=True)
ds.variables['time'][:] = hours + (pd.to_datetime(date_from) - pd.to_datetime('2000-01-01')).days*24 - 1
ds.createVariable('lat',lat.dtype, dimensions=('lat','lon'),zlib=True)
ds.variables['lat'][:] = lat
ds.createVariable('lon',lon.dtype, dimensions=('lat','lon'),zlib=True)
ds.variables['lon'][:] = lon
ds.createVariable('F_PH','int', dimensions=('time','lat','lon'),zlib=True)

ds.variables['time'].setncatts({'long_name': 'time', 'units': timeunit})
ds.variables['lat'].setncatts({'long_name': 'latitude',  'units': 'degrees_north'})
ds.variables['lon'].setncatts({'long_name': 'longitude', 'units': 'degrees_east'})
ds.variables['F_PH'].setncatts({'long_name': 'flashes_per_hour',  'units': 'counts'})

for h in range(1,nhours+1):
    ## Open the GRIB file
    timestep = cstr % h
    filename = os.path.join(path,filename_root+domain+'.'+timestep)
    print('processing '+filename)
    gr=pygrib.open(filename)
    #for g in gr:
    #    print(g.typeOfLevel, g.level, g.name)

    ## Extract the needed variables (latitude, longitude, cloud top and base height.
    #lat = gr.select(name='Latitude (-90 to +90)', typeOfLevel = 'surface') #To get the index
    #lon = gr.select(name='East Longitude (0 - 360)', typeOfLevel = 'surface') #To get the index
    #ctophgt = gr.select(name='Geopotential Height', typeOfLevel = 'cloudTop') #To get the index
    #ctoptemp = gr.select(name='Temperature', typeOfLevel = 'cloudTop') #Just in case it's needed
    #cbothgt = gr.select(name='Geopotential Height', typeOfLevel = 'cloudBase') #To get the index
    ctop = np.array(gr[465].values) #To get the actual values
    cbot = np.array(gr[463].values) #To get the actual values
    cdim = np.subtract(ctop, cbot) #To get the cloud vertical dimension
    cdim_km = np.divide(cdim, 1000) #To get the dimensions in km instead of m

    ## Calculate the lightning flashes per minute according to the Price and Rind scheme (F = 3.44*10⁽-5)*cdim_km(4.9)
    # catch negative height
    cdim_km[cdim_km<0.0] = 0.0
    # catch unrealistic high clouds, just in case
    # not sure whether it occurs, for individual hours I saw some pixels popping out with very high lightning, but still in range
    cdim_km[cdim_km>18.3] = 18.3
    power = np.power(cdim_km, 4.9)
    F_per_min = np.multiply(power, 0.0000344)
    #print(F_per_min)
    F_PH = np.multiply(F_per_min, 60) #To get hourly flashes (multiply by 60)

    ds.variables['F_PH'][h-1,:,:] = F_PH

ds.close()
