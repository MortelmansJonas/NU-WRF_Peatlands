import pygrib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

## Open the GRIB file
file = '/scratch/leuven/336/vsc33651/nu-wrf-dev/Great_Slave_Lake/noahmp36_modis_merra2/postprd/WRFPRS_d02.01'
gr=pygrib.open(file)
for g in gr:
    print(g.typeOfLevel, g.level, g.name)

## Extract the needed variables (latitude, longitude, cloud top and base height.
lat = gr.select(name='Latitude (-90 to +90)', typeOfLevel = 'surface') #To get the index
lon = gr.select(name='East Longitude (0 - 360)', typeOfLevel = 'surface') #To get the index
lats = np.array(gr[497].values) #To get the actual values
lons = np.array(gr[498].values) #To get the actual values
lons = np.subtract(lons, 360) #To go from 0 - 360 to -180 - 180
ctophgt = gr.select(name='Geopotential Height', typeOfLevel = 'cloudTop') #To get the index
ctoptemp = gr.select(name='Temperature', typeOfLevel = 'cloudTop') #Just in case it's needed
cbothgt = gr.select(name='Geopotential Height', typeOfLevel = 'cloudBase') #To get the index
ctop = np.array(gr[465].values) #To get the actual values
cbot = np.array(gr[463].values) #To get the actual values
cdim = np.subtract(ctop, cbot) #To get the cloud vertical dimension
cdim_km = np.divide(cdim, 1000) #To get the dimensions in km instead of m

## Calculate the lightning flashes per minute according to the Price and Rind scheme (F = 3.44*10⁽-5)*cdim_km(4.9)
power = np.power(cdim_km, 4.9)
F_per_min = np.multiply(power, 0.0000344)
print(F_per_min)
F_per_hour = np.multiply(F_per_min, 60) #To get hourly flashes (multiply by 60)
print(F_per_hour)

## Put the lightning flashes in a pandas Dataframe with latitude and longitue
print('dataframe')
F_PH = F_per_hour.reshape(78408) #To get all records in one long array
lats = lats.reshape(78408) #To get all records in one long array
lons = lons.reshape(78408) #To get all records in one long array
Flashes = np.zeros((78408,3))
Flashes[:,2] = F_PH
Flashes[:,0] = lats
Flashes[:,1] = lons
F_DF = pd.DataFrame(Flashes, columns=['Latitude', 'Longitude', 'Flashes_per_hour'])
print(F_DF)
print(np.amax(F_per_hour))