## Import the necessary modules
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()
from netCDF4 import Dataset
from datetime import datetime

## Read in the data
data = pd.read_csv('/scratch/leuven/336/vsc33651/nu-wrf-dev/Great_Slave_Lake/SL_Short.csv', delimiter = ',')

print(data.head(5))
print(data.columns)
print(data.dtypes)

data['Date'] = pd.to_datetime(data['Date']).astype(int) # Set as integer to later on be able to compare to the CCPL data. If it's not done, there are some issues
data['Time'] = pd.to_datetime(data['Time'])
data['Type'] = data['Type'].astype('category')
data['Latitude'] = round(data['Latitude'],2)
data['Longitude'] = round(data['Longitude'],2)
data['Count'] = 1
print(data.dtypes)


## Get a Flash density per day
flashes_per_day= data.groupby(['Date', 'Type']).agg({'Count': ['sum'],'EPC':['mean']}).reset_index()
CC = flashes_per_day[flashes_per_day.Type == 'C'].reset_index()
CG = flashes_per_day[flashes_per_day.Type == 'G'].reset_index()

## Plot it
plt.plot(CC['Date'], CC['Count'],'.', label='CC')
plt.plot(CG['Date'], CG['Count'],'.',  label='CG')
plt.xlabel("date")
plt.yscale("log")
plt.grid()
plt.ylabel("number of flashes")
plt.title('Number of flashes per day per location')
plt.legend()
plt.show()

## Get a Flash density per day per latlon
flashes= data.groupby(['Date', 'Latitude', 'Longitude', 'Type']).agg({'Count': ['sum'],'EPC':['mean']}).reset_index()
print(flashes)
CCPL = flashes[flashes.Type == 'C'].reset_index()
CGPL = flashes[flashes.Type == 'G'].reset_index()
del CCPL['index']
del CGPL['index']
CCPL.columns = ['Date', 'Latitude', 'Longitude', 'Type', 'Count', 'avg_EPC']
CGPL.columns = ['Date', 'Latitude', 'Longitude', 'Type', 'Count', 'avg_EPC']
print(CCPL)
print(CGPL)
## Plot it
plt.plot(CCPL.Date, CCPL['Count'], '.', label='CC')
plt.plot(CGPL.Date, CGPL['Count'], '.', label='CG')
plt.xlabel("date")
plt.ylabel("number of flashes")
plt.title('Number of flashes per day per location')
plt.legend()
plt.show()

## Create a netCDF file
nf= '/scratch/leuven/336/vsc33651/nu-wrf-dev/Great_Slave_Lake/Slave_lake.nc'
ds = Dataset(nf, 'w', format='NETCDF4')

# Create dimensions
time = ds.createDimension('time',None)
lat = ds.createDimension('lat', None)
lon = ds.createDimension('lon', None)

# Add NetCDF Variables
times = ds.createVariable('time', 'f4', ('time',))
lats = ds.createVariable('lat', float, ('lat',))
lons = ds.createVariable('lon', float, ('lon',))
Flashdensity_CC = ds.createVariable('Flashdensity_CC', int, ( 'time', 'lat', 'lon',), fill_value=0)
Flashdensity_CC.units = 'Amount/day'
Flashdensity_CG = ds.createVariable('Flashdensity_CG', int, ( 'time', 'lat', 'lon',), fill_value=0)
Flashdensity_CG.units = 'Amount/day'
EPC_CC = ds.createVariable('EPC_CC', float, ( 'time', 'lat', 'lon',), fill_value=0)
EPC_CC.units = 'kA'
EPC_CG = ds.createVariable('EPC_CG', float, ( 'time', 'lat', 'lon',), fill_value=0)
EPC_CG.units = 'kA'

# Create arrays of the necessary variables
ta = np.unique(flashes['Date']) # To get each day only once
Time_array = np.array(ta)
times[:] = Time_array # Assign time values

# Assign Latitude and Longitude Values
lats[:] = np.array(np.unique(flashes['Latitude']))
lons[:] = np.array(np.unique(flashes['Longitude']))

# Assign NETCDF Data Values
for dim in ds.dimensions.values():
    print(dim)
for var in ds.variables.values():
    print(var)

## Flashdensity_CC
print('Flashdensity_CC')
for i in range(0,len(Time_array)): # Loop through all the times
     for j in range(0,len(CCPL)): # Loop through all the CCPL records
         if Time_array[i] == CCPL['Date'][j]: # if the time in both arrays is the same, go to the following loop
             for k in range(0,len(lats)): # Loop through all latitudes
                 if lats[k] == CCPL['Latitude'][j]: # If the latitudes match, go to the next loop
                     for m in range(0,len(lons)): # Loop through all longitudes
                         if lons[m] == CCPL['Longitude'][j]: # If the longitudes match, put the CGPL record in Flashdensity_CG
                             Flashdensity_CC[i,k,m] = CCPL['Count'][j]

## Flashdensity_CG
print('Flashdensity_CG')
for i in range(0,len(Time_array)): # Loop through all the times
     for j in range(0,len(CGPL)): # Loop through all the CGPL records
         if Time_array[i] == CGPL['Date'][j]: # if the time in both arrays is the same, go to the following loop
             for k in range(0,len(lats)): # Loop through all latitudes
                 if lats[k] == CGPL['Latitude'][j]: # If the latitudes match, go to the next loop
                     for m in range(0,len(lons)): # Loop through all longitudes
                         if lons[m] == CGPL['Longitude'][j]: # If the longitudes match, put the CCPL record in EPC_CC
                             Flashdensity_CG[i,k,m] = CGPL['Count'][j]

## EPC_CC
print('EPC_CC')
for i in range(0,len(Time_array)): # Loop through all the times
     for j in range(0,len(CCPL)): # Loop through all the CCPL records
         if Time_array[i] == CCPL['Date'][j]: # if the time in both arrays is the same, go to the following loop
             for k in range(0,len(lats)): # Loop through all latitudes
                 if lats[k] == CCPL['Latitude'][j]: # If the latitudes match, go to the next loop
                     for m in range(0,len(lons)): # Loop through all longitudes
                         if lons[m] == CCPL['Longitude'][j]: # If the longitudes match, put the CCPL record in EPC_CG
                             EPC_CC[i,k,m] = CCPL['avg_EPC'][j]

## EPC_CG
print('EPC_CG')
for i in range(0,len(Time_array)): # Loop through all the times
     for j in range(0,len(CGPL)): # Loop through all the CGPL records
         if Time_array[i] == CGPL['Date'][j]: # if the time in both arrays is the same, go to the following loop
             for k in range(0,len(lats)): # Loop through all latitudes
                 if lats[k] == CGPL['Latitude'][j]: # If the latitudes match, go to the next loop
                     for m in range(0,len(lons)): # Loop through all longitudes
                         if lons[m] == CGPL['Longitude'][j]: # If the longitudes match, put the CCPL record in Flashdensity_CC
                             EPC_CG[i,k,m] = CGPL['avg_EPC'][j]

# For ease: convert times back to easy readable format
TA = np.unique(pd.to_datetime(flashes['Date']))
Time_array = np.array(TA)
times[:] = Time_array # Assign time values
