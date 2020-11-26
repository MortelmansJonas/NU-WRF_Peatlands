## Import the necessary modules
import pandas as pd
import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()
import netCDF4 as nc
from netCDF4 import Dataset

## Read in the data
data = pd.read_csv('/scratch/leuven/336/vsc33651/nu-wrf-dev/Great_Slave_Lake/SL_Short.csv', delimiter = ',')

print(data.head(5))
print(data.columns)
print(data.dtypes)

data['Date'] = pd.to_datetime(data['Date'])
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
print(CCPL)
CGPL = flashes[flashes.Type == 'G'].reset_index()
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
# Method 1
print('netCDF file: Method 1')
df = data.to_xarray()
print(df)
df.to_netcdf('/scratch/leuven/336/vsc33651/nu-wrf-dev/Great_Slave_Lake/Slave_lake2.nc')

# Doesn't give flashes/type -> can be split? (Doesn't work yet
print('netCDF CC')
CCPL.drop('index', axis=1)
dCC = CCPL.to_xarray()
dCC.to_netcdf('/scratch/leuven/336/vsc33651/nu-wrf-dev/Great_Slave_Lake/Slave_lake_CC.nc')

print('netCDF CG')
CGPL.drop('index', axis=1)
dCG = CGPL.to_xarray()
dCG.to_netcdf('/scratch/leuven/336/vsc33651/nu-wrf-dev/Great_Slave_Lake/Slave_lake_CG.nc')

# Method 2 (Gives netCDF file with right variables, but empty)
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
Flashdensity_CC = ds.createVariable('Flashdensity_CC', int, ( 'time', 'lat', 'lon',), fill_value=-9999.0)
Flashdensity_CC.units = 'Amount/day'
Flashdensity_CG = ds.createVariable('Flashdensity_CG', int, ( 'time', 'lat', 'lon',), fill_value=-9999.0)
Flashdensity_CG.units = 'Amount/day'
EPC_CC = ds.createVariable('EPC_CC', float, ( 'time', 'lat', 'lon',), fill_value=-9999.0)
EPC_CC.units = 'kA'
EPC_CG = ds.createVariable('EPC_CG', float, ( 'time', 'lat', 'lon',), fill_value=-9999.0)
EPC_CG.units = 'kA'

# Create arrays of the necessary variables
ta = np.array(CCPL['Date'])
Ta = np.unique(ta) # To get each day only once
Time_array = np.array(Ta)

# Assign Latitude and Longitude Values
lats[:] = np.array(np.unique(CCPL['Latitude']))
lons[:] = np.array(np.unique(CCPL['Longitude']))

# Assign time values
times[:] = Time_array

# Assign NETCDF Data Values
for dim in ds.dimensions.values():
    print(dim)
for var in ds.variables.values():
    print(var)

Flashdensity_CC
Flash_CC = []
for i in range(0,len(times)):
     for j in range(0,len(CCPL)):
         if times[i] == CCPL['Date'][j]:
             for k in range(0,len(lats)):
                 if lats[k] == CCPL['Latitude'][j]:
                     for m in range(0,len(lons)):
                         if lons[m] == CCPL['Longitude'][j]:
                             Flashdensity_CC[i,k,m] = CCPL['Count'][j]
                         Flash_CC.append(Flashdensity_CC)
Flashdensity_CC = Flash_CC

