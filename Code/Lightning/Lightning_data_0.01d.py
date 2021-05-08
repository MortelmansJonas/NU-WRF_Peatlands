# ---------------------------------------------------------------------------------------------
# MODULES
# ---------------------------------------------------------------------------------------------
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()
from netCDF4 import Dataset
import datetime as dt

# ---------------------------------------------------------------------------------------------
# LOAD DATA
# ---------------------------------------------------------------------------------------------
data = pd.read_csv('/data/leuven/336/vsc33651/projects/nu-wrf-dev/Lightning_data/Slave_lake.txt', delimiter = '\s+')

print(data.head(5))
print(data.columns)
print(data.dtypes)

data['Date'] = pd.to_datetime(data['Date'])
time = pd.to_datetime(data['Time'], format="%H:%M:%S.%f")
hour = time.dt.hour
data['Date'] = data['Date'] + pd.to_timedelta(hour, unit='h')
print(data.head(35))
data['Date'] = data['Date'].astype(int)
data['Type'] = data['Type'].astype('category')
data['Latitude'] = round(data['Latitude'],2) # 1 for d01, 2 for d02
data['Longitude'] = round(data['Longitude'],2)
data['Count'] = 1

print(data.dtypes)
print(data.head(35))

# ---------------------------------------------------------------------------------------------
# GET HOURLY FLASH DENSITY
# ---------------------------------------------------------------------------------------------
flashes_per_day= data.groupby(['Date', 'Type']).agg({'Count': ['sum'],'EPC':['mean']}).reset_index()
CC = flashes_per_day[flashes_per_day.Type == 'C'].reset_index()
CG = flashes_per_day[flashes_per_day.Type == 'G'].reset_index()

# ---------------------------------------------------------------------------------------------
# PLOT IT TO GET AN OVERVIEW
# ---------------------------------------------------------------------------------------------
plt.plot(CC['Date'], CC['Count'],'.', label='CC')
plt.plot(CG['Date'], CG['Count'],'.',  label='CG')
plt.xlabel("date")
plt.yscale("log")
plt.grid()
plt.ylabel("number of flashes")
plt.title('Number of flashes per day per location')
plt.legend()
plt.show()

# ---------------------------------------------------------------------------------------------
# GET HOURLY FLASH DENSITY AND LOCATION
# ---------------------------------------------------------------------------------------------
flashes= data.groupby(['Date', 'Latitude', 'Longitude', 'Type']).agg({'Count': ['sum'],'EPC':['mean']}).reset_index()
#print(flashes)
CCPL = flashes[flashes.Type == 'C'].reset_index()
CGPL = flashes[flashes.Type == 'G'].reset_index()
del CCPL['index']
del CGPL['index']
CCPL.columns = ['Date', 'Latitude', 'Longitude', 'Type', 'Count', 'avg_EPC']
CGPL.columns = ['Date', 'Latitude', 'Longitude', 'Type', 'Count', 'avg_EPC']
print(CCPL)
print(CGPL)

# ---------------------------------------------------------------------------------------------
# PLOT IT TO GET AN OVERVIEW
# ---------------------------------------------------------------------------------------------
plt.plot(CCPL.Date, CCPL['Count'], '.', label='CC')
plt.plot(CGPL.Date, CGPL['Count'], '.', label='CG')
plt.xlabel("date")
plt.ylabel("number of flashes")
plt.title('Number of flashes per day per location')
plt.legend()
plt.show()

# ---------------------------------------------------------------------------------------------
# CREATE NETCDF FILE
# ---------------------------------------------------------------------------------------------
nf= '/scratch/leuven/336/vsc33651/nu-wrf-dev/wrfout_nc_files/Slave_lake_daily_0.01d.nc'
ds = Dataset(nf, 'w', format='NETCDF4')

# Create dimensions
time = ds.createDimension('time',None)
lat = ds.createDimension('lat', None)
lon = ds.createDimension('lon', None)

# Add NetCDF Variables
times = ds.createVariable('time', 'f4', ('time',))
lats = ds.createVariable('lat', 'f4', ('lat',))
lons = ds.createVariable('lon', 'f4', ('lon',))
Flashdensity_CC = ds.createVariable('Flashdensity_CC', int, ( 'time', 'lat', 'lon',), fill_value=0)
Flashdensity_CC.units = 'Amount/hour'
Flashdensity_CG = ds.createVariable('Flashdensity_CG', int, ( 'time', 'lat', 'lon',), fill_value=0)
Flashdensity_CG.units = 'Amount/hour'

Flashdensity_CC[:,:,:] = 0
Flashdensity_CG[:,:,:] = 0

# Create arrays of the necessary variables
ta = np.unique(flashes['Date']) # To get each day/hour only once
Time_array = np.array(ta)
times[:] = Time_array # Assign time values
print(Time_array)
# Assign Latitude and Longitude Values
lats[:] = np.arange(np.nanmin(flashes['Latitude']),np.nanmax(flashes['Latitude'])+0.01,0.01)
lons[:] = np.arange(np.nanmin(flashes['Longitude']),np.nanmax(flashes['Longitude'])+0.01,0.01)

# Assign NETCDF Data Values
for dim in ds.dimensions.values():
    print(dim)
for var in ds.variables.values():
    print(var)

# ---------------------------------------------------------------------------------------------
# FILL IN CLOUD-TO-CLOUD DATA (CC)
# ---------------------------------------------------------------------------------------------
print('Flashdensity_CC and EPC_CC')
for i in range(0,len(CCPL)):
    print('loop 1 '+str(i) + '/' + str(len(CGPL)))
    times_indices = np.where(Time_array == CCPL['Date'][i])[0]
    lats_indices = np.where(lats == CCPL['Latitude'][i])[0]
    lons_indices = np.where(lons == CCPL['Longitude'][i])[0]
    Stack_CC = np.hstack((times_indices, lats_indices, lons_indices))
    j = Stack_CC[0]
    k = Stack_CC[1]
    l = Stack_CC[2]
    Flashdensity_CC[j,k,l] = CCPL['Count'][i]

# Just to check not everythin is still 0
print(np.amax(Flashdensity_CC))
print(np.amax(EPC_CC))

# ---------------------------------------------------------------------------------------------
# FILL IN CLOUD-TO-GROUND DATA (CG)
# ---------------------------------------------------------------------------------------------
print('Flashdensity_CG and EPC_CG')

for i in range(0,len(CGPL)):
    print('loop 2 '+str(i) + '/' + str(len(CGPL)))
    times_indices = np.where(Time_array == CGPL['Date'][i])[0]
    lats_indices = np.where(lats == CGPL['Latitude'][i])[0]
    lons_indices = np.where(lons == CGPL['Longitude'][i])[0]
    Stack_CG = np.hstack((times_indices, lats_indices, lons_indices))
    j = Stack_CG[0]
    k = Stack_CG[1]
    l = Stack_CG[2]
    Flashdensity_CG[j,k,l] = CGPL['Count'][i]

# Just to check not everything is still 0
print(np.amax(Flashdensity_CG))
print(np.amax(EPC_CG))
# For ease: convert times back to easy readable format
TA = np.unique(pd.to_datetime(flashes['Date']))
Time_array = np.array(TA)
times[:] = Time_array  # Assign time values
ds.close()
