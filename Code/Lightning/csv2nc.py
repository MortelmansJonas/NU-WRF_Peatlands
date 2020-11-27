import pandas as pd
import xarray

data = pd.read_csv('/scratch/leuven/336/vsc33651/nu-wrf-dev/Great_Slave_Lake/SL_Short.csv', delimiter = ',')
data['Date'] = pd.to_datetime(data['Date'])
data['Time'] = pd.to_datetime(data['Time'])
data['Type'] = data['Type'].astype('category')
data['Latitude'] = round(data['Latitude'],2)
data['Longitude'] = round(data['Longitude'],2)
data['Count'] = 1

flashes= data.groupby(['Date', 'Latitude', 'Longitude', 'Type']).agg({'Count': ['sum'],'EPC':['mean']}).reset_index()
CCPL = flashes[flashes.Type == 'C'].reset_index()
CGPL = flashes[flashes.Type == 'G'].reset_index()
del CCPL['index']
del CGPL['index']
CCPL.columns = ['Date', 'Latitude', 'Longitude', 'Type', 'Count', 'avg_EPC']
CGPL.columns = ['Date', 'Latitude', 'Longitude', 'Type', 'Count', 'avg_EPC']

# Cloud-to-Cloud
xr1 = xarray.Dataset.from_dataframe(CCPL)
print(xr1)
xr1['Count'].attrs={'units':'flashes/day', 'long_name':'Flashdensity'}
xr1['avg_EPC'].attrs={'units':'kA', 'long_name':'average estimated peak current'}

xr1.to_netcdf('/scratch/leuven/336/vsc33651/nu-wrf-dev/Great_Slave_Lake/vsc2nc_CC.nc')

# Cloud-to-Ground
xr2 = xarray.Dataset.from_dataframe(CGPL)
print(xr2)
xr2['Count'].attrs={'units':'flashes/day', 'long_name':'Flashdensity'}
xr2['avg_EPC'].attrs={'units':'kA', 'long_name':'average estimated peak current'}

xr2.to_netcdf('/scratch/leuven/336/vsc33651/nu-wrf-dev/Great_Slave_Lake/vsc2nc_CG.nc')