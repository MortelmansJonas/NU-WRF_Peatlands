# ---------------------------------------------------------------------------------------------
# MODULES
# ---------------------------------------------------------------------------------------------
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from sklearn.linear_model import LinearRegression as lr
from sklearn.preprocessing import PolynomialFeatures as pf
# ---------------------------------------------------------------------------------------------
# LOAD DATA
# ---------------------------------------------------------------------------------------------
ds_d01 = Dataset('/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/wrfout_nc_files/d01_all.nc','r')
ds_d02 = Dataset('/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/wrfout_nc_files/domain2_at_domain1_all.nc','r')
ds_obs = Dataset('/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/wrfout_nc_files/CLDN_at_domain1_all_v4.nc')

lat = ds_d01['lat'][:]
lon = ds_d01['lon'][:]

L = ds_obs['Flashdensity_CC'][:,:,:].data.flatten() + ds_obs['Flashdensity_CG'][:,:,:].data.flatten()
print(L)
# ---------------------------------------------------------------------------------------------
# NORMALIZE DATA
# ---------------------------------------------------------------------------------------------
Ln = L/np.nanmax(L)
LPI_n_d01 = ds_d01['LPI'][:]/np.nanmax(ds_d01['LPI'][:])
LPI_n_d02 = ds_d02['LPI'][:]/np.nanmax(ds_d02['LPI'][:])
LTG3_n_d01 = ds_d01['LTG3'][:]/np.nanmax(ds_d01['LTG3'][:])
LTG3_n_d02 = ds_d02['LTG3'][:]/np.nanmax(ds_d02['LTG3'][:])
PR92W_n_d01 = ds_d01['PR92W'][:]/np.nanmax(ds_d01['PR92W'][:])
PR92W_n_d02 = ds_d02['PR92W'][:]/np.nanmax(ds_d02['PR92W'][:])
CAPExP_R_n_d01 = ds_d01['CAPExP_R'][:]/np.nanmax(ds_d01['CAPExP_R'][:])
CAPExP_R_n_d02 = ds_d02['CAPExP_R'][:]/np.nanmax(ds_d02['CAPExP_R'][:])
CAPExP_CSI_n_d01 = ds_d01['CAPExP_CSI'][:]/np.nanmax(ds_d01['CAPExP_CSI'][:])
CAPExP_CSI_n_d02 = ds_d02['CAPExP_CSI'][:]/np.nanmax(ds_d02['CAPExP_CSI'][:])

# ---------------------------------------------------------------------------------------------
# bins
# ---------------------------------------------------------------------------------------------
# specify bin edges
my_bins = np.linspace(-3,0,20)
# DOMAIN 1
n_d01, bins_log10_d01, patches_d01 = plt.hist(np.log10(Ln[Ln!=0]),bins = my_bins)
plt.close()
n_CAPExP_R_d01, bins_CAPExP_R_log10_d01, patches_CAPExP_R_d01 = plt.hist(np.log10(CAPExP_R_n_d01[CAPExP_R_n_d01!=0]),bins=my_bins)
plt.close()
n_CAPExP_CSI_d01, bins_CAPExP_CSI_log10_d01, patches_CAPExP_CSI_d01 = plt.hist(np.log10(CAPExP_CSI_n_d01[CAPExP_CSI_n_d01!=0]),bins=my_bins)
plt.close()
n_lpi_d01, bins_lpi_log10_d01, patches_lpi_d01 = plt.hist(np.log10(LPI_n_d01[LPI_n_d01!=0]),bins=my_bins)
plt.close()
n_LTG3_d01, bins_LTG3_log10_d01, patches_LTG3_d01 = plt.hist(np.log10(LTG3_n_d01[LTG3_n_d01!=0]),bins=my_bins)
plt.close()
n_PR92W_d01, bins_PR92W_log10_d01, patches_PR92W_d01 = plt.hist(np.log10(PR92W_n_d01[PR92W_n_d01!=0]),bins=my_bins)
plt.close()

# Set every bin that doesn't contain data to nan
for i in range(len(n_d01)):
    if n_d01[i]==0:
        n_d01[i] = np.nan
    else:
        break
for i in range(len(n_d01)):
    if n_CAPExP_R_d01[i] == 0:
        n_CAPExP_R_d01[i] = np.nan
    else:
        break
for i in range(len(n_d01)):
    if n_CAPExP_CSI_d01[i] == 0:
        n_CAPExP_CSI_d01[i] = np.nan
    else:
        break
for i in range(len(n_d01)):
    if n_lpi_d01[i] == 0:
        n_lpi_d01[i] = np.nan
    else:
        break
for i in range(len(n_d01)):
    if n_LTG3_d01[i] == 0:
        n_LTG3_d01[i] = np.nan
    else:
        break
for i in range(len(n_d01)):
    if n_PR92W_d01[i] == 0:
        n_PR92W_d01[i] = np.nan
    else:
        break

#plots
plt.plot(n_PR92W_d01, n_d01, '.')
plt.xlabel('PR92W')
plt.ylabel('Obs')
plt.title('d01')
plt.show()

plt.plot(n_lpi_d01, n_d01, '.')
plt.xlabel('LPI')
plt.ylabel('Obs')
plt.title('d01')
plt.show()

plt.plot(n_LTG3_d01, n_d01, '.')
plt.xlabel('LTG3')
plt.ylabel('Obs')
plt.title('d01')
plt.show()

plt.plot(n_CAPExP_R_d01, n_d01, '.')
plt.xlabel('CAPExP_R')
plt.ylabel('Obs')
plt.title('d01')
plt.show()

plt.plot(n_CAPExP_CSI_d01, n_d01, '.')
plt.xlabel('CAPExP_CSI')
plt.ylabel('Obs')
plt.title('d01')
plt.show()

# DOMAIN 2
n_d02, bins_log10_d02, patches_d02 = plt.hist(np.log10(Ln[Ln!=0]),bins = my_bins)
plt.close()
n_CAPExP_R_d02, bins_CAPExP_R_log10_d02, patches_CAPExP_R_d02 = plt.hist(np.log10(CAPExP_R_n_d02[CAPExP_R_n_d02!=0]),bins=my_bins)
plt.close()
n_CAPExP_CSI_d02, bins_CAPExP_CSI_log10_d02, patches_CAPExP_CSI_d02 = plt.hist(np.log10(CAPExP_CSI_n_d02[CAPExP_CSI_n_d02!=0]),bins=my_bins)
plt.close()
n_lpi_d02, bins_lpi_log10_d02, patches_lpi_d02 = plt.hist(np.log10(LPI_n_d02[LPI_n_d02!=0]),bins=my_bins)
plt.close()
n_LTG3_d02, bins_LTG3_log10_d02, patches_LTG3_d02 = plt.hist(np.log10(LTG3_n_d02[LTG3_n_d02!=0]),bins=my_bins)
plt.close()
n_PR92W_d02, bins_PR92W_log10_d02, patches_PR92W_d02 = plt.hist(np.log10(PR92W_n_d02[PR92W_n_d02!=0]),bins=my_bins)
plt.close()

# Set every bin that doesn't contain data to nan
for i in range(len(n_d02)):
    if n_d02[i]==0:
        n_d02[i] = np.nan
    else:
        break
for i in range(len(n_d02)):
    if n_CAPExP_R_d02[i] == 0:
        n_CAPExP_R_d02[i] = np.nan
    else:
        break
for i in range(len(n_d02)):
    if n_CAPExP_CSI_d02[i] == 0:
        n_CAPExP_CSI_d02[i] = np.nan
    else:
        break
for i in range(len(n_d02)):
    if n_lpi_d02[i] == 0:
        n_lpi_d02[i] = np.nan
    else:
        break
for i in range(len(n_d02)):
    if n_LTG3_d02[i] == 0:
        n_LTG3_d02[i] = np.nan
    else:
        break
for i in range(len(n_d02)):
    if n_PR92W_d02[i] == 0:
        n_PR92W_d02[i] = np.nan
    else:
        break

#plots
plt.plot(n_PR92W_d02, n_d02, '.')
plt.xlabel('PR92W')
plt.ylabel('Obs')
plt.title('d02')
plt.show()

plt.plot(n_lpi_d02, n_d02, '.')
plt.xlabel('LPI')
plt.ylabel('Obs')
plt.title('d02')
plt.show()

plt.plot(n_LTG3_d02, n_d02, '.')
plt.xlabel('LTG3')
plt.ylabel('Obs')
plt.title('d02')
plt.show()

plt.plot(n_CAPExP_R_d02, n_d02, '.')
plt.xlabel('CAPExP_R')
plt.ylabel('Obs')
plt.title('d02')
plt.show()

plt.plot(n_CAPExP_CSI_d02, n_d02, '.')
plt.xlabel('CAPExP_CSI')
plt.ylabel('Obs')
plt.title('d02')
plt.show()


# ---------------------------------------------------------------------------------------------
# PERFORM LINEAR REGRESSION
# ---------------------------------------------------------------------------------------------
# DOMAIN 1
reg_PR92W_d01 = lr().fit(np.log10(n_PR92W_d01.reshape((-1,1))), np.log10(n_d01))
print('PR92W_d01')
print(reg_PR92W_d01.coef_)
print(reg_PR92W_d01.intercept_)
print(reg_PR92W_d01.score(np.log10(n_PR92W_d01.reshape((-1,1))), np.log10(n_d01)))
# linear_PR92W_d01 = np.add(reg_PR92W_d01.intercept_,np.multiply(reg_PR92W_d01.coef_,ds_d01['PR92W'][:]))
linear_PR92W_d01 = np.multiply(reg_PR92W_d01.coef_,ds_d01['PR92W'][:])
print(np.min(linear_PR92W_d01))
print(np.max(linear_PR92W_d01))

reg_LPI_d01 = lr().fit(np.log10(n_lpi_d01.reshape((-1,1))), np.log10(n_d01))
print('LPI_d01')
print(reg_LPI_d01.coef_)
print(reg_LPI_d01.intercept_)
print(reg_LPI_d01.score(np.log10(n_lpi_d01.reshape((-1,1))), np.log10(n_d01)))
# linear_LPI_d01 = np.add(reg_LPI_d01.intercept_,np.multiply(reg_LPI_d01.coef_,ds_d01['LPI'][:]))
linear_LPI_d01 =np.multiply(reg_LPI_d01.coef_,ds_d01['LPI'][:])
print(np.min(linear_LPI_d01))
print(np.max(linear_LPI_d01))

n_LTG3_d01[n_LTG3_d01==0] = np.nan # needed because log10 of 0 gives NaN
# n_LTG3_d01[np.isnan(n_LTG3_d01)] = 0

reg_LTG3_d01 = lr().fit(np.log10(n_LTG3_d01[~np.isnan(n_LTG3_d01)].reshape((-1,1))), np.log10(n_d01[~np.isnan(n_LTG3_d01)]))
print('LTG3_d01')
print(reg_LTG3_d01.coef_)
print(reg_LTG3_d01.intercept_)
print(reg_LTG3_d01.score(np.log10(n_LTG3_d01[~np.isnan(n_LTG3_d01)].reshape((-1,1))), np.log10(n_d01[~np.isnan(n_LTG3_d01)])))
# linear_LTG3_d01 = np.add(reg_LTG3_d01.intercept_,np.multiply(reg_LTG3_d01.coef_,ds_d01['LTG3'][:]))
linear_LTG3_d01 = np.multiply(reg_LTG3_d01.coef_,ds_d01['LTG3'][:])
print(np.min(linear_LTG3_d01))
print(np.max(linear_LTG3_d01))

reg_CAPExP_CSI_d01 = lr().fit(np.log10(n_CAPExP_CSI_d01.reshape((-1,1))), np.log10(n_d01))
print('CAPExP_CSI_d01')
print(reg_CAPExP_CSI_d01.coef_)
print(reg_CAPExP_CSI_d01.intercept_)
print(reg_CAPExP_CSI_d01.score(np.log10(n_CAPExP_CSI_d01.reshape((-1,1))), np.log10(n_d01)))
# linear_CAPExP_CSI_d01 = np.add(reg_CAPExP_CSI_d01.intercept_,np.multiply(reg_CAPExP_CSI_d01.coef_,ds_d01['CAPExP_CSI'][:]))
linear_CAPExP_CSI_d01 = np.multiply(reg_CAPExP_CSI_d01.coef_,ds_d01['CAPExP_CSI'][:])
print(np.min(linear_CAPExP_CSI_d01))
print(np.max(linear_CAPExP_CSI_d01))

reg_CAPExP_R_d01 = lr().fit(np.log10(n_CAPExP_R_d01.reshape((-1,1))), np.log10(n_d01))
print('CAPExP_R_d01')
print(reg_CAPExP_R_d01.coef_)
print(reg_CAPExP_R_d01.intercept_)
print(reg_CAPExP_R_d01.score(np.log10(n_CAPExP_R_d01.reshape((-1,1))), np.log10(n_d01)))
# linear_CAPExP_R_d01 = np.add(reg_CAPExP_R_d01.intercept_,np.multiply(reg_CAPExP_R_d01.coef_,ds_d01['CAPExP_R'][:]))
linear_CAPExP_R_d01 = np.multiply(reg_CAPExP_R_d01.coef_,ds_d01['CAPExP_R'][:])
print(np.min(linear_CAPExP_R_d01))
print(np.max(linear_CAPExP_R_d01))

# DOMAIN 2

reg_PR92W_d02 = lr().fit(np.log10(n_PR92W_d02.reshape((-1,1))), np.log10(n_d02))
print('PR92W_d02')
print(reg_PR92W_d02.coef_)
print(reg_PR92W_d02.intercept_)
print(reg_PR92W_d02.score(np.log10(n_PR92W_d02.reshape((-1,1))), np.log10(n_d02)))
# linear_PR92W_d02= np.add(reg_PR92W_d02.intercept_,np.multiply(reg_PR92W_d02.coef_,ds_d02['PR92W'][:]))
linear_PR92W_d02= np.multiply(reg_PR92W_d02.coef_,ds_d02['PR92W'][:])
print(np.min(linear_PR92W_d02))
print(np.max(linear_PR92W_d02))

reg_LPI_d02 = lr().fit(np.log10(n_lpi_d02.reshape((-1,1))), np.log10(n_d02))
print('LPI_d02')
print(reg_LPI_d02.coef_)
print(reg_LPI_d02.intercept_)
print(reg_LPI_d02.score(np.log10(n_lpi_d02.reshape((-1,1))), np.log10(n_d02)))
# linear_LPI_d02 = np.add(reg_LPI_d02.intercept_,np.multiply(reg_LPI_d02.coef_,ds_d02['LPI'][:]))
linear_LPI_d02 = np.multiply(reg_LPI_d02.coef_,ds_d02['LPI'][:])
print(np.min(linear_LPI_d02))
print(np.max(linear_LPI_d02))

n_LTG3_d02[n_LTG3_d02==0] = np.nan # needed because log10 of 0 gives NaN
# n_LTG3_d02[np.isnan(n_LTG3_d02)] = 0.1

reg_LTG3_d02 = lr().fit(np.log10(n_LTG3_d02[~np.isnan(n_LTG3_d02)].reshape((-1,1))), np.log10(n_d02[~np.isnan(n_LTG3_d02)]))
print('LTG3_d02')
print(reg_LTG3_d02.coef_)
print(reg_LTG3_d02.intercept_)
print(reg_LTG3_d02.score(np.log10(n_LTG3_d02[~np.isnan(n_LTG3_d02)].reshape((-1,1))), np.log10(n_d02[~np.isnan(n_LTG3_d02)])))
# linear_LTG3_d02 = np.add(reg_LTG3_d02.intercept_,np.multiply(reg_LTG3_d02.coef_,ds_d02['LTG3'][:]))
linear_LTG3_d02 = np.multiply(reg_LTG3_d02.coef_,ds_d02['LTG3'][:])
print(np.min(linear_LTG3_d02))
print(np.max(linear_LTG3_d02))

reg_CAPExP_CSI_d02 = lr().fit(np.log10(n_CAPExP_CSI_d02.reshape((-1,1))), np.log10(n_d02))
print('CAPExP_CSI_d02')
print(reg_CAPExP_CSI_d02.coef_)
print(reg_CAPExP_CSI_d02.intercept_)
print(reg_CAPExP_CSI_d02.score(np.log10(n_CAPExP_CSI_d02.reshape((-1,1))), np.log10(n_d02)))
# linear_CAPExP_CSI_d02 = np.add(reg_CAPExP_CSI_d02.intercept_,np.multiply(reg_CAPExP_CSI_d02.coef_,ds_d02['CAPExP_CSI'][:]))
linear_CAPExP_CSI_d02 = np.multiply(reg_CAPExP_CSI_d02.coef_,ds_d02['CAPExP_CSI'][:])
print(np.nanmin(linear_CAPExP_CSI_d02))
print(np.nanmax(linear_CAPExP_CSI_d02))

reg_CAPExP_R_d02 = lr().fit(np.log10(n_CAPExP_R_d02.reshape((-1,1))), np.log10(n_d02))
print('CAPExP_R_d02')
print(reg_CAPExP_R_d02.coef_)
print(reg_CAPExP_R_d02.intercept_)
print(reg_CAPExP_R_d02.score(np.log10(n_CAPExP_R_d02.reshape((-1,1))), np.log10(n_d02)))
# linear_CAPExP_R_d02 = np.add(reg_CAPExP_R_d02.intercept_,np.multiply(reg_CAPExP_R_d02.coef_,ds_d02['CAPExP_R'][:]))
linear_CAPExP_R_d02 = np.multiply(reg_CAPExP_R_d02.coef_,ds_d02['CAPExP_R'][:])
print(np.min(linear_CAPExP_R_d02))
print(np.max(linear_CAPExP_R_d02))


# ---------------------------------------------------------------------------------------------
# PLOTS TO CHECK LINEAR REGRESSION
# ---------------------------------------------------------------------------------------------
plt.scatter(linear_PR92W_d01, L)
plt.ylabel('Observations')
plt.xlabel('Calibrated PR92W')
plt.title('d01')
plt.show()

plt.scatter(linear_LPI_d01, L)
plt.ylabel('Observations')
plt.xlabel('Calibrated LPI')
plt.title('d01')
plt.show()

plt.scatter(linear_LTG3_d01, L)
plt.ylabel('Observations')
plt.xlabel('Calibrated LTG3')
plt.title('d01')
plt.show()

plt.scatter(linear_CAPExP_R_d01, L)
plt.ylabel('Observations')
plt.xlabel('Calibrated CAPExP_R')
plt.title('d01')
plt.show()

plt.scatter(linear_CAPExP_CSI_d01, L)
plt.ylabel('Observations')
plt.xlabel('Calibrated CAPExP_CSI')
plt.title('d01')
plt.show()

plt.scatter(linear_PR92W_d02, L)
plt.ylabel('Observations')
plt.xlabel('Calibrated PR92W')
plt.title('d02')
plt.show()

plt.scatter(linear_LPI_d02, L)
plt.ylabel('Observations')
plt.xlabel('Calibrated LPI')
plt.title('d02')
plt.show()

plt.scatter(linear_LTG3_d02, L)
plt.ylabel('Observations')
plt.xlabel('Calibrated LTG3')
plt.title('d02')
plt.show()

plt.scatter(linear_CAPExP_R_d02, L)
plt.ylabel('Observations')
plt.xlabel('Calibrated CAPExP_R')
plt.title('d02')
plt.show()

plt.scatter(linear_CAPExP_CSI_d02, L)
plt.ylabel('Observations')
plt.xlabel('Calibrated CAPExP_CSI')
plt.title('d02')
plt.show()
# ---------------------------------------------------------------------------------------------
# PUT CALIBRATED DATA IN NEW NETCDF FILE
# ---------------------------------------------------------------------------------------------
# Create .nc file
ds = Dataset('/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/wrfout_nc_files/data_calibrated_ax+b.nc', mode='w', format='NETCDF4')
ds.createDimension('time', 13248)
ds.createDimension('lat', 63)
ds.createDimension('lon', 109)
ds.createVariable('time','int', dimensions=('time',),zlib=True)
ds.variables['time'][:] = ds_d01['time'][:]
ds.createVariable('lat',ds_d01['lat'][:].dtype, dimensions=('lat','lon',),zlib=True)
ds.variables['lat'][:] = ds_d01['lat'][:]
ds.createVariable('lon',ds_d01['lon'][:].dtype, dimensions=('lat','lon',),zlib=True)
ds.variables['lon'][:] = ds_d01['lon'][:]
ds.createVariable('LPI_d01', 'f4', dimensions=('time','lat','lon',), zlib=True)
ds.createVariable('LTG3_d01', 'f4', dimensions=('time','lat','lon',), zlib=True)
ds.createVariable('PR92W_d01', 'f4', dimensions=('time','lat','lon',), zlib=True)
ds.createVariable('CAPExP_R_d01', 'f4', dimensions=('time','lat','lon',), zlib=True)
ds.createVariable('CAPExP_CSI_d01', 'f4', dimensions=('time','lat','lon',), zlib=True)
ds.createVariable('LPI_d02', 'f4', dimensions=('time','lat','lon',), zlib=True)
ds.createVariable('LTG3_d02', 'f4', dimensions=('time','lat','lon',), zlib=True)
ds.createVariable('PR92W_d02', 'f4', dimensions=('time','lat','lon',), zlib=True)
ds.createVariable('CAPExP_R_d02', 'f4', dimensions=('time','lat','lon',), zlib=True)
ds.createVariable('CAPExP_CSI_d02', 'f4', dimensions=('time','lat','lon',), zlib=True)
ds.createVariable('Obs', 'f4', dimensions=('time','lat','lon',), zlib=True)

# Fill it
ds['LPI_d01'][:] = linear_LPI_d01[:]
ds['LPI_d02'][:] = linear_LPI_d02[:]
ds['PR92W_d01'][:] = linear_PR92W_d01[:]
ds['PR92W_d02'][:] = linear_PR92W_d02[:]
ds['LTG3_d01'][:] = linear_LTG3_d01[:]
ds['LTG3_d02'][:] = linear_LTG3_d02[:]
ds['CAPExP_CSI_d01'][:] = linear_CAPExP_CSI_d01[:]
ds['CAPExP_CSI_d02'][:] = linear_CAPExP_CSI_d02[:]
ds['CAPExP_R_d01'][:] = linear_CAPExP_R_d01[:]
ds['CAPExP_R_d02'][:] = linear_CAPExP_R_d02[:]
ds['Obs'][:] = ds_obs['Flashdensity_CC'][:,:,:].data + ds_obs['Flashdensity_CG'][:,:,:].data
ds.close()
