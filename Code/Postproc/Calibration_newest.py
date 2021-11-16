# ---------------------------------------------------------------------------------------------
# MODULES
# ---------------------------------------------------------------------------------------------
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from sklearn.linear_model import LinearRegression as lr
import seaborn as sns

# ---------------------------------------------------------------------------------------------
# LOAD DATA
# ---------------------------------------------------------------------------------------------
ds_d01 = Dataset('/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/wrfout_nc_files/d01_all.nc','r')
ds_d02 = Dataset('/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/wrfout_nc_files/domain2_at_domain1_all.nc','r')
ds_obs = Dataset('/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/wrfout_nc_files/CLDN_at_domain1_all_v4.nc')

lat = ds_d01['lat'][:]
lon = ds_d01['lon'][:]

L = ds_obs['Flashdensity_CC'][:,:,:].data.flatten() + ds_obs['Flashdensity_CG'][:,:,:].data.flatten()

# ---------------------------------------------------------------------------------------------
# 'The lowest values of the parameterization's output are set to 0'
# ---------------------------------------------------------------------------------------------
LPI_d01 = ds_d01['LPI'][:]
LPI_d01[LPI_d01 == np.nanmin(LPI_d01)] = 0

LTG3_d01 = ds_d01['LTG3'][:]
LTG3_d01[LTG3_d01 == np.nanmin(LTG3_d01)] = 0

PR92W_d01 = ds_d01['PR92W'][:]
PR92W_d01[PR92W_d01 == np.nanmin(PR92W_d01)] = 0

CAPExP_R_d01 = ds_d01['CAPExP_R'][:]
CAPExP_R_d01[CAPExP_R_d01 == np.nanmin(CAPExP_R_d01)] = 0

CAPExP_CSI_d01 = ds_d01['CAPExP_CSI'][:]
CAPExP_CSI_d01[CAPExP_CSI_d01 == np.nanmin(CAPExP_CSI_d01)] = 0

LPI_d02 = ds_d02['LPI'][:]
LPI_d02[LPI_d02 == np.nanmin(LPI_d02)] = 0

LTG3_d02 = ds_d02['LTG3'][:]
LTG3_d02[LTG3_d02 == np.nanmin(LTG3_d02)] = 0

PR92W_d02 = ds_d02['PR92W'][:]
PR92W_d02[PR92W_d02 == np.nanmin(PR92W_d02)] = 0

CAPExP_R_d02 = ds_d02['CAPExP_R'][:]
CAPExP_R_d02[CAPExP_R_d02 == np.nanmin(CAPExP_R_d02)] = 0

CAPExP_CSI_d02 = ds_d02['CAPExP_CSI'][:]
CAPExP_CSI_d02[CAPExP_CSI_d02 == np.nanmin(CAPExP_CSI_d02)] = 0

# -------------------------------------------------------------------------------------------------------------
# 'Second, a simple linear model is built that relates the observed flash rates to the parameterization output'
# -------------------------------------------------------------------------------------------------------------
# First determine cutoff value 'c', done similarly to the sorting in 'Calibration_2.0.py'
Obs = np.sort(L[L!=0])[::-1]

# Sort model and take n highest values (with n = number of observations)
c = Obs.shape[0]

LPI_sorted_d01 = np.sort(LPI_d01.flatten())[::-1]
LPI_c_d01 = LPI_sorted_d01[0:c]
LPI_d01_new = np.where(LPI_d01>=LPI_sorted_d01[c],LPI_d01, np.nan)
print('cutoff value LPI_d01 = ', LPI_sorted_d01[c])

LTG3_sorted_d01 = np.sort(LTG3_d01.flatten())[::-1]
LTG3_c_d01 = LTG3_sorted_d01[0:c]
LTG3_d01_new = np.where(LTG3_d01>=LTG3_sorted_d01[c],LTG3_d01, np.nan)
print('cutoff value LTG3_d01 = ', LTG3_sorted_d01[c])

PR92W_sorted_d01 = np.sort(PR92W_d01.flatten())[::-1]
PR92W_c_d01 = PR92W_sorted_d01[0:c]
PR92W_d01_new = np.where(PR92W_d01>=PR92W_sorted_d01[c],PR92W_d01, np.nan)
print('cutoff value PR92W_d01 = ', PR92W_sorted_d01[c])

CAPExP_CSI_sorted_d01 = np.sort(CAPExP_CSI_d01.flatten())[::-1]
CAPExP_CSI_c_d01 = CAPExP_CSI_sorted_d01[0:c]
CAPExP_CSI_d01_new = np.where(CAPExP_CSI_d01>=CAPExP_CSI_sorted_d01[c],CAPExP_CSI_d01, np.nan)
print('cutoff value CAPExP_CSI_d01 = ', CAPExP_CSI_sorted_d01[c])

CAPExP_R_sorted_d01 = np.sort(CAPExP_R_d01.flatten())[::-1]
CAPExP_R_c_d01 = CAPExP_R_sorted_d01[0:c]
CAPExP_R_d01_new = np.where(CAPExP_R_d01>=CAPExP_R_sorted_d01[c],CAPExP_R_d01, np.nan)
print('cutoff value CAPExP_R_d01 = ', CAPExP_R_sorted_d01[c])

LPI_sorted_d02 = np.sort(LPI_d02.flatten())[::-1]
LPI_c_d02 = LPI_sorted_d02[0:c]
LPI_d02_new = np.where(LPI_d02>=LPI_sorted_d02[c],LPI_d02, np.nan)
print('cutoff value LPI_d02 = ', LPI_sorted_d02[c])

LTG3_sorted_d02 = np.sort(LTG3_d02.flatten())[::-1]
LTG3_c_d02 = LTG3_sorted_d02[0:c]
LTG3_d02_new = np.where(LTG3_d02>=LTG3_sorted_d02[c],LTG3_d02, np.nan)
print('cutoff value LTG3_d02 = ', LTG3_sorted_d02[c])

PR92W_sorted_d02 = np.sort(PR92W_d02.flatten())[::-1]
PR92W_c_d02 = PR92W_sorted_d02[0:c]
PR92W_d02_new = np.where(PR92W_d02>=PR92W_sorted_d02[c],PR92W_d02, np.nan)
print('cutoff value PR92W_d02 = ', PR92W_sorted_d02[c])

CAPExP_CSI_sorted_d02 = np.sort(CAPExP_CSI_d02.flatten())[::-1]
CAPExP_CSI_c_d02 = CAPExP_CSI_sorted_d02[0:c]
CAPExP_CSI_d02_new = np.where(CAPExP_CSI_d02>=CAPExP_CSI_sorted_d02[c],CAPExP_CSI_d02, np.nan)
print('cutoff value CAPExP_CSI_d02 = ', CAPExP_CSI_sorted_d02[c])

CAPExP_R_sorted_d02 = np.sort(CAPExP_R_d02.flatten())[::-1]
CAPExP_R_c_d02 = CAPExP_R_sorted_d02[0:c]
CAPExP_R_d02_new = np.where(CAPExP_R_d02>=CAPExP_R_sorted_d02[c],CAPExP_R_d02, np.nan)
print('cutoff value CAPExP_R_d02 = ', CAPExP_R_sorted_d02[c])

# -------------------------------------------------------------------------------------------------------------
# X_adj = aX+b if X >= c
# -------------------------------------------------------------------------------------------------------------
reg_LPI_d01 = lr().fit(LPI_c_d01.reshape((-1,1)), Obs)
LPI_d01_adj =np.add(reg_LPI_d01.intercept_, np.multiply(reg_LPI_d01.coef_,LPI_d01_new))
LPI_d01_adj[np.isnan(LPI_d01_adj)] = 0

reg_LTG3_d01 = lr().fit(LTG3_c_d01.reshape((-1,1)), Obs)
LTG3_d01_adj =np.add(reg_LTG3_d01.intercept_, np.multiply(reg_LTG3_d01.coef_,LTG3_d01_new))
LTG3_d01_adj[np.isnan(LTG3_d01_adj)] = 0

reg_PR92W_d01 = lr().fit(PR92W_c_d01.reshape((-1,1)), Obs)
PR92W_d01_adj =np.add(reg_PR92W_d01.intercept_, np.multiply(reg_PR92W_d01.coef_,PR92W_d01_new))
PR92W_d01_adj[np.isnan(PR92W_d01_adj)] = 0

reg_CAPExP_CSI_d01 = lr().fit(CAPExP_CSI_c_d01.reshape((-1,1)), Obs)
CAPExP_CSI_d01_adj =np.add(reg_CAPExP_CSI_d01.intercept_, np.multiply(reg_CAPExP_CSI_d01.coef_,CAPExP_CSI_d01_new))
CAPExP_CSI_d01_adj[np.isnan(CAPExP_CSI_d01_adj)] = 0

reg_CAPExP_R_d01 = lr().fit(CAPExP_R_c_d01.reshape((-1,1)), Obs)
CAPExP_R_d01_adj =np.add(reg_CAPExP_R_d01.intercept_, np.multiply(reg_CAPExP_R_d01.coef_,CAPExP_R_d01_new))
CAPExP_R_d01_adj[np.isnan(CAPExP_R_d01_adj)] = 0

reg_LPI_d02 = lr().fit(LPI_c_d02.reshape((-1,1)), Obs)
LPI_d02_adj =np.add(reg_LPI_d02.intercept_, np.multiply(reg_LPI_d02.coef_,LPI_d02_new))
LPI_d02_adj[np.isnan(LPI_d02_adj)] = 0

reg_LTG3_d02 = lr().fit(LTG3_c_d02.reshape((-1,1)), Obs)
LTG3_d02_adj =np.add(reg_LTG3_d02.intercept_, np.multiply(reg_LTG3_d02.coef_,LTG3_d02_new))
LTG3_d02_adj[np.isnan(LTG3_d02_adj)] = 0

reg_PR92W_d02 = lr().fit(PR92W_c_d02.reshape((-1,1)), Obs)
PR92W_d02_adj =np.add(reg_PR92W_d02.intercept_, np.multiply(reg_PR92W_d02.coef_,PR92W_d02_new))
PR92W_d02_adj[np.isnan(PR92W_d02_adj)] = 0

reg_CAPExP_CSI_d02 = lr().fit(CAPExP_CSI_c_d02.reshape((-1,1)), Obs)
CAPExP_CSI_d02_adj =np.add(reg_CAPExP_CSI_d02.intercept_, np.multiply(reg_CAPExP_CSI_d02.coef_,CAPExP_CSI_d02_new))
CAPExP_CSI_d02_adj[np.isnan(CAPExP_CSI_d02_adj)] = 0

reg_CAPExP_R_d02 = lr().fit(CAPExP_R_c_d02.reshape((-1,1)), Obs)
CAPExP_R_d02_adj =np.add(reg_CAPExP_R_d02.intercept_, np.multiply(reg_CAPExP_R_d02.coef_,CAPExP_R_d02_new))
CAPExP_R_d02_adj[np.isnan(CAPExP_R_d02_adj)] = 0


# Cutoff lowest values to match number of flashes between obs and model
Tot_flashes = np.nansum(Obs)

LPI_adj_sorted_d01 = np.sort(LPI_d01_adj.flatten())[::-1]
LTG3_adj_sorted_d01 = np.sort(LTG3_d01_adj.flatten())[::-1]
PR92W_adj_sorted_d01 = np.sort(PR92W_d01_adj.flatten())[::-1]
CAPExP_CSI_adj_sorted_d01 = np.sort(CAPExP_CSI_d01_adj.flatten())[::-1]
CAPExP_R_adj_sorted_d01 = np.sort(CAPExP_R_d01_adj.flatten())[::-1]
LPI_adj_sorted_d02 = np.sort(LPI_d02_adj.flatten())[::-1]
LTG3_adj_sorted_d02 = np.sort(LTG3_d02_adj.flatten())[::-1]
PR92W_adj_sorted_d02 = np.sort(PR92W_d02_adj.flatten())[::-1]
CAPExP_CSI_adj_sorted_d02 = np.sort(CAPExP_CSI_d02_adj.flatten())[::-1]
CAPExP_R_adj_sorted_d02 = np.sort(CAPExP_R_d02_adj.flatten())[::-1]

sum_lpi_d01 = 0
sum_ltg3_d01 = 0
sum_pr92w_d01 = 0
sum_capexp_csi_d01 = 0
sum_capexp_r_d01 = 0
sum_lpi_d02 = 0
sum_ltg3_d02 = 0
sum_pr92w_d02 = 0
sum_capexp_csi_d02 = 0
sum_capexp_r_d02 = 0
i = 0

while sum_lpi_d01 < Tot_flashes:
    sum_lpi_d01 = sum_lpi_d01 + LPI_adj_sorted_d01[i]
    i += 1
else:
    LPI_adj_sorted_d01[i:] = 0
    LPI_adj_sorted_d01.reshape((13248,63,109))
i = 0
while sum_ltg3_d01 < Tot_flashes:
    sum_ltg3_d01 = sum_ltg3_d01 + LTG3_adj_sorted_d01[i]
    i += 1
else:
    LTG3_adj_sorted_d01[i:] = 0
    LTG3_adj_sorted_d01.reshape((13248, 63, 109))
i=0
while sum_pr92w_d01 < Tot_flashes:
    sum_pr92w_d01 = sum_pr92w_d01 + PR92W_adj_sorted_d01[i]
    i += 1
else:
    PR92W_adj_sorted_d01[i:] = 0
    PR92W_adj_sorted_d01.reshape((13248, 63, 109))
i=0
while sum_capexp_csi_d01 < Tot_flashes:
    sum_capexp_csi_d01 = sum_capexp_csi_d01 + CAPExP_CSI_adj_sorted_d01[i]
    i += 1
else:
    CAPExP_CSI_adj_sorted_d01[i:] = 0
    CAPExP_CSI_adj_sorted_d01.reshape((13248, 63, 109))
i=0
while sum_capexp_r_d01 < Tot_flashes:
    sum_capexp_r_d01 = sum_capexp_r_d01 + CAPExP_R_adj_sorted_d01[i]
    i += 1
else:
    CAPExP_R_adj_sorted_d01[i:] = 0
    CAPExP_R_adj_sorted_d01.reshape((13248, 63, 109))
i=0
while sum_lpi_d02 < Tot_flashes:
    sum_lpi_d02 = sum_lpi_d02 + LPI_adj_sorted_d02[i]
    i += 1
else:
    LPI_adj_sorted_d02[i:] = 0
    LPI_adj_sorted_d02.reshape((13248, 63, 109))
i=0
while sum_ltg3_d02 < Tot_flashes:
    sum_ltg3_d02 = sum_ltg3_d02 + LTG3_adj_sorted_d02[i]
    i += 1
else:
    LTG3_adj_sorted_d02[i:] = 0
    LTG3_adj_sorted_d02.reshape((13248, 63, 109))
i=0
while sum_pr92w_d02 < Tot_flashes:
    sum_pr92w_d02 = sum_pr92w_d02 + PR92W_adj_sorted_d02[i]
    i += 1
else:
    PR92W_adj_sorted_d02[i:] = 0
    PR92W_adj_sorted_d02.reshape((13248, 63, 109))
i=0
while sum_capexp_csi_d02 < Tot_flashes:
    sum_capexp_csi_d02 = sum_capexp_csi_d02 + CAPExP_CSI_adj_sorted_d02[i]
    i += 1
else:
    CAPExP_CSI_adj_sorted_d02[i:] = 0
    CAPExP_CSI_adj_sorted_d02.reshape((13248, 63, 109))
i=0
while sum_capexp_r_d02 < Tot_flashes:
    sum_capexp_r_d02 = sum_capexp_r_d02 + CAPExP_R_adj_sorted_d02[i]
    i += 1
else:
    CAPExP_R_adj_sorted_d02[i:] = 0
    CAPExP_R_adj_sorted_d02.reshape((13248, 63, 109))

# ---------------------------------------------------------------------------------------------
# PUT CALIBRATED DATA IN NEW NETCDF FILE
# ---------------------------------------------------------------------------------------------
# Create .nc file
ds = Dataset('/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/wrfout_nc_files/data_calibrated_ax_d01.nc', mode='w', format='NETCDF4')
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
ds['LPI_d01'][:] = LPI_d01_adj[:]
ds['LPI_d02'][:] = LPI_d02_adj[:]
ds['PR92W_d01'][:] = PR92W_d01_adj[:]
ds['PR92W_d02'][:] = PR92W_d02_adj[:]
ds['LTG3_d01'][:] = LTG3_d01_adj[:]
ds['LTG3_d02'][:] = LTG3_d02_adj[:]
ds['CAPExP_CSI_d01'][:] = CAPExP_CSI_d01_adj[:]
ds['CAPExP_CSI_d02'][:] = CAPExP_CSI_d02_adj[:]
ds['CAPExP_R_d01'][:] = CAPExP_R_d01_adj[:]
ds['CAPExP_R_d02'][:] = CAPExP_R_d02_adj[:]
ds['Obs'][:] = ds_obs['Flashdensity_CC'][:,:,:].data + ds_obs['Flashdensity_CG'][:,:,:].data
ds.close()

# -------------------------------------------------------------------------------------------------------------
# Plots to check
# -------------------------------------------------------------------------------------------------------------
LPI_adj_sorted_d01[LPI_adj_sorted_d01 == 0] = np.nan
LPI_adj_sorted_d02[LPI_adj_sorted_d02 == 0] = np.nan
LTG3_adj_sorted_d01[LTG3_adj_sorted_d01 == 0] = np.nan
LTG3_adj_sorted_d02[LTG3_adj_sorted_d02 == 0] = np.nan
PR92W_adj_sorted_d01[PR92W_adj_sorted_d01 == 0] = np.nan
PR92W_adj_sorted_d02[PR92W_adj_sorted_d02 == 0] = np.nan
CAPExP_CSI_adj_sorted_d01[CAPExP_CSI_adj_sorted_d01 == 0] = np.nan
CAPExP_CSI_adj_sorted_d02[CAPExP_CSI_adj_sorted_d02 == 0] = np.nan
CAPExP_R_adj_sorted_d02[CAPExP_R_adj_sorted_d02 == 0] = np.nan
CAPExP_R_adj_sorted_d01[CAPExP_R_adj_sorted_d01 == 0] = np.nan

# ---------------------------------------------------------------------------------------------
# FREQUENCY IFO HOURLY GRID FLASH DENSITY
# ---------------------------------------------------------------------------------------------
# specify bin edges
my_bins = np.linspace(-3,2,20)
# DOMAIN 1
n, bins_log10, patches = plt.hist(np.log10(L[L!=0]),bins = my_bins)
plt.close()
n_CAPExP_R, bins_CAPExP_R_log10, patches_CAPExP_R = plt.hist(np.log10(CAPExP_R_adj_sorted_d01[CAPExP_R_adj_sorted_d01!=0]),bins=my_bins)
plt.close()
n_CAPExP_CSI, bins_CAPExP_CSI_log10, patches_CAPExP_CSI = plt.hist(np.log10(CAPExP_CSI_adj_sorted_d01[CAPExP_CSI_adj_sorted_d01!=0]),bins=my_bins)
plt.close()
n_lpi, bins_lpi_log10, patches_lpi = plt.hist(np.log10(LPI_adj_sorted_d01[LPI_adj_sorted_d01!=0]),bins=my_bins)
plt.close()
n_LTG3, bins_LTG3_log10, patches_LTG3 = plt.hist(np.log10(LTG3_adj_sorted_d01[LTG3_adj_sorted_d01!=0]),bins=my_bins)
plt.close()
n_PR92W, bins_PR92W_log10, patches_PR92W = plt.hist(np.log10(PR92W_adj_sorted_d01[PR92W_adj_sorted_d01!=0]),bins=my_bins)
plt.close()

# set bins with 0 to the left to np.nan for a line plot, lines should just end to the left when there are no values
for i in range(len(n)):
    if n[i]==0:
        n[i] = np.nan
    else:
        break
for i in range(len(n)):
    if n_CAPExP_R[i] == 0:
        n_CAPExP_R[i] = np.nan
    else:
        break
for i in range(len(n)):
    if n_CAPExP_CSI[i] == 0:
        n_CAPExP_CSI[i] = np.nan
    else:
        break
for i in range(len(n)):
    if n_lpi[i] == 0:
        n_lpi[i] = np.nan
    else:
        break
for i in range(len(n)):
    if n_LTG3[i] == 0:
        n_LTG3[i] = np.nan
    else:
        break
for i in range(len(n)):
    if n_PR92W[i] == 0:
        n_PR92W[i] = np.nan
    else:
        break

# ---------------------------------------------------------------------------------------------
# PLOTS
# ---------------------------------------------------------------------------------------------
common_xmin = 0.003
common_ymin = 1
common_ymax = 3000000

x1 =0.004

fig = plt.figure()
ax1 = plt.subplot2grid((3,2),(0,0))
ax2 = plt.subplot2grid((3,2),(2,0))
ax3 = plt.subplot2grid((3,2),(0,1))
ax4 = plt.subplot2grid((3,2),(1,0))
ax5 = plt.subplot2grid((3,2),(2,1))
ax6 = plt.subplot2grid((3,2),(1,1))

bins_log10_centered = bins_log10[0:-1]+0.5*(bins_log10[1:]-bins_log10[0:-1])
bins_centered = 10**bins_log10_centered
ax1.loglog(bins_centered,n,'grey')
ax1.set_title('CLDN')
ax1.grid(which='major', axis='both', color='lightgray')
ax1.set_ylabel('Frequency')
ax1.annotate('(a) \n ', xy=(x1, ax1.get_ylim()[1]),annotation_clip=False)

bins_CAPExP_R_log10_centered = bins_CAPExP_R_log10[0:-1]+0.5*(bins_CAPExP_R_log10[1:]-bins_CAPExP_R_log10[0:-1])
bins_CAPExP_R_centered = 10**bins_CAPExP_R_log10_centered
ax2.loglog(bins_centered,n,'grey')
ax2.loglog(bins_CAPExP_R_centered,n_CAPExP_R,'k')
ax2.set_title('CAPExP_R')
ax2.grid(which='major', axis='both', color='lightgray')
ax2.set_ylabel('Frequency')
ax2.set_xlabel('Hourly flash density (# hr$^{-1}$ km$^{-2}$)')
ax2.annotate('(e) \n ', xy=(x1, ax2.get_ylim()[1]),annotation_clip=False)

bins_lpi_log10_centered = bins_lpi_log10[0:-1]+0.5*(bins_lpi_log10[1:]-bins_lpi_log10[0:-1])
bins_lpi_centered = 10**bins_lpi_log10_centered
ax3.loglog(bins_centered,n,'grey')
ax3.loglog(bins_lpi_centered,n_lpi,'k')
ax3.set_title('LPI')
ax3.grid(which='major', axis='both', color='lightgray')
ax3.set_ylabel('Frequency')
ax3.annotate('(b) \n ', xy=(x1, ax3.get_ylim()[1]),annotation_clip=False)

bins_LTG3_log10_centered = bins_LTG3_log10[0:-1]+0.5*(bins_LTG3_log10[1:]-bins_LTG3_log10[0:-1])
bins_LTG3_centered = 10**bins_LTG3_log10_centered
ax4.loglog(bins_centered,n,'grey')
ax4.loglog(bins_LTG3_centered,n_LTG3,'k')
ax4.set_title('LT3')
ax4.grid(which='major', axis='both', color='lightgray')
ax4.set_ylabel('Frequency')
ax4.annotate('(c) \n ', xy=(x1, ax4.get_ylim()[1]),annotation_clip=False)

bins_CAPExP_CSI_log10_centered = bins_CAPExP_CSI_log10[0:-1]+0.5*(bins_CAPExP_CSI_log10[1:]-bins_CAPExP_CSI_log10[0:-1])
bins_CAPExP_CSI_centered = 10**bins_CAPExP_CSI_log10_centered
ax5.loglog(bins_centered,n,'grey')
ax5.loglog(bins_CAPExP_CSI_centered,n_CAPExP_CSI,'k')
ax5.set_title('CAPExP_CSI')
ax5.grid(which='major', axis='both', color='lightgray')
ax5.set_ylabel('Frequency')
ax5.set_xlabel('Hourly flash density (# hr$^{-1}$ km$^{-2}$)')
ax5.annotate('(f) \n ', xy=(x1, ax5.get_ylim()[1]),annotation_clip=False)

bins_PR92W_log10_centered = bins_PR92W_log10[0:-1]+0.5*(bins_PR92W_log10[1:]-bins_PR92W_log10[0:-1])
bins_PR92W_centered = 10**bins_PR92W_log10_centered
ax6.loglog(bins_centered,n,'grey')
ax6.loglog(bins_PR92W_centered,n_PR92W,'k')
ax6.set_title('PR92W')
ax6.grid(which='major', axis='both', color='lightgray')
ax6.set_ylabel('Frequency')
ax6.annotate('(d) \n ', xy=(x1, ax6.get_ylim()[1]),annotation_clip=False)

plt.suptitle('Convection-parameterized (9 km)')
plt.show()

# DOMAIN 2
n, bins_log10, patches = plt.hist(np.log10(L[L!=0]),bins = my_bins)
plt.close()
n_CAPExP_R_d02, bins_CAPExP_R_log10_d02, patches_CAPExP_R_d02 = plt.hist(np.log10(CAPExP_R_adj_sorted_d02[CAPExP_R_adj_sorted_d02!=0]),bins=my_bins)
plt.close()
n_CAPExP_CSI_d02, bins_CAPExP_CSI_log10_d02, patches_CAPExP_CSI_d02 = plt.hist(np.log10(CAPExP_CSI_adj_sorted_d02[CAPExP_CSI_adj_sorted_d02!=0]),bins=my_bins)
plt.close()
n_lpi_d02, bins_lpi_log10_d02, patches_lpi_d02 = plt.hist(np.log10(LPI_adj_sorted_d02[LPI_adj_sorted_d02!=0]),bins=my_bins)
plt.close()
n_LTG3_d02, bins_LTG3_log10_d02, patches_LTG3_d02 = plt.hist(np.log10(LTG3_adj_sorted_d02[LTG3_adj_sorted_d02!=0]),bins=my_bins)
plt.close()
n_PR92W_d02, bins_PR92W_log10_d02, patches_PR92W_d02 = plt.hist(np.log10(PR92W_adj_sorted_d02[PR92W_adj_sorted_d02!=0]),bins=my_bins)
plt.close()

# set bins with 0 to the left to np.nan for a line plot, lines should just end to the left when there are no values
for i in range(len(n)):
    if n[i]==0:
        n[i] = np.nan
    else:
        break
for i in range(len(n)):
    if n_CAPExP_R_d02[i] == 0:
        n_CAPExP_R_d02[i] = np.nan
    else:
        break
for i in range(len(n)):
    if n_CAPExP_CSI_d02[i] == 0:
        n_CAPExP_CSI_d02[i] = np.nan
    else:
        break
for i in range(len(n)):
    if n_lpi_d02[i] == 0:
        n_lpi_d02[i] = np.nan
    else:
        break
for i in range(len(n)):
    if n_LTG3_d02[i] == 0:
        n_LTG3_d02[i] = np.nan
    else:
        break
for i in range(len(n)):
    if n_PR92W_d02[i] == 0:
        n_PR92W_d02[i] = np.nan
    else:
        break

# ---------------------------------------------------------------------------------------------
# PLOTS
# ---------------------------------------------------------------------------------------------
common_xmin = 0.003
common_ymin = 1
common_ymax = 3000000

fig = plt.figure()
ax1 = plt.subplot2grid((3,2),(0,0))
ax2 = plt.subplot2grid((3,2),(2,0))
ax3 = plt.subplot2grid((3,2),(0,1))
ax4 = plt.subplot2grid((3,2),(1,0))
ax5 = plt.subplot2grid((3,2),(2,1))
ax6 = plt.subplot2grid((3,2),(1,1))

bins_log10_centered = bins_log10[0:-1]+0.5*(bins_log10[1:]-bins_log10[0:-1])
bins_centered = 10**bins_log10_centered
ax1.loglog(bins_centered,n,'grey')
ax1.set_title('CLDN')
ax1.grid(which='major', axis='both', color='lightgray')
ax1.set_ylabel('Frequency')
ax1.annotate('(a) \n ', xy=(x1, ax1.get_ylim()[1]),annotation_clip=False)

bins_CAPExP_R_log10_centered_d02 = bins_CAPExP_R_log10_d02[0:-1]+0.5*(bins_CAPExP_R_log10_d02[1:]-bins_CAPExP_R_log10_d02[0:-1])
bins_CAPExP_R_centered_d02 = 10**bins_CAPExP_R_log10_centered_d02
ax2.loglog(bins_centered,n,'grey')
ax2.loglog(bins_CAPExP_R_centered_d02,n_CAPExP_R_d02,'k')
ax2.set_title('CAPExP_R')
ax2.grid(which='major', axis='both', color='lightgray')
ax2.set_ylabel('Frequency')
ax2.set_xlabel('Hourly flash density (# hr$^{-1}$ km$^{-2}$)')
ax2.annotate('(e) \n ', xy=(x1,  ax2.get_ylim()[1]),annotation_clip=False)

bins_lpi_log10_centered_d02 = bins_lpi_log10_d02[0:-1]+0.5*(bins_lpi_log10_d02[1:]-bins_lpi_log10_d02[0:-1])
bins_lpi_centered_d02 = 10**bins_lpi_log10_centered_d02
ax3.loglog(bins_centered,n,'grey')
ax3.loglog(bins_lpi_centered_d02,n_lpi_d02,'k')
ax3.set_title('LPI')
ax3.grid(which='major', axis='both', color='lightgray')
ax3.set_ylabel('Frequency')
ax3.annotate('(b) \n ', xy=(x1, ax3.get_ylim()[1]),annotation_clip=False)

bins_LTG3_log10_centered_d02 = bins_LTG3_log10_d02[0:-1]+0.5*(bins_LTG3_log10_d02[1:]-bins_LTG3_log10_d02[0:-1])
bins_LTG3_centered_d02 = 10**bins_LTG3_log10_centered_d02
ax4.loglog(bins_centered,n,'grey')
ax4.loglog(bins_LTG3_centered_d02,n_LTG3_d02,'k')
ax4.set_title('LT3')
ax4.grid(which='major', axis='both', color='lightgray')
ax4.set_ylabel('Frequency')
ax4.annotate('(c) \n ', xy=(x1, ax4.get_ylim()[1]),annotation_clip=False)

bins_CAPExP_CSI_log10_centered_d02 = bins_CAPExP_CSI_log10_d02[0:-1]+0.5*(bins_CAPExP_CSI_log10_d02[1:]-bins_CAPExP_CSI_log10_d02[0:-1])
bins_CAPExP_CSI_centered_d02 = 10**bins_CAPExP_CSI_log10_centered_d02
ax5.loglog(bins_centered,n,'grey')
ax5.loglog(bins_CAPExP_CSI_centered_d02,n_CAPExP_CSI_d02,'k')
ax5.set_title('CAPExP_CSI')
ax5.grid(which='major', axis='both', color='lightgray')
ax5.set_ylabel('Frequency')
ax5.set_xlabel('Hourly flash density (# hr$^{-1}$ km$^{-2}$)')
ax5.annotate('(f) \n ', xy=(x1,  ax5.get_ylim()[1]),annotation_clip=False)

bins_PR92W_log10_centered_d02 = bins_PR92W_log10_d02[0:-1]+0.5*(bins_PR92W_log10_d02[1:]-bins_PR92W_log10_d02[0:-1])
bins_PR92W_centered_d02 = 10**bins_PR92W_log10_centered_d02
ax6.loglog(bins_centered,n,'grey')
ax6.loglog(bins_PR92W_centered_d02,n_PR92W_d02,'k')
ax6.set_title('PR92W')
ax6.grid(which='major', axis='both', color='lightgray')
ax6.set_ylabel('Frequency')
ax6.annotate('(d) \n ', xy=(x1, ax6.get_ylim()[1]),annotation_clip=False)

plt.suptitle('Convection-permitting (3 km)')
plt.show()
