# ---------------------------------------------------------------------------------------------
# MODULES
# ---------------------------------------------------------------------------------------------
import netCDF4 as nc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from netCDF4 import Dataset
# ---------------------------------------------------------------------------------------------
# LOAD DATA
# ---------------------------------------------------------------------------------------------
ds_d01 = Dataset('/scratch/leuven/336/vsc33651/nu-wrf-dev/domain1_all.nc','r')
ds_d02 = Dataset('/scratch/leuven/336/vsc33651/nu-wrf-dev/domain2_at_domain1_all.nc','r')
ds_obs = Dataset('/scratch/leuven/336/vsc33651/nu-wrf-dev/wrfout_nc_files/CLDN_at_domain1_all.nc')

L = ds_obs['Flashdensity_CC'][:,:,:].data.flatten() + ds_obs['Flashdensity_CG'][:,:,:].data.flatten()

# ---------------------------------------------------------------------------------------------
# NORMALIZE DATA
# ---------------------------------------------------------------------------------------------
Ln = L/np.amax(L)
print(np.amax(Ln))
LPI_n_d01 = ds_d01['LPI'][:]/np.nanmax(ds_d01['LPI'][:])
print(np.amax(LPI_n_d01))
LPI_n_d02 = ds_d02['LPI'][:]/np.nanmax(ds_d02['LPI'][:])
print(np.amax(LPI_n_d02))
LTG3_n_d01 = ds_d01['LTG3'][:]/np.nanmax(ds_d01['LTG3'][:])
print(np.amax(LTG3_n_d01))
LTG3_n_d02 = ds_d02['LTG3'][:]/np.nanmax(ds_d02['LTG3'][:])
print(np.amax(LTG3_n_d02))
PR92H_n_d01 = ds_d01['PR92H'][:]/np.nanmax(ds_d01['PR92H'][:])
print(np.amax(PR92H_n_d01))
PR92H_n_d02 = ds_d02['PR92H'][:]/np.nanmax(ds_d02['PR92H'][:])
print(np.amax(PR92H_n_d02))
PR92W_n_d01 = ds_d01['PR92W'][:]/np.nanmax(ds_d01['PR92W'][:])
print(np.amax(PR92W_n_d01))
PR92W_n_d02 = ds_d02['PR92W'][:]/np.nanmax(ds_d02['PR92W'][:])
print(np.amax(PR92W_n_d02))
CAPExP_n_d01 = ds_d01['CAPExP'][:]/np.nanmax(ds_d01['CAPExP'][:])
print(np.amax(CAPExP_n_d01))

# ---------------------------------------------------------------------------------------------
# FREQUENCY IFO HOURLY GRID FLASH DENSITY
# ---------------------------------------------------------------------------------------------
# DOMAIN 1
n, bins_log10, patches = plt.hist(np.log10(Ln[Ln!=0]),10)
plt.close()
n_CAPExP, bins_CAPExP_log10, patches_CAPExP = plt.hist(np.log10(CAPExP_n_d01[CAPExP_n_d01!=0]),10)
plt.close()
n_lpi, bins_lpi_log10, patches_lpi = plt.hist(np.log10(LPI_n_d01[LPI_n_d01!=0]),10)
plt.close()
n_LTG3, bins_LTG3_log10, patches_LTG3 = plt.hist(np.log10(LTG3_n_d01[LTG3_n_d01!=0]),10)
plt.close()
n_PR92H, bins_PR92H_log10, patches_PR92H = plt.hist(np.log10(PR92H_n_d01[PR92H_n_d01!=0]),10)
plt.close()
n_PR92W, bins_PR92W_log10, patches_PR92W = plt.hist(np.log10(PR92W_n_d01[PR92W_n_d01!=0]),10)
plt.close()

fig = plt.figure()
plt.title('Convection-parameterized (9 km)')
ax1 = plt.subplot2grid((3,2),(0,0))
ax2 = plt.subplot2grid((3,2),(0,1))
ax3 = plt.subplot2grid((3,2),(1,0))
ax4 = plt.subplot2grid((3,2),(1,1))
ax5 = plt.subplot2grid((3,2),(2,0))
ax6 = plt.subplot2grid((3,2),(2,1))

bins_log10_centered = bins_log10[0:-1]+0.5*(bins_log10[1:]-bins_log10[0:-1])
bins_centered = 10**bins_log10_centered
ax1.loglog(bins_centered,n,'grey')
ax1.set_title('CLDN')
# ax1.set_ylabel('Frequency')
# ax1.set_xlabel('Normalized hourly grid flash density (#/km$[2}$)')

bins_CAPExP_log10_centered = bins_CAPExP_log10[0:-1]+0.5*(bins_CAPExP_log10[1:]-bins_CAPExP_log10[0:-1])
bins_CAPExP_centered = 10**bins_CAPExP_log10_centered
ax2.loglog(bins_centered,n,'grey')
ax2.loglog(bins_CAPExP_centered,n_CAPExP,'k')
ax2.set_title('CAPExP')
# ax2.set_ylabel('Frequency')
# ax2.set_xlabel('Normalized hourly grid flash density (#/km$[2}$)')

bins_lpi_log10_centered = bins_lpi_log10[0:-1]+0.5*(bins_lpi_log10[1:]-bins_lpi_log10[0:-1])
bins_lpi_centered = 10**bins_lpi_log10_centered
ax3.loglog(bins_centered,n,'grey')
ax3.loglog(bins_lpi_centered,n_lpi,'k')
ax3.set_title('LPI')
# ax3.set_ylabel('Frequency')
# ax3.set_xlabel('Normalized hourly grid flash density (#/km$[2}$)')

bins_LTG3_log10_centered = bins_LTG3_log10[0:-1]+0.5*(bins_LTG3_log10[1:]-bins_LTG3_log10[0:-1])
bins_LTG3_centered = 10**bins_LTG3_log10_centered
ax4.loglog(bins_centered,n,'grey')
ax4.loglog(bins_LTG3_centered,n_LTG3,'k')
ax4.set_title('LT3')
# ax4.set_ylabel('Frequency')
# ax4.set_xlabel('Normalized hourly grid flash density (#/km$[2}$)')

bins_PR92H_log10_centered = bins_PR92H_log10[0:-1]+0.5*(bins_PR92H_log10[1:]-bins_PR92H_log10[0:-1])
bins_PR92H_centered = 10**bins_PR92H_log10_centered
ax5.loglog(bins_centered,n,'grey')
ax5.loglog(bins_PR92H_centered,n_PR92H,'k')
ax5.set_title('PR92H')
# ax5.set_ylabel('Frequency')
# ax5.set_xlabel('Normalized hourly grid flash density (#/km$[2}$)')

bins_PR92W_log10_centered = bins_PR92W_log10[0:-1]+0.5*(bins_PR92W_log10[1:]-bins_PR92W_log10[0:-1])
bins_PR92W_centered = 10**bins_PR92W_log10_centered
ax6.loglog(bins_centered,n,'grey')
ax6.loglog(bins_PR92W_centered,n_PR92W,'k')
ax6.set_title('PR92W')
# ax6.set_ylabel('Frequency')
# ax6.set_xlabel('Normalized hourly grid flash density (#/km$[2}$)')
plt.show()

# DOMAIN 2
n, bins_log10, patches = plt.hist(np.log10(Ln[Ln!=0]),10)
plt.close()
n_CAPExP, bins_CAPExP_log10, patches_CAPExP = plt.hist(np.log10(CAPExP_n_d01[CAPExP_n_d01!=0]),10)
plt.close()
n_lpi_d02, bins_lpi_log10_d02, patches_lpi_d02 = plt.hist(np.log10(LPI_n_d02[LPI_n_d02!=0]),10)
plt.close()
n_LTG3_d02, bins_LTG3_log10_d02, patches_LTG3_d02 = plt.hist(np.log10(LTG3_n_d02[LTG3_n_d02!=0]),10)
plt.close()
n_PR92H_d02, bins_PR92H_log10_d02, patches_PR92H_d02 = plt.hist(np.log10(PR92H_n_d02[PR92H_n_d02!=0]),10)
plt.close()
n_PR92W_d02, bins_PR92W_log10_d02, patches_PR92W_d02 = plt.hist(np.log10(PR92W_n_d02[PR92W_n_d02!=0]),10)
plt.close()

fig = plt.figure()
plt.title('Convection-parameterized (9 km)')
ax1 = plt.subplot2grid((3,2),(0,0))
ax2 = plt.subplot2grid((3,2),(0,1))
ax3 = plt.subplot2grid((3,2),(1,0))
ax4 = plt.subplot2grid((3,2),(1,1))
ax5 = plt.subplot2grid((3,2),(2,0))
ax6 = plt.subplot2grid((3,2),(2,1))

bins_log10_centered = bins_log10[0:-1]+0.5*(bins_log10[1:]-bins_log10[0:-1])
bins_centered = 10**bins_log10_centered
ax1.loglog(bins_centered,n,'grey')
ax1.set_title('CLDN')
# ax1.set_ylabel('Frequency')
# ax1.set_xlabel('Normalized hourly grid flash density (#/km$[2}$)')

bins_CAPExP_log10_centered = bins_CAPExP_log10[0:-1]+0.5*(bins_CAPExP_log10[1:]-bins_CAPExP_log10[0:-1])
bins_CAPExP_centered = 10**bins_CAPExP_log10_centered
ax2.loglog(bins_centered,n,'grey')
ax2.loglog(bins_CAPExP_centered,n_CAPExP,'k')
ax2.set_title('CAPExP')
# ax2.set_ylabel('Frequency')
# ax2.set_xlabel('Normalized hourly grid flash density (#/km$[2}$)')

bins_lpi_log10_centered_d02 = bins_lpi_log10_d02[0:-1]+0.5*(bins_lpi_log10_d02[1:]-bins_lpi_log10_d02[0:-1])
bins_lpi_centered_d02 = 10**bins_lpi_log10_centered_d02
ax3.loglog(bins_centered,n,'grey')
ax3.loglog(bins_lpi_centered_d02,n_lpi_d02,'k')
ax3.set_title('LPI')
# ax3.set_ylabel('Frequency')
# ax3.set_xlabel('Normalized hourly grid flash density (#/km$[2}$)')

bins_LTG3_log10_centered_d02 = bins_LTG3_log10_d02[0:-1]+0.5*(bins_LTG3_log10_d02[1:]-bins_LTG3_log10_d02[0:-1])
bins_LTG3_centered_d02 = 10**bins_LTG3_log10_centered_d02
ax4.loglog(bins_centered,n,'grey')
ax4.loglog(bins_LTG3_centered_d02,n_LTG3_d02,'k')
ax4.set_title('LT3')
# ax4.set_ylabel('Frequency')
# ax4.set_xlabel('Normalized hourly grid flash density (#/km$[2}$)')

bins_PR92H_log10_centered_d02 = bins_PR92H_log10_d02[0:-1]+0.5*(bins_PR92H_log10_d02[1:]-bins_PR92H_log10_d02[0:-1])
bins_PR92H_centered_d02 = 10**bins_PR92H_log10_centered_d02
ax5.loglog(bins_centered,n,'grey')
ax5.loglog(bins_PR92H_centered_d02,n_PR92H_d02,'k')
ax5.set_title('PR92H')
# ax5.set_ylabel('Frequency')
# ax5.set_xlabel('Normalized hourly grid flash density (#/km$[2}$)')

bins_PR92W_log10_centered_d02 = bins_PR92W_log10_d02[0:-1]+0.5*(bins_PR92W_log10_d02[1:]-bins_PR92W_log10_d02[0:-1])
bins_PR92W_centered_d02 = 10**bins_PR92W_log10_centered_d02
ax6.loglog(bins_centered,n,'grey')
ax6.loglog(bins_PR92W_centered_d02,n_PR92W_d02,'k')
ax6.set_title('PR92W')
# ax6.set_ylabel('Frequency')
# ax6.set_xlabel('Normalized hourly grid flash density (#/km$[2}$)')
plt.show()
