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
print('cutoff value LPI_d01 = ', LPI_sorted_d01[c])

LTG3_sorted_d01 = np.sort(LTG3_d01.flatten())[::-1]
LTG3_c_d01 = LTG3_sorted_d01[0:c]
print('cutoff value LTG3_d01 = ', LTG3_sorted_d01[c])

PR92W_sorted_d01 = np.sort(PR92W_d01.flatten())[::-1]
PR92W_c_d01 = PR92W_sorted_d01[0:c]
print('cutoff value PR92W_d01 = ', PR92W_sorted_d01[c])

CAPExP_CSI_sorted_d01 = np.sort(CAPExP_CSI_d01.flatten())[::-1]
CAPExP_CSI_c_d01 = CAPExP_CSI_sorted_d01[0:c]
print('cutoff value CAPExP_CSI_d01 = ', CAPExP_CSI_sorted_d01[c])

CAPExP_R_sorted_d01 = np.sort(CAPExP_R_d01.flatten())[::-1]
CAPExP_R_c_d01 = CAPExP_R_sorted_d01[0:c]
print('cutoff value CAPExP_R_d01 = ', CAPExP_R_sorted_d01[c])

LPI_sorted_d02 = np.sort(LPI_d02.flatten())[::-1]
LPI_c_d02 = LPI_sorted_d02[0:c]
print('cutoff value LPI_d02 = ', LPI_sorted_d02[c])

LTG3_sorted_d02 = np.sort(LTG3_d02.flatten())[::-1]
LTG3_c_d02 = LTG3_sorted_d02[0:c]
print('cutoff value LTG3_d02 = ', LTG3_sorted_d02[c])

PR92W_sorted_d02 = np.sort(PR92W_d02.flatten())[::-1]
PR92W_c_d02 = PR92W_sorted_d02[0:c]
print('cutoff value PR92W_d02 = ', PR92W_sorted_d02[c])

CAPExP_CSI_sorted_d02 = np.sort(CAPExP_CSI_d02.flatten())[::-1]
CAPExP_CSI_c_d02 = CAPExP_CSI_sorted_d02[0:c]
print('cutoff value CAPExP_CSI_d02 = ', CAPExP_CSI_sorted_d02[c])

CAPExP_R_sorted_d02 = np.sort(CAPExP_R_d02.flatten())[::-1]
CAPExP_R_c_d02 = CAPExP_R_sorted_d02[0:c]
print('cutoff value CAPExP_R_d02 = ', CAPExP_R_sorted_d02[c])

# -------------------------------------------------------------------------------------------------------------
# X_adj = aX+b if X >= c
# -------------------------------------------------------------------------------------------------------------
reg_LPI_d01 = lr().fit(LPI_c_d01.reshape((-1,1)), Obs)
LPI_d01_adj =np.add(reg_LPI_d01.intercept_, np.multiply(reg_LPI_d01.coef_,LPI_d01[LPI_d01>=LPI_sorted_d01[c]]))

reg_LTG3_d01 = lr().fit(LTG3_c_d01.reshape((-1,1)), Obs)
LTG3_d01_adj =np.add(reg_LTG3_d01.intercept_, np.multiply(reg_LTG3_d01.coef_,LTG3_d01[LTG3_d01>=LTG3_sorted_d01[c]]))

reg_PR92W_d01 = lr().fit(PR92W_c_d01.reshape((-1,1)), Obs)
PR92W_d01_adj =np.add(reg_PR92W_d01.intercept_, np.multiply(reg_PR92W_d01.coef_,PR92W_d01[PR92W_d01>=PR92W_sorted_d01[c]]))

reg_CAPExP_CSI_d01 = lr().fit(CAPExP_CSI_c_d01.reshape((-1,1)), Obs)
CAPExP_CSI_d01_adj =np.add(reg_CAPExP_CSI_d01.intercept_, np.multiply(reg_CAPExP_CSI_d01.coef_,CAPExP_CSI_d01[CAPExP_CSI_d01>=CAPExP_CSI_sorted_d01[c]]))

reg_CAPExP_R_d01 = lr().fit(CAPExP_R_c_d01.reshape((-1,1)), Obs)
CAPExP_R_d01_adj =np.add(reg_CAPExP_R_d01.intercept_, np.multiply(reg_CAPExP_R_d01.coef_,CAPExP_R_d01[CAPExP_R_d01>=CAPExP_R_sorted_d01[c]]))

reg_LPI_d02 = lr().fit(LPI_c_d02.reshape((-1,1)), Obs)
LPI_d02_adj =np.add(reg_LPI_d02.intercept_, np.multiply(reg_LPI_d02.coef_,LPI_d02[LPI_d02>=LPI_sorted_d02[c]]))

reg_LTG3_d02 = lr().fit(LTG3_c_d02.reshape((-1,1)), Obs)
LTG3_d02_adj =np.add(reg_LTG3_d02.intercept_, np.multiply(reg_LTG3_d02.coef_,LTG3_d02[LTG3_d02>=LTG3_sorted_d02[c]]))

reg_PR92W_d02 = lr().fit(PR92W_c_d02.reshape((-1,1)), Obs)
PR92W_d02_adj =np.add(reg_PR92W_d02.intercept_, np.multiply(reg_PR92W_d02.coef_,PR92W_d02[PR92W_d02>=PR92W_sorted_d02[c]]))

reg_CAPExP_CSI_d02 = lr().fit(CAPExP_CSI_c_d02.reshape((-1,1)), Obs)
CAPExP_CSI_d02_adj =np.add(reg_CAPExP_CSI_d02.intercept_, np.multiply(reg_CAPExP_CSI_d02.coef_,CAPExP_CSI_d02[CAPExP_CSI_d02>=CAPExP_CSI_sorted_d02[c]]))

reg_CAPExP_R_d02 = lr().fit(CAPExP_R_c_d02.reshape((-1,1)), Obs)
CAPExP_R_d02_adj =np.add(reg_CAPExP_R_d02.intercept_, np.multiply(reg_CAPExP_R_d02.coef_,CAPExP_R_d02[CAPExP_R_d02>=CAPExP_R_sorted_d02[c]]))

# -------------------------------------------------------------------------------------------------------------
# Plots to check
# -------------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------
# exclude low values that likely don't produce lightning
# Remove lightning predictions with less than 1 flash per hour per storm (e.g. 100 sqkm)
cond_LTG3_d01 = LTG3_d01_adj < 0.01
cond_LTG3_d02 = LTG3_d02_adj < 0.01
cond_PR92W_d01 = PR92W_d01_adj < 0.01
cond_PR92W_d02 = PR92W_d02_adj < 0.01
cond_CAPExP_CSI_d01 = CAPExP_CSI_d01_adj < 0.01
cond_CAPExP_CSI_d02 = CAPExP_CSI_d01_adj < 0.01
cond_CAPExP_R_d01 = CAPExP_R_d01_adj < 0.01
cond_CAPExP_R_d02 = CAPExP_R_d02_adj < 0.01
# for LPI everything smaller than 0.001
cond_LPI_d01 = LPI_d01_adj<0.001
cond_LPI_d02 = LPI_d02_adj<0.001
# set to nan
LPI_d01_adj[cond_LPI_d01] = np.nan
LPI_d02_adj[cond_LPI_d02] = np.nan
LTG3_d01_adj[cond_LTG3_d01] = np.nan
LTG3_d02_adj[cond_LTG3_d02] = np.nan
PR92W_d01_adj[cond_PR92W_d01] = np.nan
PR92W_d02_adj[cond_PR92W_d02] = np.nan
CAPExP_CSI_d01_adj[cond_CAPExP_CSI_d01] = np.nan
CAPExP_CSI_d01_adj[cond_CAPExP_CSI_d02] = np.nan
CAPExP_R_d01_adj[cond_CAPExP_R_d01] = np.nan
CAPExP_R_d02_adj[cond_CAPExP_R_d02] = np.nan

# ---------------------------------------------------------------------------------------------
# FREQUENCY IFO HOURLY GRID FLASH DENSITY
# ---------------------------------------------------------------------------------------------
# specify bin edges
my_bins = np.linspace(-3,2,20)
# DOMAIN 1
n, bins_log10, patches = plt.hist(np.log10(L[L!=0]),bins = my_bins)
plt.close()
n_CAPExP_R, bins_CAPExP_R_log10, patches_CAPExP_R = plt.hist(np.log10(CAPExP_R_d01_adj[CAPExP_R_d01_adj!=0]),bins=my_bins)
plt.close()
n_CAPExP_CSI, bins_CAPExP_CSI_log10, patches_CAPExP_CSI = plt.hist(np.log10(CAPExP_CSI_d01_adj[CAPExP_CSI_d01_adj!=0]),bins=my_bins)
plt.close()
n_lpi, bins_lpi_log10, patches_lpi = plt.hist(np.log10(LPI_d01_adj[LPI_d01_adj!=0]),bins=my_bins)
plt.close()
n_LTG3, bins_LTG3_log10, patches_LTG3 = plt.hist(np.log10(LTG3_d01_adj[LTG3_d01_adj!=0]),bins=my_bins)
plt.close()
n_PR92W, bins_PR92W_log10, patches_PR92W = plt.hist(np.log10(PR92W_d01_adj[PR92W_d01_adj!=0]),bins=my_bins)
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
# ax1.set_xlim([common_xmin,1])
# ax1.set_ylim([common_ymin,common_ymax])
ax1.grid(which='major', axis='both', color='lightgray')
ax1.set_ylabel('Frequency')
ax1.annotate('(a) \n ', xy=(x1, ax1.get_ylim()[1]),annotation_clip=False)

bins_CAPExP_R_log10_centered = bins_CAPExP_R_log10[0:-1]+0.5*(bins_CAPExP_R_log10[1:]-bins_CAPExP_R_log10[0:-1])
bins_CAPExP_R_centered = 10**bins_CAPExP_R_log10_centered
ax2.loglog(bins_centered,n,'grey')
ax2.loglog(bins_CAPExP_R_centered,n_CAPExP_R,'k')
ax2.set_title('CAPExP_R')
# ax2.set_xlim([common_xmin,1])
# ax2.set_ylim([common_ymin,common_ymax])
ax2.grid(which='major', axis='both', color='lightgray')
ax2.set_ylabel('Frequency')
ax2.set_xlabel('Hourly flash density (# hr$^{-1}$ km$^{-2}$)')
ax2.annotate('(e) \n ', xy=(x1, ax2.get_ylim()[1]),annotation_clip=False)

bins_lpi_log10_centered = bins_lpi_log10[0:-1]+0.5*(bins_lpi_log10[1:]-bins_lpi_log10[0:-1])
bins_lpi_centered = 10**bins_lpi_log10_centered
ax3.loglog(bins_centered,n,'grey')
ax3.loglog(bins_lpi_centered,n_lpi,'k')
ax3.set_title('LPI*')
# ax3.set_xlim([common_xmin,1])
# ax3.set_ylim([common_ymin,common_ymax])
ax3.grid(which='major', axis='both', color='lightgray')
ax3.set_ylabel('Frequency')
ax3.annotate('(b) \n ', xy=(x1, ax3.get_ylim()[1]),annotation_clip=False)

bins_LTG3_log10_centered = bins_LTG3_log10[0:-1]+0.5*(bins_LTG3_log10[1:]-bins_LTG3_log10[0:-1])
bins_LTG3_centered = 10**bins_LTG3_log10_centered
ax4.loglog(bins_centered,n,'grey')
ax4.loglog(bins_LTG3_centered,n_LTG3,'k')
ax4.set_title('LT3')
# ax4.set_xlim([common_xmin,1])
# ax4.set_ylim([common_ymin,common_ymax])
ax4.grid(which='major', axis='both', color='lightgray')
ax4.set_ylabel('Frequency')
ax4.annotate('(c) \n ', xy=(x1, ax4.get_ylim()[1]),annotation_clip=False)

bins_CAPExP_CSI_log10_centered = bins_CAPExP_CSI_log10[0:-1]+0.5*(bins_CAPExP_CSI_log10[1:]-bins_CAPExP_CSI_log10[0:-1])
bins_CAPExP_CSI_centered = 10**bins_CAPExP_CSI_log10_centered
ax5.loglog(bins_centered,n,'grey')
ax5.loglog(bins_CAPExP_CSI_centered,n_CAPExP_CSI,'k')
ax5.set_title('CAPExP_CSI')
# ax5.set_xlim([common_xmin,1])
# ax5.set_ylim([common_ymin,common_ymax])
ax5.grid(which='major', axis='both', color='lightgray')
ax5.set_ylabel('Frequency')
ax5.set_xlabel('Hourly flash density (# hr$^{-1}$ km$^{-2}$) \n * Hourly LPI (J/kg)')
ax5.annotate('(f) \n ', xy=(x1, ax5.get_ylim()[1]),annotation_clip=False)

bins_PR92W_log10_centered = bins_PR92W_log10[0:-1]+0.5*(bins_PR92W_log10[1:]-bins_PR92W_log10[0:-1])
bins_PR92W_centered = 10**bins_PR92W_log10_centered
ax6.loglog(bins_centered,n,'grey')
ax6.loglog(bins_PR92W_centered,n_PR92W,'k')
ax6.set_title('PR92W')
# ax6.set_xlim([common_xmin,1])
# ax6.set_ylim([common_ymin,common_ymax])
ax6.grid(which='major', axis='both', color='lightgray')
ax6.set_ylabel('Frequency')
ax6.annotate('(d) \n ', xy=(x1, ax6.get_ylim()[1]),annotation_clip=False)

plt.suptitle('Convection-parameterized (9 km)')
plt.show()

# DOMAIN 2
n, bins_log10, patches = plt.hist(np.log10(L[L!=0]),bins=my_bins)
plt.close()
n_CAPExP_R, bins_CAPExP_R_log10, patches_CAPExP_R = plt.hist(np.log10(CAPExP_R_d02_adj[CAPExP_R_d02_adj!=0]),bins=my_bins)
plt.close()
n_CAPExP_CSI, bins_CAPExP_CSI_log10, patches_CAPExP_CSI = plt.hist(np.log10(CAPExP_CSI_d02_adj[CAPExP_CSI_d02_adj!=0]),bins=my_bins)
plt.close()
n_lpi_d02, bins_lpi_log10_d02, patches_lpi_d02 = plt.hist(np.log10(LPI_d02_adj[LPI_d02_adj!=0]),bins=my_bins)
plt.close()
n_LTG3_d02, bins_LTG3_log10_d02, patches_LTG3_d02 = plt.hist(np.log10(LTG3_d02_adj[LTG3_d02_adj!=0]),bins=my_bins)
plt.close()
n_PR92W_d02, bins_PR92W_log10_d02, patches_PR92W_d02 = plt.hist(np.log10(PR92W_d02_adj[PR92W_d02_adj!=0]),bins=my_bins)
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
# ax1.set_xlim([common_xmin,1])
# ax1.set_ylim([common_ymin,common_ymax])
ax1.annotate('(a) \n ', xy=(x1, ax1.get_ylim()[1]),annotation_clip=False)

bins_CAPExP_R_log10_centered = bins_CAPExP_R_log10[0:-1]+0.5*(bins_CAPExP_R_log10[1:]-bins_CAPExP_R_log10[0:-1])
bins_CAPExP_R_centered = 10**bins_CAPExP_R_log10_centered
ax2.loglog(bins_centered,n,'grey')
ax2.loglog(bins_CAPExP_R_centered,n_CAPExP_R,'k')
ax2.set_title('CAPExP_R')
# ax2.set_xlim([common_xmin,1])
# ax2.set_ylim([common_ymin,common_ymax])
ax2.grid(which='major', axis='both', color='lightgray')
ax2.set_ylabel('Frequency')
ax2.set_xlabel('Hourly flash density (# hr$^{-1}$ km$^{-2}$)')
ax2.annotate('(e) \n ', xy=(x1,  ax2.get_ylim()[1]),annotation_clip=False)

bins_lpi_log10_centered_d02 = bins_lpi_log10_d02[0:-1]+0.5*(bins_lpi_log10_d02[1:]-bins_lpi_log10_d02[0:-1])
bins_lpi_centered_d02 = 10**bins_lpi_log10_centered_d02
ax3.loglog(bins_centered,n,'grey')
ax3.loglog(bins_lpi_centered_d02,n_lpi_d02,'k')
ax3.set_title('LPI*')
ax3.grid(which='major', axis='both', color='lightgray')
ax3.set_ylabel('Frequency')
# ax3.set_xlim([common_xmin,1])
# ax3.set_ylim([common_ymin,common_ymax])
ax3.annotate('(b) \n ', xy=(x1, ax3.get_ylim()[1]),annotation_clip=False)

bins_LTG3_log10_centered_d02 = bins_LTG3_log10_d02[0:-1]+0.5*(bins_LTG3_log10_d02[1:]-bins_LTG3_log10_d02[0:-1])
bins_LTG3_centered_d02 = 10**bins_LTG3_log10_centered_d02
ax4.loglog(bins_centered,n,'grey')
ax4.loglog(bins_LTG3_centered_d02,n_LTG3_d02,'k')
ax4.set_title('LT3')
ax4.grid(which='major', axis='both', color='lightgray')
ax4.set_ylabel('Frequency')
# ax4.set_xlim([common_xmin,1])
# ax4.set_ylim([common_ymin,common_ymax])
ax4.annotate('(c) \n ', xy=(x1, ax4.get_ylim()[1]),annotation_clip=False)

bins_CAPExP_CSI_log10_centered = bins_CAPExP_CSI_log10[0:-1]+0.5*(bins_CAPExP_CSI_log10[1:]-bins_CAPExP_CSI_log10[0:-1])
bins_CAPExP_CSI_centered = 10**bins_CAPExP_CSI_log10_centered
ax5.loglog(bins_centered,n,'grey')
ax5.loglog(bins_CAPExP_CSI_centered,n_CAPExP_CSI,'k')
ax5.set_title('CAPExP_CSI')
# ax5.set_xlim([common_xmin,1])
# ax5.set_ylim([common_ymin,common_ymax])
ax5.grid(which='major', axis='both', color='lightgray')
ax5.set_ylabel('Frequency')
ax5.set_xlabel('Hourly flash density (-) \n * Hourly LPI (J/kg)')
ax5.annotate('(f) \n ', xy=(x1,  ax5.get_ylim()[1]),annotation_clip=False)

bins_PR92W_log10_centered_d02 = bins_PR92W_log10_d02[0:-1]+0.5*(bins_PR92W_log10_d02[1:]-bins_PR92W_log10_d02[0:-1])
bins_PR92W_centered_d02 = 10**bins_PR92W_log10_centered_d02
ax6.loglog(bins_centered,n,'grey')
ax6.loglog(bins_PR92W_centered_d02,n_PR92W_d02,'k')
ax6.set_title('PR92W')
ax6.grid(which='major', axis='both', color='lightgray')
ax6.set_ylabel('Frequency')
# ax6.set_xlim([common_xmin,1])
# ax6.set_ylim([common_ymin,common_ymax])
ax6.annotate('(d) \n ', xy=(x1, ax6.get_ylim()[1]),annotation_clip=False)

plt.suptitle('Convection-permitting (3 km)')
plt.show()
