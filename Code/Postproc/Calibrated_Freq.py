# ---------------------------------------------------------------------------------------------
# MODULES
# ---------------------------------------------------------------------------------------------
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
# ---------------------------------------------------------------------------------------------
# LOAD DATA
# ---------------------------------------------------------------------------------------------
ds = Dataset('/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/wrfout_nc_files/data_calibrated_ax.nc','r')

# ---------------------------------------------------------------------------------------------
# NORMALIZE DATA
# ---------------------------------------------------------------------------------------------
Ln = ds['Obs'][:]
LPI_n_d01 = ds['LPI_d01'][:]
LPI_n_d02 = ds['LPI_d02'][:]
LTG3_n_d01 = ds['LTG3_d01'][:]
LTG3_n_d02 = ds['LTG3_d02'][:]
PR92W_n_d01 = ds['PR92W_d01'][:]
PR92W_n_d02 = ds['PR92W_d02'][:]
CAPExP_R_n_d01 = ds['CAPExP_R_d01'][:]
CAPExP_R_n_d02 = ds['CAPExP_R_d02'][:]
CAPExP_CSI_n_d01 = ds['CAPExP_CSI_d01'][:]
CAPExP_CSI_n_d02 = ds['CAPExP_CSI_d02'][:]

# ---------------------------------------------------------------------------------------------
# exclude low values that likely don't produce lightning
# Remove lightning predictions with less than 1 flash per hour per storm (e.g. 100 sqkm)
cond_LTG3_d01 = ds['LTG3_d01'][:] < 0.01
cond_LTG3_d02 = ds['LTG3_d02'][:] < 0.01
cond_PR92W_d01 = ds['PR92W_d01'][:] < 0.01
cond_PR92W_d02 = ds['PR92W_d02'][:] < 0.01
cond_CAPExP_CSI_d01 = ds['CAPExP_CSI_d01'][:] < 0.01
cond_CAPExP_CSI_d02 = ds['CAPExP_CSI_d02'][:] < 0.01
cond_CAPExP_R_d01 = ds['CAPExP_R_d01'][:] < 0.01
cond_CAPExP_R_d02 = ds['CAPExP_R_d02'][:] < 0.01
# for LPI everything smaller than 0.001
cond_LPI_d01 = LPI_n_d01<0.001
cond_LPI_d02 = LPI_n_d02<0.001
# set to nan
LPI_n_d01[cond_LPI_d01] = np.nan
LPI_n_d02[cond_LPI_d02] = np.nan
LTG3_n_d01[cond_LTG3_d01] = np.nan
LTG3_n_d02[cond_LTG3_d02] = np.nan
PR92W_n_d01[cond_PR92W_d01] = np.nan
PR92W_n_d02[cond_PR92W_d02] = np.nan
CAPExP_CSI_n_d01[cond_CAPExP_CSI_d01] = np.nan
CAPExP_CSI_n_d02[cond_CAPExP_CSI_d02] = np.nan
CAPExP_R_n_d01[cond_CAPExP_R_d01] = np.nan
CAPExP_R_n_d02[cond_CAPExP_R_d02] = np.nan

# ---------------------------------------------------------------------------------------------
# FREQUENCY IFO HOURLY GRID FLASH DENSITY
# ---------------------------------------------------------------------------------------------
# specify bin edges
my_bins = np.linspace(-3,2,20)
# DOMAIN 1
n, bins_log10, patches = plt.hist(np.log10(Ln[Ln!=0]),bins = my_bins)
plt.close()
n_CAPExP_R, bins_CAPExP_R_log10, patches_CAPExP_R = plt.hist(np.log10(CAPExP_R_n_d01[CAPExP_R_n_d01!=0]),bins=my_bins)
plt.close()
n_CAPExP_CSI, bins_CAPExP_CSI_log10, patches_CAPExP_CSI = plt.hist(np.log10(CAPExP_CSI_n_d01[CAPExP_CSI_n_d01!=0]),bins=my_bins)
plt.close()
n_lpi, bins_lpi_log10, patches_lpi = plt.hist(np.log10(LPI_n_d01[LPI_n_d01!=0]),bins=my_bins)
plt.close()
n_LTG3, bins_LTG3_log10, patches_LTG3 = plt.hist(np.log10(LTG3_n_d01[LTG3_n_d01!=0]),bins=my_bins)
plt.close()
n_PR92W, bins_PR92W_log10, patches_PR92W = plt.hist(np.log10(PR92W_n_d01[PR92W_n_d01!=0]),bins=my_bins)
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
n, bins_log10, patches = plt.hist(np.log10(Ln[Ln!=0]),bins=my_bins)
plt.close()
n_CAPExP_R, bins_CAPExP_R_log10, patches_CAPExP_R = plt.hist(np.log10(CAPExP_R_n_d02[CAPExP_R_n_d02!=0]),bins=my_bins)
plt.close()
n_CAPExP_CSI, bins_CAPExP_CSI_log10, patches_CAPExP_CSI = plt.hist(np.log10(CAPExP_CSI_n_d02[CAPExP_CSI_n_d02!=0]),bins=my_bins)
plt.close()
n_lpi_d02, bins_lpi_log10_d02, patches_lpi_d02 = plt.hist(np.log10(LPI_n_d02[LPI_n_d02!=0]),bins=my_bins)
plt.close()
n_LTG3_d02, bins_LTG3_log10_d02, patches_LTG3_d02 = plt.hist(np.log10(LTG3_n_d02[LTG3_n_d02!=0]),bins=my_bins)
plt.close()
n_PR92W_d02, bins_PR92W_log10_d02, patches_PR92W_d02 = plt.hist(np.log10(PR92W_n_d02[PR92W_n_d02!=0]),bins=my_bins)
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
