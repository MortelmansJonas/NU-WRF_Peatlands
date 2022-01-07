# ---------------------------------------------------------------------------------------------
# MODULES
# ---------------------------------------------------------------------------------------------
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import seaborn as sns
import pandas as pd
# ---------------------------------------------------------------------------------------------
# LOAD DATA
# ---------------------------------------------------------------------------------------------
ds_d01_G = Dataset('/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/wrfout_nc_files/d01_seasonality.nc', 'r')
ds_d02_G = Dataset('/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/wrfout_nc_files/d02_seasonality.nc', 'r')
ds_d01_T = Dataset('/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/wrfout_nc_files/d01_seasonality_Thompson.nc', 'r')
ds_d02_T = Dataset('/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/wrfout_nc_files/d02_seasonality_Thompson.nc', 'r')

Time = pd.date_range('2015-06-01', '2015-08-30', freq='2W')

LPI_T_d01 = ds_d01_T['LPI'][:]
LPI_T_d02 = ds_d02_T['LPI'][:]
LTG3_T_d01 = ds_d01_T['LTG3'][:]
LTG3_T_d02 = ds_d02_T['LTG3'][:]
PR92W_T_d01 = ds_d01_T['PR92W'][:]
PR92W_T_d02 = ds_d02_T['PR92W'][:]
CAPExP_R_T_d01 = ds_d01_T['CAPExP_R'][:]
CAPExP_R_T_d02 = ds_d02_T['CAPExP_R'][:]
CAPExP_CSI_T_d01 = ds_d01_T['CAPExP_CSI'][:]
CAPExP_CSI_T_d02 = ds_d02_T['CAPExP_CSI'][:]
OBS_d01 = ds_d01_T['Obs'][:]
OBS_d02 = ds_d02_T['Obs'][:]
LPI_G_d01 = ds_d01_G['LPI'][:]
LPI_G_d02 = ds_d02_G['LPI'][:]
LTG3_G_d01 = ds_d01_G['LTG3'][:]
LTG3_G_d02 = ds_d02_G['LTG3'][:]
PR92W_G_d01 = ds_d01_G['PR92W'][:]
PR92W_G_d02 = ds_d02_G['PR92W'][:]
CAPExP_R_G_d01 = ds_d01_G['CAPExP_R'][:]
CAPExP_R_G_d02 = ds_d02_G['CAPExP_R'][:]
CAPExP_CSI_G_d01 = ds_d01_G['CAPExP_CSI'][:]
CAPExP_CSI_G_d02 = ds_d02_G['CAPExP_CSI'][:]

# ---------------------------------------------------------------------------------------------
# SEASONALITY
# ---------------------------------------------------------------------------------------------
biweekly_LPI_d01_T = np.zeros((7))
biweekly_LPI_d02_T = np.zeros((7))
biweekly_LTG3_d01_T = np.zeros((7))
biweekly_LTG3_d02_T = np.zeros((7))
biweekly_PR92W_d01_T = np.zeros((7))
biweekly_PR92W_d02_T = np.zeros((7))
biweekly_CAPExP_R_d01_T = np.zeros((7))
biweekly_CAPExP_R_d02_T = np.zeros((7))
biweekly_CAPExP_CSI_d01_T = np.zeros((7))
biweekly_CAPExP_CSI_d02_T = np.zeros((7))
biweekly_Obs_d01 = np.zeros((7))
biweekly_Obs_d02 = np.zeros((7))
biweekly_LPI_d01_G = np.zeros((7))
biweekly_LPI_d02_G = np.zeros((7))
biweekly_LTG3_d01_G = np.zeros((7))
biweekly_LTG3_d02_G = np.zeros((7))
biweekly_PR92W_d01_G = np.zeros((7))
biweekly_PR92W_d02_G = np.zeros((7))
biweekly_CAPExP_R_d01_G = np.zeros((7))
biweekly_CAPExP_R_d02_G = np.zeros((7))
biweekly_CAPExP_CSI_d01_G = np.zeros((7))
biweekly_CAPExP_CSI_d02_G = np.zeros((7))

for i in range(0,7):
    j = (i+1)*2
    k=i*2
    biweekly_LPI_d01_T[i] = np.nanmean(LPI_T_d01[k:j])
    biweekly_LPI_d02_T[i] = np.nanmean(LPI_T_d02[k:j])
    biweekly_LTG3_d01_T[i] = np.nanmean(LTG3_T_d01[k:j])
    biweekly_LTG3_d02_T[i] = np.nanmean(LTG3_T_d02[k:j])
    biweekly_PR92W_d01_T[i] = np.nanmean(PR92W_T_d01[k:j])
    biweekly_PR92W_d02_T[i] = np.nanmean(PR92W_T_d02[k:j])
    biweekly_CAPExP_R_d01_T[i] = np.nanmean(CAPExP_R_T_d01[k:j])
    biweekly_CAPExP_R_d02_T[i] = np.nanmean(CAPExP_R_T_d02[k:j])
    biweekly_CAPExP_CSI_d01_T[i] = np.nanmean(CAPExP_CSI_T_d01[k:j])
    biweekly_CAPExP_CSI_d02_T[i] = np.nanmean(CAPExP_CSI_T_d02[k:j])
    biweekly_Obs_d01[i] = np.nanmean(OBS_d01[k:j])
    biweekly_Obs_d02[i] = np.nanmean(OBS_d02[k:j])
    biweekly_LPI_d01_G[i] = np.nanmean(LPI_G_d01[k:j])
    biweekly_LPI_d02_G[i] = np.nanmean(LPI_G_d02[k:j])
    biweekly_LTG3_d01_G[i] = np.nanmean(LTG3_G_d01[k:j])
    biweekly_LTG3_d02_G[i] = np.nanmean(LTG3_G_d02[k:j])
    biweekly_PR92W_d01_G[i] = np.nanmean(PR92W_G_d01[k:j])
    biweekly_PR92W_d02_G[i] = np.nanmean(PR92W_G_d02[k:j])
    biweekly_CAPExP_R_d01_G[i] = np.nanmean(CAPExP_R_G_d01[k:j])
    biweekly_CAPExP_R_d02_G[i] = np.nanmean(CAPExP_R_G_d02[k:j])
    biweekly_CAPExP_CSI_d01_G[i] = np.nanmean(CAPExP_CSI_G_d01[k:j])
    biweekly_CAPExP_CSI_d02_G[i] = np.nanmean(CAPExP_CSI_G_d02[k:j])
LPI_T_d01 = biweekly_LPI_d01_T
LPI_T_d02 = biweekly_LPI_d02_T
LTG3_T_d01 = biweekly_LTG3_d01_T
LTG3_T_d02 = biweekly_LTG3_d02_T
PR92W_T_d01 = biweekly_PR92W_d01_T
PR92W_T_d02 = biweekly_PR92W_d02_T
CAPExP_R_T_d01 = biweekly_CAPExP_R_d01_T
CAPExP_R_T_d02 = biweekly_CAPExP_R_d02_T
CAPExP_CSI_T_d01 = biweekly_CAPExP_CSI_d01_T
CAPExP_CSI_T_d02 = biweekly_CAPExP_CSI_d01_T
OBS_d01 = biweekly_Obs_d01
OBS_d02 = biweekly_Obs_d02
LPI_G_d01 = biweekly_LPI_d01_G
LPI_G_d02 = biweekly_LPI_d02_G
LTG3_G_d01 = biweekly_LTG3_d01_G
LTG3_G_d02 = biweekly_LTG3_d02_G
PR92W_G_d01 = biweekly_PR92W_d01_G
PR92W_G_d02 = biweekly_PR92W_d02_G
CAPExP_R_G_d01 = biweekly_CAPExP_R_d01_G
CAPExP_R_G_d02 = biweekly_CAPExP_R_d02_G
CAPExP_CSI_G_d01 = biweekly_CAPExP_CSI_d01_G
CAPExP_CSI_G_d02 = biweekly_CAPExP_CSI_d01_G

# ---------------------------------------------------------------------------------------------
# PLOTS
# ---------------------------------------------------------------------------------------------
palette = sns.color_palette('colorblind')
fig, axes = plt.subplots(2,2, figsize=(15.3, 8.27))
fig.suptitle('Seasonality', fontsize = 32)

axes[0][1].plot(Time,OBS_d01, color='lightgrey', linewidth = 6)
axes[0][1].plot(Time,LPI_T_d01,marker='d', ms=7, color=palette[0])
axes[0][1].plot(Time,LTG3_T_d01,marker='s', ms=7, color=palette[1])
axes[0][1].plot(Time,PR92W_T_d01,marker='v', ms=7, color=palette[2])
axes[0][1].plot(Time,CAPExP_R_T_d01,marker='X', ms=7, color=palette[3])
axes[0][1].plot(Time,CAPExP_CSI_T_d01,marker='P', ms=7, color=palette[4])
axes[0][1].grid(which='major', axis='both', color='lightgray')
axes[0][1].set_title('THOM - 9 km', fontsize=28)
axes[0][1].ticklabel_format(axis='y',style='sci', scilimits=(0,0))
axes[0][1].set_ylim(0,0.02)
axes[0][1].tick_params(axis='both', which='major',labelsize=16)
axes[0][1].set_xticks(ticks=[pd.to_datetime('2015-06', format='%Y-%m'), pd.to_datetime('2015-07', format='%Y-%m'),
                     pd.to_datetime('2015-08', format='%Y-%m'), pd.to_datetime('2015-09', format='%Y-%m')])
axes[0][1].set_xticklabels(['JUN', 'JUL', 'AUG', 'SEP'])
axes[0][1].set_yticks(ticks=[0, 0.005, 0.01, 0.015, 0.02])
axes[0][1].set_yticklabels(['0.0', '0.5', '1.0', '1.5', '2.0'])
axes[0][1].annotate('(1e-2)', xy=(0.005,  axes[0][1].get_ylim()[1]),annotation_clip=False, fontsize=16)

axes[1][1].plot(Time,OBS_d02,  label='CLDN', color='lightgrey', linewidth = 6)
axes[1][1].plot(Time,LPI_T_d02,marker='d', ms=7, label='LPI', color=palette[0])
axes[1][1].plot(Time,LTG3_T_d02,marker='s', ms=7, label='LT3', color=palette[1])
axes[1][1].plot(Time,PR92W_T_d02,marker='v', ms=7, label='PR92W', color=palette[2])
axes[1][1].plot(Time,CAPExP_R_T_d02,marker='X', ms=7, label='CAPEXP_R', color=palette[3])
axes[1][1].plot(Time,CAPExP_CSI_T_d02,marker='P', ms=7, label='CAPEXP_CSI', color=palette[4])
axes[1][1].grid(which='major', axis='both', color='lightgray')
axes[1][1].legend(bbox_to_anchor=(0.7, -0.3), ncol=6,fontsize=16)
axes[1][1].ticklabel_format(axis='y',style='sci', scilimits=(0,0))
axes[1][1].set_ylim(0,0.02)
axes[1][1].tick_params(axis='both', which='major',labelsize=16)
axes[1][1].set_xticks(ticks=[pd.to_datetime('2015-06', format='%Y-%m'), pd.to_datetime('2015-07', format='%Y-%m'),
                     pd.to_datetime('2015-08', format='%Y-%m'), pd.to_datetime('2015-09', format='%Y-%m')])
axes[1][1].set_xticklabels(['JUN', 'JUL', 'AUG', 'SEP'])
axes[1][1].set_yticks(ticks=[0, 0.005, 0.01, 0.015, 0.02])
axes[1][1].set_yticklabels(['0.0', '0.5', '1.0', '1.5', '2.0'])
axes[1][1].annotate('(1e-2)', xy=(0.0005,  axes[1][1].get_ylim()[1]),annotation_clip=False, fontsize=16)
axes[1][1].set_title('THOM - 3 km', fontsize=28)

axes[0][0].plot(Time,OBS_d01,color='lightgrey', linewidth = 6)
axes[0][0].plot(Time,LPI_G_d01,marker='d', ms=7, color=palette[0])
axes[0][0].plot(Time,LTG3_G_d01,marker='s', ms=7, color=palette[1])
axes[0][0].plot(Time,PR92W_G_d01,marker='v', ms=7, color=palette[2])
axes[0][0].plot(Time,CAPExP_R_G_d01,marker='X', ms=7, color=palette[3])
axes[0][0].plot(Time,CAPExP_CSI_G_d01,marker='P', ms=7, color=palette[4])
axes[0][0].grid(which='major', axis='both', color='lightgray')
axes[0][0].set_title('G4ICE - 9 km', fontsize=28)
axes[0][0].set_ylabel('Flash density \n (# day$^{-1}$ km$^{-2}$)', fontsize=24)
axes[0][0].ticklabel_format(axis='y',style='sci', scilimits=(0,0))
axes[0][0].set_ylim(0,0.02)
axes[0][0].tick_params(axis='both', which='major',labelsize=16)
axes[0][0].set_xticks(ticks=[pd.to_datetime('2015-06', format='%Y-%m'), pd.to_datetime('2015-07', format='%Y-%m'),
                     pd.to_datetime('2015-08', format='%Y-%m'), pd.to_datetime('2015-09', format='%Y-%m')])
axes[0][0].set_xticklabels(['JUN', 'JUL', 'AUG', 'SEP'])
axes[0][0].set_yticks(ticks=[0, 0.005, 0.01, 0.015, 0.02])
axes[0][0].set_yticklabels(['0.0', '0.5', '1.0', '1.5', '2.0'])
axes[0][0].annotate('(1e-2)', xy=(0.0005,  axes[0][0].get_ylim()[1]),annotation_clip=False, fontsize=16)

axes[1][0].plot(Time,OBS_d02, color='lightgrey', linewidth = 6)
axes[1][0].plot(Time,LPI_G_d02,marker='d', ms=7, color=palette[0])
axes[1][0].plot(Time,LTG3_G_d02,marker='s', ms=7, color=palette[1])
axes[1][0].plot(Time,PR92W_G_d02,marker='v', ms=7, color=palette[2])
axes[1][0].plot(Time,CAPExP_R_G_d02,marker='X', ms=7, color=palette[3])
axes[1][0].plot(Time,CAPExP_CSI_G_d02,marker='P', ms=7, color=palette[4])
axes[1][0].grid(which='major', axis='both', color='lightgray')
axes[1][0].ticklabel_format(axis='y',style='sci', scilimits=(0,0))
axes[1][0].set_ylim(0,0.02)
axes[1][0].tick_params(axis='both', which='major',labelsize=16)
axes[1][0].set_ylabel('Flash density\n(# day$^{-1}$ km$^{-2}$)', fontsize=24)
axes[1][0].set_title('G4ICE - 3 km', fontsize=28)
axes[1][0].set_xticks(ticks=[pd.to_datetime('2015-06', format='%Y-%m'), pd.to_datetime('2015-07', format='%Y-%m'),
                     pd.to_datetime('2015-08', format='%Y-%m'), pd.to_datetime('2015-09', format='%Y-%m')])
axes[1][0].set_xticklabels(['JUN', 'JUL', 'AUG', 'SEP'])
axes[1][0].set_yticks(ticks=[0, 0.005, 0.01, 0.015, 0.02])
axes[1][0].set_yticklabels(['0.0', '0.5', '1.0', '1.5', '2.0'])
axes[1][0].annotate('(1e-2)\n', xy=(0.0005,  axes[1][0].get_ylim()[1]),annotation_clip=False, fontsize=16)
plt.show()
