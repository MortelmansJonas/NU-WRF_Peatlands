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
ds_d01_T = Dataset('/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/NU-WRF/wrfout_nc_files/d01_diurnal_Thompson.nc', 'r')
ds_d02_T = Dataset('/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/NU-WRF/wrfout_nc_files/d02_diurnal_Thompson.nc', 'r')
ds_d01_G = Dataset('/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/NU-WRF/wrfout_nc_files/d01_diurnal.nc', 'r')
ds_d02_G = Dataset('/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/NU-WRF/wrfout_nc_files/d02_diurnal.nc', 'r')

Time = np.linspace(0,24, num=24, endpoint=False)

LPI_T_d01 = ds_d01_T['LPI'][:]
LPI_T_d02 = ds_d02_T['LPI'][:]
LTG3_T_d01 = ds_d01_T['LTG3'][:]
LTG3_T_d02 = ds_d02_T['LTG3'][:]
PR92W_T_d01 = ds_d01_T['PR92W'][:]
PR92W_T_d02 = ds_d02_T['PR92W'][:]
CAPExP_R_T_d01 = ds_d01_T['CAPExP_R'][:]
CAPExP_R_T_d02 = ds_d02_T['CAPExP_R'][:]
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

# ---------------------------------------------------------------------------------------------
# PLOT D01
# ---------------------------------------------------------------------------------------------
palette = sns.color_palette('colorblind')
fig, axes = plt.subplots(2,2, figsize=(15.3, 8.27), dpi=150)
fig.suptitle('Diurnal Cycle', fontsize = 20)

axes[0][1].plot(Time,OBS_d01, color='lightgrey', linewidth = 6)
axes[0][1].plot(Time,LPI_T_d01,marker='d', ms=5, color=palette[0])
axes[0][1].plot(Time,LTG3_T_d01,marker='s', ms=5, color=palette[1])
axes[0][1].plot(Time,PR92W_T_d01,marker='v', ms=5, color=palette[4])
axes[0][1].plot(Time,CAPExP_R_T_d01,marker='X', ms=5, color=palette[3])
axes[0][1].grid(which='major', axis='both', color='lightgray')
axes[0][1].set_title('THOM - 9 km', fontsize=14)
axes[0][1].ticklabel_format(axis='y',style='sci', scilimits=(0,0))
# axes[0][1].set_ylim(0,0.02)
axes[0][1].set_ylim(0,0.0014)
axes[0][1].tick_params(axis='both', which='major', labelsize=12)
axes[0][1].set_yticks(ticks=[0, 0.0002, 0.0004, 0.0006, 0.0008, 0.001, 0.0012, 0.0014])

axes[1][1].plot(Time,OBS_d02,  label='CLDN', color='lightgrey', linewidth = 6)
axes[1][1].plot(Time,LPI_T_d02,marker='d', ms=5, label='LPI', color=palette[0])
axes[1][1].plot(Time,LTG3_T_d02,marker='s', ms=5, label='LT', color=palette[1])
axes[1][1].plot(Time,PR92W_T_d02,marker='v', ms=5, label='PR92W', color=palette[4])
axes[1][1].plot(Time,CAPExP_R_T_d02,marker='X', ms=5, label='CAPExP', color=palette[3])
axes[1][1].grid(which='major', axis='both', color='lightgray')
axes[1][1].legend(bbox_to_anchor=(0.5, -0.3), ncol=5, fontsize=12)
axes[1][1].ticklabel_format(axis='y',style='sci', scilimits=(0,0))
axes[1][1].set_ylim(0,0.0014)
axes[1][1].tick_params(axis='both', which='major', labelsize=12)
axes[1][1].set_xlabel('Hour of the day (local time)', fontsize=12)
axes[1][1].set_yticks(ticks=[0, 0.0002, 0.0004, 0.0006, 0.0008, 0.001, 0.0012, 0.0014])
axes[1][1].set_title('THOM - 3 km', fontsize=14)

axes[0][0].plot(Time,OBS_d01,color='lightgrey', linewidth = 6)
axes[0][0].plot(Time,LPI_G_d01,marker='d', ms=5, color=palette[0])
axes[0][0].plot(Time,LTG3_G_d01,marker='s', ms=5, color=palette[1])
axes[0][0].plot(Time,PR92W_G_d01,marker='v', ms=5, color=palette[4])
axes[0][0].plot(Time,CAPExP_R_G_d01,marker='X', ms=5, color=palette[3])
axes[0][0].grid(which='major', axis='both', color='lightgray')
axes[0][0].set_title('G4ICE - 9 km', fontsize=14)
axes[0][0].set_ylabel('Flash density \n (# hr$^{-1}$ km$^{-2}$)', fontsize=12)
axes[0][0].ticklabel_format(axis='y',style='sci', scilimits=(0,0))
axes[0][0].set_ylim(0,0.0014)
axes[0][0].tick_params(axis='both', which='major', labelsize=12)
axes[0][0].set_yticks(ticks=[0, 0.0002, 0.0004, 0.0006, 0.0008, 0.001, 0.0012, 0.0014])

axes[1][0].plot(Time,OBS_d02, color='lightgrey', linewidth = 6)
axes[1][0].plot(Time,LPI_G_d02,marker='d', ms=5, color=palette[0])
axes[1][0].plot(Time,LTG3_G_d02,marker='s', ms=5, color=palette[1])
axes[1][0].plot(Time,PR92W_G_d02,marker='v', ms=5, color=palette[4])
axes[1][0].plot(Time,CAPExP_R_G_d02,marker='X', ms=5, color=palette[3])
axes[1][0].grid(which='major', axis='both', color='lightgray')
axes[1][0].ticklabel_format(axis='y',style='sci', scilimits=(0,0))
axes[1][0].set_ylim(0,0.0014)
axes[1][0].tick_params(axis='both', which='major', labelsize=12)
axes[1][0].set_ylabel('Flash density\n(# hr$^{-1}$ km$^{-2}$)', fontsize=12)
axes[1][0].set_title('G4ICE - 3 km', fontsize=14)
axes[1][0].set_xlabel('Hour of the day (local time)', fontsize=12)
axes[1][0].set_yticks(ticks=[0, 0.0002, 0.0004, 0.0006, 0.0008, 0.001, 0.0012, 0.0014])
plt.show()
