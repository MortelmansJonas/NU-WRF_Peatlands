# ---------------------------------------------------------------------------------------------
# MODULES
# ---------------------------------------------------------------------------------------------
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import seaborn as sns

# ---------------------------------------------------------------------------------------------
# LOAD DATA
# ---------------------------------------------------------------------------------------------
ds_d01_G = Dataset('/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/wrfout_nc_files/data_calibrated_ax_d01.nc', 'r')
ds_d01_T = Dataset('/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/wrfout_nc_files/data_calibrated_ax_d01_Thompson.nc', 'r')

lat_d01 = ds_d01_G['lat'][:]
lon_d01 = ds_d01_G['lon'][:]

# ---------------------------------------------------------------------------------------------
# CALCULATE TEMPORAL MEAN
# ---------------------------------------------------------------------------------------------
obs = np.nanmean(ds_d01_T['Obs'][:],0)
CAPExP_R_d01_T = np.nanmean(ds_d01_T['CAPExP_R_d01'][:],0)
CAPExP_CSI_d01_T = np.nanmean(ds_d01_T['CAPExP_CSI_d01'][:],0)
LPI_d01_T = np.nanmean(ds_d01_T['LPI_d01'][:],0)
LT3_d01_T = np.nanmean(ds_d01_T['LTG3_d01'][:],0)
PR92W_d01_T = np.nanmean(ds_d01_T['PR92W_d01'][:],0)

CAPExP_R_d02_T = np.nanmean(ds_d01_T['CAPExP_R_d02'][:],0)
CAPExP_CSI_d02_T = np.nanmean(ds_d01_T['CAPExP_CSI_d02'][:],0)
LPI_d02_T = np.nanmean(ds_d01_T['LPI_d02'][:],0)
LT3_d02_T = np.nanmean(ds_d01_T['LTG3_d02'][:],0)
PR92W_d02_T = np.nanmean(ds_d01_T['PR92W_d02'][:],0)

CAPExP_R_d01_G = np.nanmean(ds_d01_G['CAPExP_R_d01'][:],0)
CAPExP_CSI_d01_G = np.nanmean(ds_d01_G['CAPExP_CSI_d01'][:],0)
LPI_d01_G = np.nanmean(ds_d01_G['LPI_d01'][:],0)
LT3_d01_G = np.nanmean(ds_d01_G['LTG3_d01'][:],0)
PR92W_d01_G = np.nanmean(ds_d01_G['PR92W_d01'][:],0)

CAPExP_R_d02_G = np.nanmean(ds_d01_G['CAPExP_R_d02'][:],0)
CAPExP_CSI_d02_G = np.nanmean(ds_d01_G['CAPExP_CSI_d02'][:],0)
LPI_d02_G = np.nanmean(ds_d01_G['LPI_d02'][:],0)
LT3_d02_G = np.nanmean(ds_d01_G['LTG3_d02'][:],0)
PR92W_d02_G = np.nanmean(ds_d01_G['PR92W_d02'][:],0)

# prepare fields for panda correlation function, flatten 2D data.
d_d01_T = {'CLDN': obs[:].flatten(),\
           'LPI': LPI_d01_T[:].flatten(), \
           'LT3': LT3_d01_T[:].flatten(),\
           'PR92W': PR92W_d01_T[:].flatten(), \
           'CAPExP_R': CAPExP_R_d01_T[:].flatten(), \
           'CAPExP_CSI': CAPExP_CSI_d01_T[:].flatten()}
d_d02_T = {'CLDN': obs[:].flatten(),\
           'LPI': LPI_d02_T[:].flatten(), \
           'LT3': LT3_d02_T[:].flatten(),\
           'PR92W': PR92W_d02_T[:].flatten(), \
           'CAPExP_R': CAPExP_R_d02_T[:].flatten(), \
           'CAPExP_CSI': CAPExP_CSI_d02_T[:].flatten()}
d_d01_G = {'CLDN': obs[:].flatten(),\
           'LPI': LPI_d01_G[:].flatten(), \
           'LT3': LT3_d01_G[:].flatten(),\
           'PR92W': PR92W_d01_G[:].flatten(), \
           'CAPExP_R': CAPExP_R_d01_G[:].flatten(), \
           'CAPExP_CSI': CAPExP_CSI_d01_G[:].flatten()}
d_d02_G = {'CLDN': obs[:].flatten(),\
           'LPI': LPI_d02_G[:].flatten(), \
           'LT3': LT3_d02_G[:].flatten(),\
           'PR92W': PR92W_d02_G[:].flatten(), \
           'CAPExP_R': CAPExP_R_d02_G[:].flatten(), \
           'CAPExP_CSI': CAPExP_CSI_d02_G[:].flatten()}

df_d01_T = pd.DataFrame(data=d_d01_T)
df_d02_T = pd.DataFrame(data=d_d02_T)
df_d01_G = pd.DataFrame(data=d_d01_G)
df_d02_G = pd.DataFrame(data=d_d02_G)
# Pearson:
corr_matrix_d01_T = df_d01_T.corr()
corr_matrix_d02_T = df_d02_T.corr()
corr_matrix_d01_G = df_d01_G.corr()
corr_matrix_d02_G = df_d02_G.corr()

corr_table = np.zeros((4,5))
corr_table[0,:] = corr_matrix_d01_T.iloc[1:,0]
corr_table[1,:] = corr_matrix_d02_T.iloc[1:,0]
corr_table[2,:] = corr_matrix_d01_G.iloc[1:,0]
corr_table[3,:] = corr_matrix_d02_G.iloc[1:,0]

print(corr_table)
corr_table_df = pd.DataFrame(corr_table, index=['THOM - 9 km', 'THOM - 3 km', 'G4ICE - 9 km', 'G4ICE - 3 km'],
                             columns=['LPI', 'LT3', 'PR92W', 'CAPExP_R', 'CAPExP_CSI'])

f, ax = plt.subplots(figsize=(11, 15))
heatmap = sns.heatmap(np.round(corr_table_df,4), square = True, linewidths = .5, cmap = 'Blues',\
cbar_kws = {'shrink': .4, 'ticks' : [0, 0.5, 1]},vmin = 0,vmax = 1,annot = True,annot_kws = {'size': 20})

#add the column names as labels
ax.set_yticklabels(corr_table_df.index, rotation = 0, fontsize=20)
ax.set_xticklabels(corr_table_df.columns, fontsize=20, rotation=30)
ax.set_title('Spatial correlation (-) of long-term average', fontsize=28)

cbar = ax.collections[0].colorbar
# here set the labelsize by 14
cbar.ax.tick_params(labelsize=20)

sns.set_style({'xtick.bottom': True}, {'ytick.left': True})
plt.show()
