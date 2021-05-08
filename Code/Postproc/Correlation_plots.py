## Load Modules
import numpy as np
import pandas as pd
from netCDF4 import Dataset
from matplotlib import pyplot as plt
import seaborn as sns

# adjust directory
outpath = '/scratch/leuven/336/vsc33651/nu-wrf-dev/Plots'
figname = 'correlation_for_temporal_mean.png'
figname_diurnal = 'correlation_for_diurnal.png'
figname_seasonal = 'correlation_for_seasonality.png'

# d01 model
infile = '/scratch/leuven/336/vsc33651/nu-wrf-dev/domain1_all_v4.nc'
ds01 = Dataset(infile, 'r')
ds01_diurnal = Dataset('/scratch/leuven/336/vsc33651/nu-wrf-dev/domain1_diurnal.nc', 'r')
ds01_seasonal = Dataset('/scratch/leuven/336/vsc33651/nu-wrf-dev/domain1_seasonality.nc', 'r')

# d02 model
infile = '/scratch/leuven/336/vsc33651/nu-wrf-dev/domain2_at_domain1_all_v4.nc'
ds02 = Dataset(infile, 'r')
ds02_diurnal = Dataset('/scratch/leuven/336/vsc33651/nu-wrf-dev/domain2_diurnal.nc', 'r')
ds02_seasonal = Dataset('/scratch/leuven/336/vsc33651/nu-wrf-dev/domain2_seasonality.nc', 'r')

# obs CLDN
obsfile = '/scratch/leuven/336/vsc33651/nu-wrf-dev/CLDN_at_domain1_all_v4.nc'
dsobs = Dataset(obsfile, 'r')

# calculate temporal mean
obs = dsobs['Flashdensity_CC'][:].data + dsobs['Flashdensity_CG'][:].data
obs = np.nanmean(obs,0)
# temporal mean for ds01
CAPExP = ds01['CAPExP'][:].data
CAPExP_ds01 = np.nanmean(CAPExP,0)
LPI = ds01['LPI'][:].data
LPI_ds01 = np.nanmean(LPI,0)
LT3 = ds01['LTG3'][:].data
LT3_ds01 = np.nanmean(LT3,0)
PR92W = ds01['PR92W'][:].data
PR92W_ds01 = np.nanmean(PR92W,0)
# temporal mean for ds02
CAPExP = ds02['CAPExP'][:].data
CAPExP_ds02 = np.nanmean(CAPExP,0)
LPI = ds02['LPI'][:].data
LPI_ds02 = np.nanmean(LPI,0)
LT3 = ds02['LTG3'][:].data
LT3_ds02 = np.nanmean(LT3,0)
PR92W = ds02['PR92W'][:].data
PR92W_ds02 = np.nanmean(PR92W,0)
# Calculate spatial mean for diurnal cycle d01
obs_ds01_diurnal = ds01_diurnal['Obs'][:].data
CAPExP_ds01_diurnal = ds01_diurnal['CAPExP'][:].data
LPI_ds01_diurnal = ds01_diurnal['LPI'][:].data
LT3_ds01_diurnal = ds01_diurnal['LTG3'][:].data
PR92W_ds01_diurnal = ds01_diurnal['PR92W'][:].data
# Calculate spatial mean for diurnal cycle d02
obs_ds02_diurnal = ds02_diurnal['Obs'][:].data
LPI_ds02_diurnal = ds02_diurnal['LPI'][:].data
LT3_ds02_diurnal = ds02_diurnal['LTG3'][:].data
PR92W_ds02_diurnal = ds02_diurnal['PR92W'][:].data
# Calculate spatial mean for seasonality d01
obs_ds01_seasonal = ds01_seasonal['Obs'][:].data
CAPExP_ds01_seasonal = ds01_seasonal['CAPExP'][:].data
LPI_ds01_seasonal = ds01_seasonal['LPI'][:].data
LT3_ds01_seasonal = ds01_seasonal['LTG3'][:].data
PR92W_ds01_seasonal = ds01_seasonal['PR92W'][:].data
# Calculate spatial mean for seasonality d02
obs_ds02_seasonal = ds02_seasonal['Obs'][:].data
LPI_ds02_seasonal = ds02_seasonal['LPI'][:].data
LT3_ds02_seasonal = ds02_seasonal['LTG3'][:].data
PR92W_ds02_seasonal = ds02_seasonal['PR92W'][:].data

# prepare fields for panda correlation function, flatten 2D data.
d = {'CLDN': obs.flatten(),\
     'CAPExP 9km': CAPExP_ds01.flatten(), \
     'LPI 9km': LPI_ds01.flatten(), \
     'LT3 9km': LT3_ds01.flatten(), \
     'PR92W 9km': PR92W_ds01.flatten(), \
     'LPI 3km': LPI_ds02.flatten(), \
     'LT3 3km': LT3_ds02.flatten(), \
     'PR92W 3km': PR92W_ds02.flatten()\
     }
d_diurnal = {'CLDN': obs_ds01_diurnal.flatten(),\
     'CAPExP 9km': CAPExP_ds01_diurnal.flatten(), \
     'LPI 9km': LPI_ds01_diurnal.flatten(), \
     'LT3 9km': LT3_ds01_diurnal.flatten(), \
     'PR92W 9km': PR92W_ds01_diurnal.flatten(), \
     'LPI 3km': LPI_ds02_diurnal.flatten(), \
     'LT3 3km': LT3_ds02_diurnal.flatten(), \
     'PR92W 3km': PR92W_ds02_diurnal.flatten()\
     }
d_seasonal = {'CLDN': obs_ds01_seasonal,\
     'CAPExP 9km': CAPExP_ds01_seasonal, \
     'LPI 9km': LPI_ds01_seasonal, \
     'LT3 9km': LT3_ds01_seasonal, \
     'PR92W 9km': PR92W_ds01_seasonal, \
     'LPI 3km': LPI_ds02_seasonal, \
     'LT3 3km': LT3_ds02_seasonal, \
     'PR92W 3km': PR92W_ds02_seasonal\
     }
df = pd.DataFrame(data=d)
df_diurnal = pd.DataFrame(data=d_diurnal)
df_seasonal = pd.DataFrame(data=d_seasonal)
# panda dataframe is ready, apply correlation matrix function,
# default method is pearson
#corr_matrix = df.corr(method='spearman')
# Pearson:
corr_matrix = df.corr()
corr_matrix_diurnal = df_diurnal.corr()
corr_matrix_seasonal = df_seasonal.corr()
# hide upper triangle
mask = np.zeros_like(corr_matrix, dtype=np.bool)
mask[np.triu_indices_from(mask)]= True
mask_diurnal = np.zeros_like(corr_matrix_diurnal, dtype=np.bool)
mask_diurnal[np.triu_indices_from(mask_diurnal)]= True
mask_seasonal = np.zeros_like(corr_matrix_seasonal, dtype=np.bool)
mask_seasonal[np.triu_indices_from(mask_seasonal)]= True

#f, ax = plt.subplots(figsize=(11, 15))
f, ax = plt.subplots(figsize=(11*0.8, 15*0.8))

heatmap = sns.heatmap(corr_matrix, mask = mask, square = True, linewidths = .5, cmap = 'coolwarm',\
cbar_kws = {'shrink': .4, 'ticks' : [0, 0.5, 1]},vmin = 0,vmax = 1,annot = True,annot_kws = {'size': 12})

#add the column names as labels
ax.set_yticklabels(corr_matrix.columns, rotation = 0)
ax.set_xticklabels(corr_matrix.columns)
ax.set_title('Spatial correlation (-) of long-term average')

sns.set_style({'xtick.bottom': True}, {'ytick.left': True})
heatmap.get_figure().savefig(outpath+figname, bbox_inches='tight')

# Diurnal
f_diurnal, ax_diurnal = plt.subplots(figsize=(11*0.8, 15*0.8))

heatmap_diurnal = sns.heatmap(corr_matrix_diurnal, mask = mask_diurnal, square = True, linewidths = .5, cmap = 'coolwarm',\
cbar_kws = {'shrink': .4, 'ticks' : [0, 0.5, 1]},vmin = 0,vmax = 1,annot = True,annot_kws = {'size': 12})

#add the column names as labels
ax_diurnal.set_yticklabels(corr_matrix_diurnal.columns, rotation = 0)
ax_diurnal.set_xticklabels(corr_matrix_diurnal.columns)
ax_diurnal.set_title('Temporal correlation (-) of diurnal cycle')

sns.set_style({'xtick.bottom': True}, {'ytick.left': True})
heatmap_diurnal.get_figure().savefig(outpath+figname_diurnal, bbox_inches='tight')

# Seasonal
f_seasonal, ax_seasonal = plt.subplots(figsize=(11*0.8, 15*0.8))

heatmap_seasonal = sns.heatmap(corr_matrix_seasonal, mask = mask_seasonal, square = True, linewidths = .5, cmap = 'coolwarm',\
cbar_kws = {'shrink': .4, 'ticks' : [0, 0.5, 1]},vmin = 0,vmax = 1,annot = True,annot_kws = {'size': 12})

#add the column names as labels
ax_seasonal.set_yticklabels(corr_matrix_seasonal.columns, rotation = 0)
ax_seasonal.set_xticklabels(corr_matrix_seasonal.columns)
ax_seasonal.set_title('Temporal correlation (-) of seasonality')


sns.set_style({'xtick.bottom': True}, {'ytick.left': True})
heatmap_seasonal.get_figure().savefig(outpath+figname_seasonal, bbox_inches='tight')
