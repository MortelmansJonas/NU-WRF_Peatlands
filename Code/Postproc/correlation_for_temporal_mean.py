## Load Modules
import numpy as np
import pandas as pd
import os
from netCDF4 import Dataset
from wrf import getvar, latlon_coords
import datetime as dt
from pyproj import Proj, transform
#cmd = os.popen('module load CDO/1.9.7.1-intel-2018a')
#cmd = os.popen('/usr/bin/modulecmd python load CDO/1.9.7.1-intel-2018a')
#exec(cmd)
#import python as mod
#mod.module('load', 'program/1.2.3')
import sys
#sys.path.append('/usr/share/Modules/init')
#exec(open('/usr/share/Modules/init/python.py').read())
#import nctoolkit as nc
from pyresample.geometry import SwathDefinition
from pyresample.kd_tree import resample_nearest
from pyresample import bilinear, geometry
from matplotlib import pyplot as plt
import netCDF4 as nc
import xarray as xr
from scipy import ndimage
from dateutil.relativedelta import *
import seaborn as sns

# adjust directory
outpath = '/staging/leuven/stg_00024/OUTPUT/michelb/FIG_tmp/'
figname = 'correlation_for_temporal_mean.png'

# d01 model
infile = '/scratch/leuven/336/vsc33651/nu-wrf-dev/wrfout_nc_files/domain1_all.nc'
ds01 = Dataset(infile, 'r')
# d02 model
infile = '/scratch/leuven/317/vsc31786/nu-wrf-dev/wrfout_nc_files/domain2_at_domain1_all.nc'
ds02 = Dataset(infile, 'r')
# obs CLDN
obsfile = '/scratch/leuven/317/vsc31786/nu-wrf-dev/wrfout_nc_files/CLDN_at_domain1_all.nc'
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
PR92H = ds01['PR92H'][:].data
PR92H_ds01 = np.nanmean(PR92H,0)
PR92W = ds01['PR92W'][:].data
PR92W_ds01 = np.nanmean(PR92W,0)
# temporal mean for ds02
CAPExP = ds02['CAPExP'][:].data
CAPExP_ds02 = np.nanmean(CAPExP,0)
LPI = ds02['LPI'][:].data
LPI_ds02 = np.nanmean(LPI,0)
LT3 = ds02['LTG3'][:].data
LT3_ds02 = np.nanmean(LT3,0)
PR92H = ds02['PR92H'][:].data
PR92H_ds02 = np.nanmean(PR92H,0)
PR92W = ds02['PR92W'][:].data
PR92W_ds02 = np.nanmean(PR92W,0)

# prepare fields for panda correlation function, flatten 2D data. 
d = {'CLDN': obs.flatten(),\
     'CAPExP 9km': CAPExP_ds01.flatten(), \
     'LPI 9km': LPI_ds01.flatten(), \
     'LT3 9km': LT3_ds01.flatten(), \
     'PR92H 9km': PR92H_ds01.flatten(), \
     'PR92W 9km': PR92W_ds01.flatten(), \
     'LPI 3km': LPI_ds02.flatten(), \
     'LT3 3km': LT3_ds02.flatten(), \
     'PR92H 3km': PR92H_ds02.flatten(), \
     'PR92W 3km': PR92W_ds02.flatten()\
     }

df = pd.DataFrame(data=d)
# panda dataframe is ready, apply correlation matrix function, 
# default method is pearson
#corr_matrix = df.corr(method='spearman')
# Pearson:
corr_matrix = df.corr()

# hide upper triangle
mask = np.zeros_like(corr_matrix, dtype=np.bool)
mask[np.triu_indices_from(mask)]= True

#f, ax = plt.subplots(figsize=(11, 15))
f, ax = plt.subplots(figsize=(11*0.8, 15*0.8))

heatmap = sns.heatmap(corr_matrix, mask = mask, square = True, linewidths = .5, cmap = 'coolwarm',\
cbar_kws = {'shrink': .4, 'ticks' : [-1, -.5, 0, 0.5, 1]},vmin = -1,vmax = 1,annot = True,annot_kws = {'size': 12})

#add the column names as labels
ax.set_yticklabels(corr_matrix.columns, rotation = 0)
ax.set_xticklabels(corr_matrix.columns)
ax.set_title('Spatial correlation (-) of long-term average')

sns.set_style({'xtick.bottom': True}, {'ytick.left': True})
heatmap.get_figure().savefig(outpath+figname, bbox_inches='tight')

