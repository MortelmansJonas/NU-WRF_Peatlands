# ---------------------------------------------------------------------------------------------
# MODULES
# ---------------------------------------------------------------------------------------------
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import seaborn as sns
import math
from mpl_toolkits.basemap import Basemap


# ---------------------------------------------------------------------------------------------
# LOAD DATA
# ---------------------------------------------------------------------------------------------
ds_d01_G = Dataset('/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/NU-WRF/wrfout_nc_files/d01_diurnal_spat.nc', 'r')
ds_d02_G = Dataset('/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/NU-WRF/wrfout_nc_files/d02_diurnal_spa.nc', 'r')
ds_d01_T = Dataset('/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/NU-WRF/wrfout_nc_files/d01_diurnal_spat_Thompson.nc', 'r')
ds_d02_T = Dataset('/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/NU-WRF/wrfout_nc_files/d02_diurnal_spat_Thompson.nc', 'r')

lat_d01 = ds_d01_T['lat'][:]
lon_d01 = ds_d01_T['lon'][:]

obs = ds_d01_T['Obs'][:]
CAPExP_R_d01_T = ds_d01_T['CAPExP_R'][:]
LPI_d01_T = ds_d01_T['LPI'][:]
LT3_d01_T = ds_d01_T['LTG3'][:]
PR92W_d01_T = ds_d01_T['PR92W'][:]

CAPExP_R_d02_T = ds_d02_T['CAPExP_R'][:]
LPI_d02_T = ds_d02_T['LPI'][:]
LT3_d02_T = ds_d02_T['LTG3'][:]
PR92W_d02_T = ds_d02_T['PR92W'][:]

CAPExP_R_d01_G = ds_d01_G['CAPExP_R'][:]
LPI_d01_G = ds_d01_G['LPI'][:]
LT3_d01_G = ds_d01_G['LTG3'][:]
PR92W_d01_G = ds_d01_G['PR92W'][:]

CAPExP_R_d02_G = ds_d02_G['CAPExP_R'][:]
LPI_d02_G = ds_d02_G['LPI'][:]
LT3_d02_G = ds_d02_G['LTG3'][:]
PR92W_d02_G = ds_d02_G['PR92W'][:]
# ---------------------------------------------------------------------------------------------
# CALCULATE TEMPORAL R FOR EACH PIXEL
# ---------------------------------------------------------------------------------------------
R_LPI_d01_T=np.zeros((8,14))
R_LT3_d01_T=np.zeros((8,14))
R_PR92W_d01_T=np.zeros((8,14))
R_CAPExP_R_d01_T=np.zeros((8,14))
R_LPI_d02_T=np.zeros((8,14))
R_LT3_d02_T=np.zeros((8,14))
R_PR92W_d02_T=np.zeros((8,14))
R_CAPExP_R_d02_T=np.zeros((8,14))
R_LPI_d01_G=np.zeros((8,14))
R_LT3_d01_G=np.zeros((8,14))
R_PR92W_d01_G=np.zeros((8,14))
R_CAPExP_R_d01_G=np.zeros((8,14))
R_LPI_d02_G=np.zeros((8,14))
R_LT3_d02_G=np.zeros((8,14))
R_PR92W_d02_G=np.zeros((8,14))
R_CAPExP_R_d02_G=np.zeros((8,14))

for i in range(0,8):
    print(str(i))
    for j in range(0,14):
        d_LPI_d01_T = {'CLDN': obs[:,i,j],\
                     'LPI 9km': LPI_d01_T[:,i,j]}
        d_LT3_d01_T = {'CLDN': obs[:,i,j], \
                     'LT3 9 km': LT3_d01_T[:,i,j]}
        d_PR92W_d01_T = {'CLDN': obs[:,i,j], \
                     'PR92W 9km': PR92W_d01_T[:,i,j]}
        d_CAPExP_R_d01_T = {'CLDN': obs[:,i,j], \
                     'CAPExP_R 9km': CAPExP_R_d01_T[:,i,j]}
        d_LPI_d02_T = {'CLDN': obs[:,i,j],\
                     'LPI 9km': LPI_d02_T[:,i,j]}
        d_LT3_d02_T = {'CLDN': obs[:,i,j], \
                     'LT3 9km': LT3_d02_T[:,i,j]}
        d_PR92W_d02_T = {'CLDN': obs[:,i,j], \
                     'PR92W 9km': PR92W_d02_T[:,i,j]}
        d_CAPExP_R_d02_T = {'CLDN': obs[:,i,j], \
                     'CAPExP_R 9km': CAPExP_R_d02_T[:,i,j]}
        d_LPI_d01_G = {'CLDN': obs[:,i,j],\
                     'LPI 9km': LPI_d01_G[:,i,j]}
        d_LT3_d01_G = {'CLDN': obs[:,i,j], \
                     'LT3 9 km': LT3_d01_G[:,i,j]}
        d_PR92W_d01_G = {'CLDN': obs[:,i,j], \
                     'PR92W 9km': PR92W_d01_G[:,i,j]}
        d_CAPExP_R_d01_G = {'CLDN': obs[:,i,j], \
                     'CAPExP_R 9km': CAPExP_R_d01_G[:,i,j]}
        d_LPI_d02_G = {'CLDN': obs[:,i,j],\
                     'LPI 9km': LPI_d02_G[:,i,j]}
        d_LT3_d02_G = {'CLDN': obs[:,i,j], \
                     'LT3 9km': LT3_d02_G[:,i,j]}
        d_PR92W_d02_G = {'CLDN': obs[:,i,j], \
                     'PR92W 9km': PR92W_d02_G[:,i,j]}
        d_CAPExP_R_d02_G = {'CLDN': obs[:,i,j], \
                     'CAPExP_R 9km': CAPExP_R_d02_G[:,i,j]}

        df_LPI_d01_T = pd.DataFrame(data=d_LPI_d01_T)
        df_LT3_d01_T = pd.DataFrame(data=d_LT3_d01_T)
        df_PR92W_d01_T = pd.DataFrame(data=d_PR92W_d01_T)
        df_CAPExP_R_d01_T = pd.DataFrame(data=d_CAPExP_R_d01_T)
        df_LPI_d02_T = pd.DataFrame(data=d_LPI_d02_T)
        df_LT3_d02_T = pd.DataFrame(data=d_LT3_d02_T)
        df_PR92W_d02_T = pd.DataFrame(data=d_PR92W_d02_T)
        df_CAPExP_R_d02_T = pd.DataFrame(data=d_CAPExP_R_d02_T)
        corr_LPI_d01_T = df_LPI_d01_T.corr()
        corr_LT3_d01_T = df_LT3_d01_T.corr()
        corr_PR92W_d01_T = df_PR92W_d01_T.corr()
        corr_CAPExP_R_d01_T = df_CAPExP_R_d01_T.corr()
        corr_LPI_d02_T = df_LPI_d02_T.corr()
        corr_LT3_d02_T = df_LT3_d02_T.corr()
        corr_PR92W_d02_T = df_PR92W_d02_T.corr()
        corr_CAPExP_R_d02_T = df_CAPExP_R_d02_T.corr()

        df_LPI_d01_G = pd.DataFrame(data=d_LPI_d01_G)
        df_LT3_d01_G = pd.DataFrame(data=d_LT3_d01_G)
        df_PR92W_d01_G = pd.DataFrame(data=d_PR92W_d01_G)
        df_CAPExP_R_d01_G = pd.DataFrame(data=d_CAPExP_R_d01_G)
        df_LPI_d02_G = pd.DataFrame(data=d_LPI_d02_G)
        df_LT3_d02_G = pd.DataFrame(data=d_LT3_d02_G)
        df_PR92W_d02_G = pd.DataFrame(data=d_PR92W_d02_G)
        df_CAPExP_R_d02_G = pd.DataFrame(data=d_CAPExP_R_d02_G)
        corr_LPI_d01_G = df_LPI_d01_G.corr()
        corr_LT3_d01_G = df_LT3_d01_G.corr()
        corr_PR92W_d01_G = df_PR92W_d01_G.corr()
        corr_CAPExP_R_d01_G = df_CAPExP_R_d01_G.corr()
        corr_LPI_d02_G = df_LPI_d02_G.corr()
        corr_LT3_d02_G = df_LT3_d02_G.corr()
        corr_PR92W_d02_G = df_PR92W_d02_G.corr()
        corr_CAPExP_R_d02_G = df_CAPExP_R_d02_G.corr()

        R_LPI_d01_T[i,j] = corr_LPI_d01_T.iloc[1,0]
        R_LT3_d01_T[i, j] = corr_LT3_d01_T.iloc[1, 0]
        R_PR92W_d01_T[i, j] = corr_PR92W_d01_T.iloc[1, 0]
        R_CAPExP_R_d01_T[i, j] = corr_CAPExP_R_d01_T.iloc[1, 0]
        R_LPI_d02_T[i,j] = corr_LPI_d02_T.iloc[1,0]
        R_LT3_d02_T[i, j] = corr_LT3_d02_T.iloc[1, 0]
        R_PR92W_d02_T[i, j] = corr_PR92W_d02_T.iloc[1, 0]
        R_CAPExP_R_d02_T[i, j] = corr_CAPExP_R_d02_T.iloc[1, 0]
        R_LPI_d01_G[i,j] = corr_LPI_d01_G.iloc[1,0]
        R_LT3_d01_G[i, j] = corr_LT3_d01_G.iloc[1, 0]
        R_PR92W_d01_G[i, j] = corr_PR92W_d01_G.iloc[1, 0]
        R_CAPExP_R_d01_G[i, j] = corr_CAPExP_R_d01_G.iloc[1, 0]
        R_LPI_d02_G[i,j] = corr_LPI_d02_G.iloc[1,0]
        R_LT3_d02_G[i, j] = corr_LT3_d02_G.iloc[1, 0]
        R_PR92W_d02_G[i, j] = corr_PR92W_d02_G.iloc[1, 0]
        R_CAPExP_R_d02_G[i, j] = corr_CAPExP_R_d02_G.iloc[1, 0]

length = len(R_LPI_d01_G.flatten())

df = pd.DataFrame(index=range((20*length)), columns =['sim', 'Parameterization', 'mean'])

df.iloc[0:(4*length),0] = 'THOM -\n9 km'
df.iloc[(4*length):(8*length),0] = 'THOM -\n3 km'
df.iloc[(8*length):(12*length),0] = 'G4ICE -\n9 km'
df.iloc[(12*length):(16*length),0] = 'G4ICE -\n3 km'
lpi_ind = np.asarray([np.arange(0,length),np.arange((4*length),(5*length)),np.arange((8*length),(9*length)),np.arange((12*length),(13*length))])
df.iloc[lpi_ind,1] = 'LPI'
lt3_ind = np.asarray([np.arange(length,(2*length)),np.arange((5*length),(6*length)),np.arange((9*length),(10*length)),np.arange((13*length),(14*length))])
df.iloc[lt3_ind,1] = 'LT'
pr92w_ind = np.asarray([np.arange(2*length,(3*length)),np.arange((6*length),(7*length)),np.arange((10*length),(11*length)),np.arange((14*length),(15*length))])
df.iloc[pr92w_ind,1] = 'PR92W'
capexp_r_ind = np.asarray([np.arange(3*length,(4*length)),np.arange((7*length),(8*length)),np.arange((11*length),(12*length)),np.arange((15*length),(16*length))])
df.iloc[capexp_r_ind,1] = 'CAPExP'

df.iloc[0:length,2] = R_LPI_d01_T.flatten()
df.iloc[length:(2*length),2] = R_LT3_d01_T.flatten()
df.iloc[(2*length):(3*length),2] = R_PR92W_d01_T.flatten()
df.iloc[(3*length):(4*length),2] = R_CAPExP_R_d01_T.flatten()

df.iloc[(4*length):(5*length),2] = R_LPI_d02_T.flatten()
df.iloc[(5*length):(6*length),2] = R_LT3_d02_T.flatten()
df.iloc[(6*length):(7*length),2] = R_PR92W_d02_T.flatten()
df.iloc[(7*length):(8*length),2] = R_CAPExP_R_d02_T.flatten()

df.iloc[(8*length):(9*length),2] = R_LPI_d01_G.flatten()
df.iloc[(9*length):(10*length),2] = R_LT3_d01_G.flatten()
df.iloc[(10*length):(11*length),2] = R_PR92W_d01_G.flatten()
df.iloc[(11*length):(12*length),2] = R_CAPExP_R_d01_G.flatten()

df.iloc[(12*length):(13*length),2] = R_LPI_d02_G.flatten()
df.iloc[(13*length):(14*length),2] = R_LT3_d02_G.flatten()
df.iloc[(14*length):(15*length),2] = R_PR92W_d02_G.flatten()
df.iloc[(15*length):(16*length),2] = R_CAPExP_R_d02_G.flatten()

colorblind = sns.color_palette('colorblind')

fig, axes = plt.subplots(1,1, figsize=(15.3, 8.27), dpi=150)
palette = [colorblind[0], colorblind[1], colorblind[4], colorblind[3]]
ax = sns.boxplot(x='sim', y='mean', hue='Parameterization', ax=axes, data=df, palette=palette)
axes.set_title('Diurnal Cycle', fontsize=14)
axes.legend([],[], frameon=False)
axes.set_ylabel('R (-)', fontsize=12)
axes.set_xlabel('')
plt.legend(bbox_to_anchor=(0.75, -0.1), ncol=5,fontsize=12)
axes.tick_params(axis='both', which='major', labelsize=12)
plt.show()
