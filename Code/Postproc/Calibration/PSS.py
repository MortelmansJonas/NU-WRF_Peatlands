#!/usr/bin/env python
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
ds_d01_T = Dataset('/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/NU-WRF/NU-WRF/wrfout_nc_files/d01_all_Thompson.nc','r')
ds_d02_T = Dataset('/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/NU-WRF/NU-WRF/wrfout_nc_files/domain2_at_domain1_all_Thompson.nc','r')
ds_obs = Dataset('/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/NU-WRF/NU-WRF/wrfout_nc_files/CLDN_at_domain1_all_v4.nc')
ds_d01_G = Dataset('/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/NU-WRF/NU-WRF/wrfout_nc_files/d01_all.nc','r')
ds_d02_G = Dataset('/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/NU-WRF/NU-WRF/wrfout_nc_files/domain2_at_domain1_all.nc','r')

lat = ds_d01_T['lat'][:]
lon = ds_d01_T['lon'][:]

L = ds_obs['Flashdensity_CC'][:,:,:].data.flatten() + ds_obs['Flashdensity_CG'][:,:,:].data.flatten()

# ---------------------------------------------------------------------------------------------
# 'The lowest values of the parameterization's output are set to 0'
# ---------------------------------------------------------------------------------------------
LPI_d01_T = ds_d01_T['LPI'][:]
LPI_d01_T[LPI_d01_T == np.nanmin(LPI_d01_T)] = 0

LTG3_d01_T = ds_d01_T['LTG3'][:]
LTG3_d01_T[LTG3_d01_T == np.nanmin(LTG3_d01_T)] = 0

PR92W_d01_T = ds_d01_T['PR92W'][:]
PR92W_d01_T[PR92W_d01_T == np.nanmin(PR92W_d01_T)] = 0

CAPExP_R_d01_T = ds_d01_T['CAPExP_R'][:]
CAPExP_R_d01_T[CAPExP_R_d01_T == np.nanmin(CAPExP_R_d01_T)] = 0

LPI_d02_T = ds_d02_T['LPI'][:]
LPI_d02_T[LPI_d02_T == np.nanmin(LPI_d02_T)] = 0

LTG3_d02_T = ds_d02_T['LTG3'][:]
LTG3_d02_T[LTG3_d02_T == np.nanmin(LTG3_d02_T)] = 0

PR92W_d02_T = ds_d02_T['PR92W'][:]
PR92W_d02_T[PR92W_d02_T == np.nanmin(PR92W_d02_T)] = 0

CAPExP_R_d02_T = ds_d02_T['CAPExP_R'][:]
CAPExP_R_d02_T[CAPExP_R_d02_T == np.nanmin(CAPExP_R_d02_T)] = 0

LPI_d01_G = ds_d01_G['LPI'][:]
LPI_d01_G[LPI_d01_G == np.nanmin(LPI_d01_G)] = 0

LTG3_d01_G = ds_d01_G['LTG3'][:]
LTG3_d01_G[LTG3_d01_G == np.nanmin(LTG3_d01_G)] = 0

PR92W_d01_G = ds_d01_G['PR92W'][:]
PR92W_d01_G[PR92W_d01_G == np.nanmin(PR92W_d01_G)] = 0

CAPExP_R_d01_G = ds_d01_G['CAPExP_R'][:]
CAPExP_R_d01_G[CAPExP_R_d01_G == np.nanmin(CAPExP_R_d01_G)] = 0

LPI_d02_G = ds_d02_G['LPI'][:]
LPI_d02_G[LPI_d02_G == np.nanmin(LPI_d02_G)] = 0

LTG3_d02_G = ds_d02_G['LTG3'][:]
LTG3_d02_G[LTG3_d02_G == np.nanmin(LTG3_d02_G)] = 0

PR92W_d02_G = ds_d02_G['PR92W'][:]
PR92W_d02_G[PR92W_d02_G == np.nanmin(PR92W_d02_G)] = 0

CAPExP_R_d02_G = ds_d02_G['CAPExP_R'][:]
CAPExP_R_d02_G[CAPExP_R_d02_G == np.nanmin(CAPExP_R_d02_G)] = 0

# -------------------------------------------------------------------------------------------------------------
# 'Second, a simple linear model is built that relates the observed flash rates to the parameterization output'
# -------------------------------------------------------------------------------------------------------------
# First determine cutoff value 'c', done similarly to the sorting in 'Calibration_2.0.py'
Obs = np.sort(L[L!=0])[::-1]

# Sort model and take n highest values (with n = number of observations)
c = Obs.shape[0]
Tot_flashes = np.nansum(Obs)

LPI_sorted_d01_T = np.sort(LPI_d01_T.flatten())[::-1]
LPI_c_d01_T = LPI_sorted_d01_T[0:c]
LPI_d01_T_new = np.where(LPI_d01_T>=LPI_sorted_d01_T[c],LPI_d01_T, np.nan)
print('cutoff value LPI_d01_T = ', LPI_sorted_d01_T[c])

LTG3_sorted_d01_T = np.sort(LTG3_d01_T.flatten())[::-1]
# exclude 0 values from calibration:
# get first zero value
c2 = np.where(LTG3_sorted_d01_T==0)[0][0]
LTG3_c_d01_T = LTG3_sorted_d01_T[0:c2]
LTG3_d01_T_new = np.where(LTG3_d01_T>=LTG3_sorted_d01_T[c2],LTG3_d01_T, np.nan)
print('cutoff value LTG3_d01_T = ', LTG3_sorted_d01_T[c2])

PR92W_sorted_d01_T = np.sort(PR92W_d01_T.flatten())[::-1]
PR92W_c_d01_T = PR92W_sorted_d01_T[0:c]
PR92W_d01_T_new = np.where(PR92W_d01_T>=PR92W_sorted_d01_T[c],PR92W_d01_T, np.nan)
print('cutoff value PR92W_d01_T = ', PR92W_sorted_d01_T[c])

CAPExP_R_sorted_d01_T = np.sort(CAPExP_R_d01_T.flatten())[::-1]
CAPExP_R_c_d01_T = CAPExP_R_sorted_d01_T[0:c]
CAPExP_R_d01_T_new = np.where(CAPExP_R_d01_T>=CAPExP_R_sorted_d01_T[c],CAPExP_R_d01_T, np.nan)
print('cutoff value CAPExP_R_d01_T = ', CAPExP_R_sorted_d01_T[c])

LPI_sorted_d02_T = np.sort(LPI_d02_T.flatten())[::-1]
LPI_c_d02_T = LPI_sorted_d02_T[0:c]
LPI_d02_T_new = np.where(LPI_d02_T>=LPI_sorted_d02_T[c],LPI_d02_T, np.nan)
print('cutoff value LPI_d02_T = ', LPI_sorted_d02_T[c])

LTG3_sorted_d02_T = np.sort(LTG3_d02_T.flatten())[::-1]
# exclude 0 values from calibration:
# get first zero value
c3 = np.where(LTG3_sorted_d02_T==0)[0][0]
LTG3_c_d02_T = LTG3_sorted_d02_T[0:c3]
LTG3_d02_T_new = np.where(LTG3_d02_T>=LTG3_sorted_d02_T[c3],LTG3_d02_T, np.nan)
print('cutoff value LTG3_d02_T = ', LTG3_sorted_d02_T[c3])

PR92W_sorted_d02_T = np.sort(PR92W_d02_T.flatten())[::-1]
PR92W_c_d02_T = PR92W_sorted_d02_T[0:c]
PR92W_d02_T_new = np.where(PR92W_d02_T>=PR92W_sorted_d02_T[c],PR92W_d02_T, np.nan)
print('cutoff value PR92W_d02_T = ', PR92W_sorted_d02_T[c])

CAPExP_R_sorted_d02_T = np.sort(CAPExP_R_d02_T.flatten())[::-1]
CAPExP_R_c_d02_T = CAPExP_R_sorted_d02_T[0:c]
CAPExP_R_d02_T_new = np.where(CAPExP_R_d02_T>=CAPExP_R_sorted_d02_T[c],CAPExP_R_d02_T, np.nan)
print('cutoff value CAPExP_R_d02_T = ', CAPExP_R_sorted_d02_T[c])

LPI_sorted_d01_G = np.sort(LPI_d01_G.flatten())[::-1]
LPI_c_d01_G = LPI_sorted_d01_G[0:c]
LPI_d01_G_new = np.where(LPI_d01_G>=LPI_sorted_d01_G[c],LPI_d01_G, np.nan)
print('cutoff value LPI_d01_G = ', LPI_sorted_d01_G[c])

LTG3_sorted_d01_G = np.sort(LTG3_d01_G.flatten())[::-1]
# exclude 0 values from calibration:
# get first zero value
c2 = np.where(LTG3_sorted_d01_G==0)[0][0]
LTG3_c_d01_G = LTG3_sorted_d01_G[0:c2]
LTG3_d01_G_new = np.where(LTG3_d01_G>=LTG3_sorted_d01_G[c2],LTG3_d01_G, np.nan)
print('cutoff value LTG3_d01_G = ', LTG3_sorted_d01_G[c2])

PR92W_sorted_d01_G = np.sort(PR92W_d01_G.flatten())[::-1]
PR92W_c_d01_G = PR92W_sorted_d01_G[0:c]
PR92W_d01_G_new = np.where(PR92W_d01_G>=PR92W_sorted_d01_G[c],PR92W_d01_G, np.nan)
print('cutoff value PR92W_d01_G = ', PR92W_sorted_d01_G[c])

CAPExP_R_sorted_d01_G = np.sort(CAPExP_R_d01_G.flatten())[::-1]
CAPExP_R_c_d01_G = CAPExP_R_sorted_d01_G[0:c]
CAPExP_R_d01_G_new = np.where(CAPExP_R_d01_G>=CAPExP_R_sorted_d01_G[c],CAPExP_R_d01_G, np.nan)
print('cutoff value CAPExP_R_d01_G = ', CAPExP_R_sorted_d01_G[c])

LPI_sorted_d02_G = np.sort(LPI_d02_G.flatten())[::-1]
LPI_c_d02_G = LPI_sorted_d02_G[0:c]
LPI_d02_G_new = np.where(LPI_d02_G>=LPI_sorted_d02_G[c],LPI_d02_G, np.nan)
print('cutoff value LPI_d02_G = ', LPI_sorted_d02_G[c])

LTG3_sorted_d02_G = np.sort(LTG3_d02_G.flatten())[::-1]
# exclude 0 values from calibration:
# get first zero value
c3 = np.where(LTG3_sorted_d02_G==0)[0][0]
LTG3_c_d02_G = LTG3_sorted_d02_G[0:c3]
LTG3_d02_G_new = np.where(LTG3_d02_G>=LTG3_sorted_d02_G[c3],LTG3_d02_G, np.nan)
print('cutoff value LTG3_d02_G = ', LTG3_sorted_d02_G[c3])

PR92W_sorted_d02_G = np.sort(PR92W_d02_G.flatten())[::-1]
PR92W_c_d02_G = PR92W_sorted_d02_G[0:c]
PR92W_d02_G_new = np.where(PR92W_d02_G>=PR92W_sorted_d02_G[c],PR92W_d02_G, np.nan)
print('cutoff value PR92W_d02_G = ', PR92W_sorted_d02_G[c])

CAPExP_R_sorted_d02_G = np.sort(CAPExP_R_d02_G.flatten())[::-1]
CAPExP_R_c_d02_G = CAPExP_R_sorted_d02_G[0:c]
CAPExP_R_d02_G_new = np.where(CAPExP_R_d02_G>=CAPExP_R_sorted_d02_G[c],CAPExP_R_d02_G, np.nan)
print('cutoff value CAPExP_R_d02_G = ', CAPExP_R_sorted_d02_G[c])
# -------------------------------------------------------------------------------------------------------------
# X_adj = aX+b if X >= c
# -------------------------------------------------------------------------------------------------------------
print('calibration')
reg_LPI_d01_T = lr().fit(LPI_c_d01_T.reshape((-1,1)), Obs)
LPI_d01_T_adj =np.add(reg_LPI_d01_T.intercept_, np.multiply(reg_LPI_d01_T.coef_,LPI_d01_T_new))
LPI_d01_T_adj[np.isnan(LPI_d01_T_adj)] = 0

# reg_LTG3_d01_T = lr().fit(LTG3_c_d01_T.reshape((-1,1)), Obs)
# cut obs if needed
c4=np.min(np.array([Obs.shape[0],LTG3_c_d01_T.shape[0]]))
reg_LTG3_d01_T = lr().fit(LTG3_c_d01_T.reshape((-1,1)), Obs[0:c4])
LTG3_d01_T_adj =np.add(reg_LTG3_d01_T.intercept_, np.multiply(reg_LTG3_d01_T.coef_,LTG3_d01_T_new))
# zero should stay zero
LTG3_d01_T_adj[LTG3_d01_T_new==0] = 0
LTG3_d01_T_adj[np.isnan(LTG3_d01_T_adj)] = 0
# scale adj lightning if insufficient number of lightning:
if np.nansum(LTG3_d01_T_adj)<Tot_flashes:
    LTG3_d01_T_adj = Tot_flashes/np.nansum(LTG3_d01_T_adj) * LTG3_d01_T_adj

reg_PR92W_d01_T = lr().fit(PR92W_c_d01_T.reshape((-1,1)), Obs)
PR92W_d01_T_adj =np.add(reg_PR92W_d01_T.intercept_, np.multiply(reg_PR92W_d01_T.coef_,PR92W_d01_T_new))
PR92W_d01_T_adj[np.isnan(PR92W_d01_T_adj)] = 0

reg_CAPExP_R_d01_T = lr().fit(CAPExP_R_c_d01_T.reshape((-1,1)), Obs)
CAPExP_R_d01_T_adj =np.add(reg_CAPExP_R_d01_T.intercept_, np.multiply(reg_CAPExP_R_d01_T.coef_,CAPExP_R_d01_T_new))
CAPExP_R_d01_T_adj[np.isnan(CAPExP_R_d01_T_adj)] = 0

reg_LPI_d02_T = lr().fit(LPI_c_d02_T.reshape((-1,1)), Obs)
LPI_d02_T_adj =np.add(reg_LPI_d02_T.intercept_, np.multiply(reg_LPI_d02_T.coef_,LPI_d02_T_new))
LPI_d02_T_adj[np.isnan(LPI_d02_T_adj)] = 0

# reg_LTG3_d02_T = lr().fit(LTG3_c_d02_T.reshape((-1,1)), Obs)
# cut obs if needed
c5=np.min(np.array([Obs.shape[0],LTG3_c_d02_T.shape[0]]))
reg_LTG3_d02_T = lr().fit(LTG3_c_d02_T.reshape((-1,1)), Obs[0:c5])
LTG3_d02_T_adj =np.add(reg_LTG3_d02_T.intercept_, np.multiply(reg_LTG3_d02_T.coef_,LTG3_d02_T_new))
# zero should stay zero
LTG3_d02_T_adj[LTG3_d02_T_new==0] = 0
LTG3_d02_T_adj[np.isnan(LTG3_d02_T_adj)] = 0
# scale adj lightning if insufficient number of lightning:
if np.nansum(LTG3_d02_T_adj)<Tot_flashes:
    LTG3_d02_T_adj = Tot_flashes/np.nansum(LTG3_d02_T_adj) * LTG3_d02_T_adj

reg_PR92W_d02_T = lr().fit(PR92W_c_d02_T.reshape((-1,1)), Obs)
PR92W_d02_T_adj =np.add(reg_PR92W_d02_T.intercept_, np.multiply(reg_PR92W_d02_T.coef_,PR92W_d02_T_new))
PR92W_d02_T_adj[np.isnan(PR92W_d02_T_adj)] = 0

reg_CAPExP_R_d02_T = lr().fit(CAPExP_R_c_d02_T.reshape((-1,1)), Obs)
CAPExP_R_d02_T_adj =np.add(reg_CAPExP_R_d02_T.intercept_, np.multiply(reg_CAPExP_R_d02_T.coef_,CAPExP_R_d02_T_new))
CAPExP_R_d02_T_adj[np.isnan(CAPExP_R_d02_T_adj)] = 0

reg_LPI_d01_G = lr().fit(LPI_c_d01_G.reshape((-1,1)), Obs)
LPI_d01_G_adj =np.add(reg_LPI_d01_G.intercept_, np.multiply(reg_LPI_d01_G.coef_,LPI_d01_G_new))
LPI_d01_G_adj[np.isnan(LPI_d01_G_adj)] = 0

# reg_LTG3_d01_G = lr().fit(LTG3_c_d01_G.reshape((-1,1)), Obs)
# cut obs if needed
c4=np.min(np.array([Obs.shape[0],LTG3_c_d01_G.shape[0]]))
reg_LTG3_d01_G = lr().fit(LTG3_c_d01_G.reshape((-1,1)), Obs[0:c4])
LTG3_d01_G_adj =np.add(reg_LTG3_d01_G.intercept_, np.multiply(reg_LTG3_d01_G.coef_,LTG3_d01_G_new))
# zero should stay zero
LTG3_d01_G_adj[LTG3_d01_G_new==0] = 0
LTG3_d01_G_adj[np.isnan(LTG3_d01_G_adj)] = 0
# scale adj lightning if insufficient number of lightning:
if np.nansum(LTG3_d01_G_adj)<Tot_flashes:
    LTG3_d01_G_adj = Tot_flashes/np.nansum(LTG3_d01_G_adj) * LTG3_d01_G_adj

reg_PR92W_d01_G = lr().fit(PR92W_c_d01_G.reshape((-1,1)), Obs)
PR92W_d01_G_adj =np.add(reg_PR92W_d01_G.intercept_, np.multiply(reg_PR92W_d01_G.coef_,PR92W_d01_G_new))
PR92W_d01_G_adj[np.isnan(PR92W_d01_G_adj)] = 0

reg_CAPExP_R_d01_G = lr().fit(CAPExP_R_c_d01_G.reshape((-1,1)), Obs)
CAPExP_R_d01_G_adj =np.add(reg_CAPExP_R_d01_G.intercept_, np.multiply(reg_CAPExP_R_d01_G.coef_,CAPExP_R_d01_G_new))
CAPExP_R_d01_G_adj[np.isnan(CAPExP_R_d01_G_adj)] = 0

reg_LPI_d02_G = lr().fit(LPI_c_d02_G.reshape((-1,1)), Obs)
LPI_d02_G_adj =np.add(reg_LPI_d02_G.intercept_, np.multiply(reg_LPI_d02_G.coef_,LPI_d02_G_new))
LPI_d02_G_adj[np.isnan(LPI_d02_G_adj)] = 0

# reg_LTG3_d02_G = lr().fit(LTG3_c_d02_G.reshape((-1,1)), Obs)
# cut obs if needed
c5=np.min(np.array([Obs.shape[0],LTG3_c_d02_G.shape[0]]))
reg_LTG3_d02_G = lr().fit(LTG3_c_d02_G.reshape((-1,1)), Obs[0:c5])
LTG3_d02_G_adj =np.add(reg_LTG3_d02_G.intercept_, np.multiply(reg_LTG3_d02_G.coef_,LTG3_d02_G_new))
# zero should stay zero
LTG3_d02_G_adj[LTG3_d02_G_new==0] = 0
LTG3_d02_G_adj[np.isnan(LTG3_d02_G_adj)] = 0
# scale adj lightning if insufficient number of lightning:
if np.nansum(LTG3_d02_G_adj)<Tot_flashes:
    LTG3_d02_G_adj = Tot_flashes/np.nansum(LTG3_d02_G_adj) * LTG3_d02_G_adj

reg_PR92W_d02_G = lr().fit(PR92W_c_d02_G.reshape((-1,1)), Obs)
PR92W_d02_G_adj =np.add(reg_PR92W_d02_G.intercept_, np.multiply(reg_PR92W_d02_G.coef_,PR92W_d02_G_new))
PR92W_d02_G_adj[np.isnan(PR92W_d02_G_adj)] = 0

reg_CAPExP_R_d02_G = lr().fit(CAPExP_R_c_d02_G.reshape((-1,1)), Obs)
CAPExP_R_d02_G_adj =np.add(reg_CAPExP_R_d02_G.intercept_, np.multiply(reg_CAPExP_R_d02_G.coef_,CAPExP_R_d02_G_new))
CAPExP_R_d02_G_adj[np.isnan(CAPExP_R_d02_G_adj)] = 0

print('sorting')
LPI_adj_sorted_d01_T = np.sort(LPI_d01_T_adj.flatten())[::-1]
LTG3_adj_sorted_d01_T = np.sort(LTG3_d01_T_adj.flatten())[::-1]
PR92W_adj_sorted_d01_T = np.sort(PR92W_d01_T_adj.flatten())[::-1]
CAPExP_R_adj_sorted_d01_T = np.sort(CAPExP_R_d01_T_adj.flatten())[::-1]
LPI_adj_sorted_d02_T = np.sort(LPI_d02_T_adj.flatten())[::-1]
LTG3_adj_sorted_d02_T = np.sort(LTG3_d02_T_adj.flatten())[::-1]
PR92W_adj_sorted_d02_T = np.sort(PR92W_d02_T_adj.flatten())[::-1]
CAPExP_R_adj_sorted_d02_T = np.sort(CAPExP_R_d02_T_adj.flatten())[::-1]

LPI_adj_sorted_d01_G = np.sort(LPI_d01_G_adj.flatten())[::-1]
LTG3_adj_sorted_d01_G = np.sort(LTG3_d01_G_adj.flatten())[::-1]
PR92W_adj_sorted_d01_G = np.sort(PR92W_d01_G_adj.flatten())[::-1]
CAPExP_R_adj_sorted_d01_G = np.sort(CAPExP_R_d01_G_adj.flatten())[::-1]
LPI_adj_sorted_d02_G = np.sort(LPI_d02_G_adj.flatten())[::-1]
LTG3_adj_sorted_d02_G = np.sort(LTG3_d02_G_adj.flatten())[::-1]
PR92W_adj_sorted_d02_G = np.sort(PR92W_d02_G_adj.flatten())[::-1]
CAPExP_R_adj_sorted_d02_G = np.sort(CAPExP_R_d02_G_adj.flatten())[::-1]

print('sum')
sum_lpi_d01_T = 0
sum_ltg3_d01_T = 0
sum_pr92w_d01_T = 0
sum_capexp_r_d01_T = 0
sum_lpi_d02_T = 0
sum_ltg3_d02_T = 0
sum_pr92w_d02_T = 0
sum_capexp_r_d02_T = 0
sum_lpi_d01_G = 0
sum_ltg3_d01_G = 0
sum_pr92w_d01_G = 0
sum_capexp_r_d01_G = 0
sum_lpi_d02_G = 0
sum_ltg3_d02_G = 0
sum_pr92w_d02_G = 0
sum_capexp_r_d02_G = 0
i = 0

print('while')
while sum_lpi_d01_T < Tot_flashes:
    sum_lpi_d01_T = sum_lpi_d01_T + LPI_adj_sorted_d01_T[i]
    i += 1
else:
    LPI_adj_sorted_d01_T[i:] = 0
    LPI_adj_sorted_d01_T.reshape((13248,63,109))
i = 0
while (sum_ltg3_d01_T < Tot_flashes) & (i < len(LTG3_adj_sorted_d01_T)):
    sum_ltg3_d01_T = sum_ltg3_d01_T + LTG3_adj_sorted_d01_T[i]
    i += 1
else:
    LTG3_adj_sorted_d01_T[i:] = 0
    LTG3_adj_sorted_d01_T.reshape((13248, 63, 109))
i=0
while sum_pr92w_d01_T < Tot_flashes:
    sum_pr92w_d01_T = sum_pr92w_d01_T + PR92W_adj_sorted_d01_T[i]
    i += 1
else:
    PR92W_adj_sorted_d01_T[i:] = 0
    PR92W_adj_sorted_d01_T.reshape((13248, 63, 109))
i=0
while sum_capexp_r_d01_T < Tot_flashes:
    sum_capexp_r_d01_T = sum_capexp_r_d01_T + CAPExP_R_adj_sorted_d01_T[i]
    i += 1
else:
    CAPExP_R_adj_sorted_d01_T[i:] = 0
    CAPExP_R_adj_sorted_d01_T.reshape((13248, 63, 109))
i=0
while sum_lpi_d02_T < Tot_flashes:
    sum_lpi_d02_T = sum_lpi_d02_T + LPI_adj_sorted_d02_T[i]
    i += 1
else:
    LPI_adj_sorted_d02_T[i:] = 0
    LPI_adj_sorted_d02_T.reshape((13248, 63, 109))
i=0
while (sum_ltg3_d02_T < Tot_flashes) & (i < len(LTG3_adj_sorted_d02_T)):
    sum_ltg3_d02_T = sum_ltg3_d02_T + LTG3_adj_sorted_d02_T[i]
    i += 1
else:
    LTG3_adj_sorted_d02_T[i:] = 0
    LTG3_adj_sorted_d02_T.reshape((13248, 63, 109))
i=0
while sum_pr92w_d02_T < Tot_flashes:
    sum_pr92w_d02_T = sum_pr92w_d02_T + PR92W_adj_sorted_d02_T[i]
    i += 1
else:
    PR92W_adj_sorted_d02_T[i:] = 0
    PR92W_adj_sorted_d02_T.reshape((13248, 63, 109))
i=0
while sum_capexp_r_d02_T < Tot_flashes:
    sum_capexp_r_d02_T = sum_capexp_r_d02_T + CAPExP_R_adj_sorted_d02_T[i]
    i += 1
else:
    CAPExP_R_adj_sorted_d02_T[i:] = 0
    CAPExP_R_adj_sorted_d02_T.reshape((13248, 63, 109))
i=0
while (sum_lpi_d01_G < Tot_flashes) & (i < len(LPI_adj_sorted_d01_G)):
    sum_lpi_d01_G = sum_lpi_d01_G + LPI_adj_sorted_d01_G[i]
    i += 1
else:
    LPI_adj_sorted_d01_G[i:] = 0
    LPI_adj_sorted_d01_G.reshape((13248,63,109))
i = 0
while (sum_ltg3_d01_G < Tot_flashes) & (i < len(LTG3_adj_sorted_d01_G)):
    sum_ltg3_d01_G = sum_ltg3_d01_G + LTG3_adj_sorted_d01_G[i]
    i += 1
else:
    LTG3_adj_sorted_d01_G[i:] = 0
    LTG3_adj_sorted_d01_G.reshape((13248, 63, 109))
i=0
while sum_pr92w_d01_G < Tot_flashes:
    sum_pr92w_d01_G = sum_pr92w_d01_G + PR92W_adj_sorted_d01_G[i]
    i += 1
else:
    PR92W_adj_sorted_d01_G[i:] = 0
    PR92W_adj_sorted_d01_G.reshape((13248, 63, 109))
i=0
while sum_capexp_r_d01_G < Tot_flashes:
    sum_capexp_r_d01_G = sum_capexp_r_d01_G + CAPExP_R_adj_sorted_d01_G[i]
    i += 1
else:
    CAPExP_R_adj_sorted_d01_G[i:] = 0
    CAPExP_R_adj_sorted_d01_G.reshape((13248, 63, 109))
i=0
while sum_lpi_d02_G < Tot_flashes:
    sum_lpi_d02_G = sum_lpi_d02_G + LPI_adj_sorted_d02_G[i]
    i += 1
else:
    LPI_adj_sorted_d02_G[i:] = 0
    LPI_adj_sorted_d02_G.reshape((13248, 63, 109))
i=0
while (sum_ltg3_d02_G < Tot_flashes) & (i < len(LTG3_adj_sorted_d02_G)):
    sum_ltg3_d02_G = sum_ltg3_d02_G + LTG3_adj_sorted_d02_G[i]
    i += 1
else:
    LTG3_adj_sorted_d02_G[i:] = 0
    LTG3_adj_sorted_d02_G.reshape((13248, 63, 109))
i=0
while sum_pr92w_d02_G < Tot_flashes:
    sum_pr92w_d02_G = sum_pr92w_d02_G + PR92W_adj_sorted_d02_G[i]
    i += 1
else:
    PR92W_adj_sorted_d02_G[i:] = 0
    PR92W_adj_sorted_d02_G.reshape((13248, 63, 109))
i=0
while sum_capexp_r_d02_G < Tot_flashes:
    sum_capexp_r_d02_G = sum_capexp_r_d02_G + CAPExP_R_adj_sorted_d02_G[i]
    i += 1
else:
    CAPExP_R_adj_sorted_d02_G[i:] = 0
    CAPExP_R_adj_sorted_d02_G.reshape((13248, 63, 109))

# -------------------------------------------------------------------------------------------------------------
# Plots to check
# -------------------------------------------------------------------------------------------------------------
print('plots')
LPI_adj_sorted_d01_T[LPI_adj_sorted_d01_T == 0] = np.nan
LPI_adj_sorted_d02_T[LPI_adj_sorted_d02_T == 0] = np.nan
LTG3_adj_sorted_d01_T[LTG3_adj_sorted_d01_T == 0] = np.nan
LTG3_adj_sorted_d02_T[LTG3_adj_sorted_d02_T == 0] = np.nan
PR92W_adj_sorted_d01_T[PR92W_adj_sorted_d01_T == 0] = np.nan
PR92W_adj_sorted_d02_T[PR92W_adj_sorted_d02_T == 0] = np.nan
CAPExP_R_adj_sorted_d02_T[CAPExP_R_adj_sorted_d02_T == 0] = np.nan
CAPExP_R_adj_sorted_d01_T[CAPExP_R_adj_sorted_d01_T == 0] = np.nan
LPI_adj_sorted_d01_G[LPI_adj_sorted_d01_G == 0] = np.nan
LPI_adj_sorted_d02_G[LPI_adj_sorted_d02_G == 0] = np.nan
LTG3_adj_sorted_d01_G[LTG3_adj_sorted_d01_G == 0] = np.nan
LTG3_adj_sorted_d02_G[LTG3_adj_sorted_d02_G == 0] = np.nan
PR92W_adj_sorted_d01_G[PR92W_adj_sorted_d01_G == 0] = np.nan
PR92W_adj_sorted_d02_G[PR92W_adj_sorted_d02_G == 0] = np.nan
CAPExP_R_adj_sorted_d02_G[CAPExP_R_adj_sorted_d02_G == 0] = np.nan
CAPExP_R_adj_sorted_d01_G[CAPExP_R_adj_sorted_d01_G == 0] = np.nan
# ---------------------------------------------------------------------------------------------
# FREQUENCY IFO HOURLY GRID FLASH DENSITY
# ---------------------------------------------------------------------------------------------
print('frequency')
# specify bin edges
my_bins = np.linspace(-1,2,20)
# DOMAIN 1
n, bins_log10, patches = plt.hist(np.log10(L[L!=0]),bins = my_bins)
plt.close()
n_CAPExP_R_T, bins_CAPExP_R_log10_T, patches_CAPExP_R_T = plt.hist(np.log10(CAPExP_R_adj_sorted_d01_T[CAPExP_R_adj_sorted_d01_T!=0]),bins=my_bins)
plt.close()
n_lpi_T, bins_lpi_log10_T, patches_lpi_T = plt.hist(np.log10(LPI_adj_sorted_d01_T[LPI_adj_sorted_d01_T!=0]),bins=my_bins)
plt.close()
n_LTG3_T, bins_LTG3_log10_T, patches_LTG3_T = plt.hist(np.log10(LTG3_adj_sorted_d01_T[LTG3_adj_sorted_d01_T!=0]),bins=my_bins)
plt.close()
n_PR92W_T, bins_PR92W_log10_T, patches_PR92W_T = plt.hist(np.log10(PR92W_adj_sorted_d01_T[PR92W_adj_sorted_d01_T!=0]),bins=my_bins)
plt.close()
n_CAPExP_R_G, bins_CAPExP_R_log10_G, patches_CAPExP_R_G = plt.hist(np.log10(CAPExP_R_adj_sorted_d01_G[CAPExP_R_adj_sorted_d01_G!=0]),bins=my_bins)
plt.close()
n_lpi_G, bins_lpi_log10_G, patches_lpi_G = plt.hist(np.log10(LPI_adj_sorted_d01_G[LPI_adj_sorted_d01_G!=0]),bins=my_bins)
plt.close()
n_LTG3_G, bins_LTG3_log10_G, patches_LTG3_G = plt.hist(np.log10(LTG3_adj_sorted_d01_G[LTG3_adj_sorted_d01_G!=0]),bins=my_bins)
plt.close()
n_PR92W_G, bins_PR92W_log10_G, patches_PR92W_G = plt.hist(np.log10(PR92W_adj_sorted_d01_G[PR92W_adj_sorted_d01_G!=0]),bins=my_bins)
plt.close()

# set bins with 0 to the left to np.nan for a line plot, lines should just end to the left when there are no values
print('0 to nan')
for i in range(len(n)):
    if n[i]==0:
        n[i] = np.nan
    else:
        break
for i in range(len(n)):
    if n_CAPExP_R_T[i] == 0:
        n_CAPExP_R_T[i] = np.nan
    else:
        break
for i in range(len(n)):
    if n_lpi_T[i] == 0:
        n_lpi_T[i] = np.nan
    else:
        break
for i in range(len(n)):
    if n_LTG3_T[i] == 0:
        n_LTG3_T[i] = np.nan
    else:
        break
for i in range(len(n)):
    if n_PR92W_T[i] == 0:
        n_PR92W_T[i] = np.nan
    else:
        break
for i in range(len(n)):
    if n_CAPExP_R_G[i] == 0:
        n_CAPExP_R_G[i] = np.nan
    else:
        break
for i in range(len(n)):
    if n_lpi_G[i] == 0:
        n_lpi_G[i] = np.nan
    else:
        break
for i in range(len(n)):
    if n_LTG3_G[i] == 0:
        n_LTG3_G[i] = np.nan
    else:
        break
for i in range(len(n)):
    if n_PR92W_G[i] == 0:
        n_PR92W_G[i] = np.nan
    else:
        break

print('avg bins')
bins_log10_centered = bins_log10[0:-1]+0.5*(bins_log10[1:]-bins_log10[0:-1])
bins_centered = 10**bins_log10_centered

bins_CAPExP_R_log10_centered_T = bins_CAPExP_R_log10_T[0:-1]+0.5*(bins_CAPExP_R_log10_T[1:]-bins_CAPExP_R_log10_T[0:-1])
bins_CAPExP_R_centered_T = 10**bins_CAPExP_R_log10_centered_T

bins_lpi_log10_centered_T = bins_lpi_log10_T[0:-1]+0.5*(bins_lpi_log10_T[1:]-bins_lpi_log10_T[0:-1])
bins_lpi_centered_T = 10**bins_lpi_log10_centered_T

bins_LTG3_log10_centered_T = bins_LTG3_log10_T[0:-1]+0.5*(bins_LTG3_log10_T[1:]-bins_LTG3_log10_T[0:-1])
bins_LTG3_centered_T = 10**bins_LTG3_log10_centered_T

bins_PR92W_log10_centered_T = bins_PR92W_log10_T[0:-1]+0.5*(bins_PR92W_log10_T[1:]-bins_PR92W_log10_T[0:-1])
bins_PR92W_centered_T = 10**bins_PR92W_log10_centered_T

bins_CAPExP_R_log10_centered_G = bins_CAPExP_R_log10_G[0:-1]+0.5*(bins_CAPExP_R_log10_G[1:]-bins_CAPExP_R_log10_G[0:-1])
bins_CAPExP_R_centered_G = 10**bins_CAPExP_R_log10_centered_G

bins_lpi_log10_centered_G = bins_lpi_log10_G[0:-1]+0.5*(bins_lpi_log10_G[1:]-bins_lpi_log10_G[0:-1])
bins_lpi_centered_G = 10**bins_lpi_log10_centered_G

bins_LTG3_log10_centered_G = bins_LTG3_log10_G[0:-1]+0.5*(bins_LTG3_log10_G[1:]-bins_LTG3_log10_G[0:-1])
bins_LTG3_centered_G = 10**bins_LTG3_log10_centered_G

bins_PR92W_log10_centered_G = bins_PR92W_log10_G[0:-1]+0.5*(bins_PR92W_log10_G[1:]-bins_PR92W_log10_G[0:-1])
bins_PR92W_centered_G = 10**bins_PR92W_log10_centered_G

# DOMAIN 2
print('frequency d02')
n, bins_log10, patches = plt.hist(np.log10(L[L!=0]),bins = my_bins)
plt.close()
n_CAPExP_R_d02_T, bins_CAPExP_R_log10_d02_T, patches_CAPExP_R_d02_T = plt.hist(np.log10(CAPExP_R_adj_sorted_d02_T[CAPExP_R_adj_sorted_d02_T!=0]),bins=my_bins)
plt.close()
n_lpi_d02_T, bins_lpi_log10_d02_T, patches_lpi_d02_T = plt.hist(np.log10(LPI_adj_sorted_d02_T[LPI_adj_sorted_d02_T!=0]),bins=my_bins)
plt.close()
n_LTG3_d02_T, bins_LTG3_log10_d02_T, patches_LTG3_d02_T = plt.hist(np.log10(LTG3_adj_sorted_d02_T[LTG3_adj_sorted_d02_T!=0]),bins=my_bins)
plt.close()
n_PR92W_d02_T, bins_PR92W_log10_d02_T, patches_PR92W_d02_T = plt.hist(np.log10(PR92W_adj_sorted_d02_T[PR92W_adj_sorted_d02_T!=0]),bins=my_bins)
plt.close()
n_CAPExP_R_d02_G, bins_CAPExP_R_log10_d02_G, patches_CAPExP_R_d02_G = plt.hist(np.log10(CAPExP_R_adj_sorted_d02_G[CAPExP_R_adj_sorted_d02_G!=0]),bins=my_bins)
plt.close()
n_lpi_d02_G, bins_lpi_log10_d02_G, patches_lpi_d02_G = plt.hist(np.log10(LPI_adj_sorted_d02_G[LPI_adj_sorted_d02_G!=0]),bins=my_bins)
plt.close()
n_LTG3_d02_G, bins_LTG3_log10_d02_G, patches_LTG3_d02_G = plt.hist(np.log10(LTG3_adj_sorted_d02_G[LTG3_adj_sorted_d02_G!=0]),bins=my_bins)
plt.close()
n_PR92W_d02_G, bins_PR92W_log10_d02_G, patches_PR92W_d02_G = plt.hist(np.log10(PR92W_adj_sorted_d02_G[PR92W_adj_sorted_d02_G!=0]),bins=my_bins)
plt.close()

# set bins with 0 to the left to np.nan for a line plot, lines should just end to the left when there are no values
print('0 to nan d02')
for i in range(len(n)):
    if n[i]==0:
        n[i] = np.nan
    else:
        break
for i in range(len(n)):
    if n_CAPExP_R_d02_T[i] == 0:
        n_CAPExP_R_d02_T[i] = np.nan
    else:
        break
for i in range(len(n)):
    if n_lpi_d02_T[i] == 0:
        n_lpi_d02_T[i] = np.nan
    else:
        break
for i in range(len(n)):
    if n_LTG3_d02_T[i] == 0:
        n_LTG3_d02_T[i] = np.nan
    else:
        break
for i in range(len(n)):
    if n_PR92W_d02_T[i] == 0:
        n_PR92W_d02_T[i] = np.nan
    else:
        break
for i in range(len(n)):
    if n_CAPExP_R_d02_G[i] == 0:
        n_CAPExP_R_d02_G[i] = np.nan
    else:
        break
for i in range(len(n)):
    if n_lpi_d02_G[i] == 0:
        n_lpi_d02_G[i] = np.nan
    else:
        break
for i in range(len(n)):
    if n_LTG3_d02_G[i] == 0:
        n_LTG3_d02_G[i] = np.nan
    else:
        break
for i in range(len(n)):
    if n_PR92W_d02_G[i] == 0:
        n_PR92W_d02_G[i] = np.nan
    else:
        break

print('avg bins d02')
bins_CAPExP_R_log10_centered_d02_T = bins_CAPExP_R_log10_d02_T[0:-1]+0.5*(bins_CAPExP_R_log10_d02_T[1:]-bins_CAPExP_R_log10_d02_T[0:-1])
bins_CAPExP_R_centered_d02_T = 10**bins_CAPExP_R_log10_centered_d02_T

bins_lpi_log10_centered_d02_T = bins_lpi_log10_d02_T[0:-1]+0.5*(bins_lpi_log10_d02_T[1:]-bins_lpi_log10_d02_T[0:-1])
bins_lpi_centered_d02_T = 10**bins_lpi_log10_centered_d02_T

bins_LTG3_log10_centered_d02_T = bins_LTG3_log10_d02_T[0:-1]+0.5*(bins_LTG3_log10_d02_T[1:]-bins_LTG3_log10_d02_T[0:-1])
bins_LTG3_centered_d02_T = 10**bins_LTG3_log10_centered_d02_T

bins_PR92W_log10_centered_d02_T = bins_PR92W_log10_d02_T[0:-1]+0.5*(bins_PR92W_log10_d02_T[1:]-bins_PR92W_log10_d02_T[0:-1])
bins_PR92W_centered_d02_T = 10**bins_PR92W_log10_centered_d02_T

bins_CAPExP_R_log10_centered_d02_G = bins_CAPExP_R_log10_d02_G[0:-1]+0.5*(bins_CAPExP_R_log10_d02_G[1:]-bins_CAPExP_R_log10_d02_G[0:-1])
bins_CAPExP_R_centered_d02_G = 10**bins_CAPExP_R_log10_centered_d02_G

bins_lpi_log10_centered_d02_G = bins_lpi_log10_d02_G[0:-1]+0.5*(bins_lpi_log10_d02_G[1:]-bins_lpi_log10_d02_G[0:-1])
bins_lpi_centered_d02_G = 10**bins_lpi_log10_centered_d02_G

bins_LTG3_log10_centered_d02_G = bins_LTG3_log10_d02_G[0:-1]+0.5*(bins_LTG3_log10_d02_G[1:]-bins_LTG3_log10_d02_G[0:-1])
bins_LTG3_centered_d02_G = 10**bins_LTG3_log10_centered_d02_G

bins_PR92W_log10_centered_d02_G = bins_PR92W_log10_d02_G[0:-1]+0.5*(bins_PR92W_log10_d02_G[1:]-bins_PR92W_log10_d02_G[0:-1])
bins_PR92W_centered_d02_G = 10**bins_PR92W_log10_centered_d02_G

print('drop 0')
n[n==0] = np.nan
n_lpi_T[n_lpi_T==0] = np.nan
n_LTG3_T[n_LTG3_T==0] = np.nan
n_PR92W_T[n_PR92W_T==0] = np.nan
n_CAPExP_R_T[n_CAPExP_R_T == 0] = np.nan

n_lpi_d02_T[n_lpi_d02_T==0] = np.nan
n_LTG3_d02_T[n_LTG3_d02_T==0] = np.nan
n_PR92W_d02_T[n_PR92W_d02_T==0] = np.nan
n_CAPExP_R_d02_T[n_CAPExP_R_d02_T == 0] = np.nan

n_lpi_G[n_lpi_G==0] = np.nan
n_LTG3_G[n_LTG3_G==0] = np.nan
n_PR92W_G[n_PR92W_G==0] = np.nan
n_CAPExP_R_G[n_CAPExP_R_G == 0] = np.nan

n_lpi_d02_G[n_lpi_d02_G==0] = np.nan
n_LTG3_d02_G[n_LTG3_d02_G==0] = np.nan
n_PR92W_d02_G[n_PR92W_d02_G==0] = np.nan
n_CAPExP_R_d02_G[n_CAPExP_R_d02_G == 0] = np.nan
# ---------------------------------------------------------------------------------------------
# CONVERT TO PDF
# ---------------------------------------------------------------------------------------------
# PDF
print('PDF')
pdf_n = n/np.nansum(n)
pdf_lpi_T = n_lpi_T/np.nansum(n_lpi_T)
pdf_LTG3_T = n_LTG3_T/np.nansum(n_LTG3_T)
pdf_PR92W_T = n_PR92W_T/np.nansum(n_PR92W_T)
pdf_CAPExP_R_T = n_CAPExP_R_T/np.nansum(n_CAPExP_R_T)

pdf_lpi_d02_T = n_lpi_d02_T/np.nansum(n_lpi_d02_T)
pdf_LTG3_d02_T = n_LTG3_d02_T/np.nansum(n_LTG3_d02_T)
pdf_PR92W_d02_T = n_PR92W_d02_T/np.nansum(n_PR92W_d02_T)
pdf_CAPExP_R_d02_T = n_CAPExP_R_d02_T/np.nansum(n_CAPExP_R_d02_T)

pdf_lpi_G = n_lpi_G/np.nansum(n_lpi_G)
pdf_LTG3_G = n_LTG3_G/np.nansum(n_LTG3_G)
pdf_PR92W_G = n_PR92W_G/np.nansum(n_PR92W_G)
pdf_CAPExP_R_G = n_CAPExP_R_G/np.nansum(n_CAPExP_R_G)

pdf_lpi_d02_G = n_lpi_d02_G/np.nansum(n_lpi_d02_G)
pdf_LTG3_d02_G = n_LTG3_d02_G/np.nansum(n_LTG3_d02_G)
pdf_PR92W_d02_G = n_PR92W_d02_G/np.nansum(n_PR92W_d02_G)
pdf_CAPExP_R_d02_G = n_CAPExP_R_d02_G/np.nansum(n_CAPExP_R_d02_G)

PSS_LPI_T =  np.nansum(np.minimum(pdf_lpi_T, pdf_n))
PSS_LTG3_T =  np.nansum(np.minimum(pdf_LTG3_T, pdf_n))
PSS_PR92W_T =  np.nansum(np.minimum(pdf_PR92W_T, pdf_n))
PSS_CAPExP_R_T =  np.nansum(np.minimum(pdf_CAPExP_R_T, pdf_n))

PSS_LPI_d02_T =  np.nansum(np.minimum(pdf_lpi_d02_T, pdf_n))
PSS_LTG3_d02_T =  np.nansum(np.minimum(pdf_LTG3_d02_T, pdf_n))
PSS_PR92W_d02_T =  np.nansum(np.minimum(pdf_PR92W_d02_T, pdf_n))
PSS_CAPExP_R_d02_T =  np.nansum(np.minimum(pdf_CAPExP_R_d02_T, pdf_n))

PSS_LPI_G =  np.nansum(np.minimum(pdf_lpi_G, pdf_n))
PSS_LTG3_G =  np.nansum(np.minimum(pdf_LTG3_G, pdf_n))
PSS_PR92W_G =  np.nansum(np.minimum(pdf_PR92W_G, pdf_n))
PSS_CAPExP_R_G =  np.nansum(np.minimum(pdf_CAPExP_R_G, pdf_n))

PSS_LPI_d02_G =  np.nansum(np.minimum(pdf_lpi_d02_G, pdf_n))
PSS_LTG3_d02_G =  np.nansum(np.minimum(pdf_LTG3_d02_G, pdf_n))
PSS_PR92W_d02_G =  np.nansum(np.minimum(pdf_PR92W_d02_G, pdf_n))
PSS_CAPExP_R_d02_G =  np.nansum(np.minimum(pdf_CAPExP_R_d02_G, pdf_n))

print('PSS_LPI_d01_T = ' + str(PSS_LPI_T))
print('PSS_LTG3_d01_T = ' + str(PSS_LTG3_T))
print('PSS_PR92W_d01_T = ' + str(PSS_PR92W_T))
print('PSS_CAPExP_R_d01_T = ' + str(PSS_CAPExP_R_T))

print('PSS_LPI_d02_T = ' + str(PSS_LPI_d02_T))
print('PSS_LTG3_d02_T = ' + str(PSS_LTG3_d02_T))
print('PSS_PR92W_d02_T = ' + str(PSS_PR92W_d02_T))
print('PSS_CAPExP_R_d02_T = ' + str(PSS_CAPExP_R_d02_T))

print('PSS_LPI_d01_G = ' + str(PSS_LPI_G))
print('PSS_LTG3_d01_G = ' + str(PSS_LTG3_G))
print('PSS_PR92W_d01_G = ' + str(PSS_PR92W_G))
print('PSS_CAPExP_R_d01_G = ' + str(PSS_CAPExP_R_G))

print('PSS_LPI_d02_G = ' + str(PSS_LPI_d02_G))
print('PSS_LTG3_d02_G = ' + str(PSS_LTG3_d02_G))
print('PSS_PR92W_d02_G = ' + str(PSS_PR92W_d02_G))
print('PSS_CAPExP_R_d02_G = ' + str(PSS_CAPExP_R_d02_G))
