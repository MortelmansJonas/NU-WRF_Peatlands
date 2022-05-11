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
# Thompson
ds_d01_T = Dataset('/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/NU-WRF/wrfout_nc_files/d01_all_Thompson.nc','r')
ds_d02_T = Dataset('/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/NU-WRF/wrfout_nc_files/domain2_at_domain1_all_Thompson.nc','r')
# Goddard
ds_d01_G = Dataset('/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/NU-WRF/wrfout_nc_files/d01_all.nc','r')
ds_d02_G = Dataset('/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/NU-WRF/wrfout_nc_files/domain2_at_domain1_all.nc','r')
# Observations
ds_obs = Dataset('/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/NU-WRF/wrfout_nc_files/CLDN_at_domain1_all_v4.nc')
L = ds_obs['Flashdensity_CC'][:,:,:].data.flatten() + ds_obs['Flashdensity_CG'][:,:,:].data.flatten()

# Extract latitudes and longitudes
lat = ds_d01_T['lat'][:]
lon = ds_d01_T['lon'][:]

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

CAPExP_CSI_d01_T = ds_d01_T['CAPExP_CSI'][:]
CAPExP_CSI_d01_T[CAPExP_CSI_d01_T == np.nanmin(CAPExP_CSI_d01_T)] = 0

LPI_d02_T = ds_d02_T['LPI'][:]
LPI_d02_T[LPI_d02_T == np.nanmin(LPI_d02_T)] = 0

LTG3_d02_T = ds_d02_T['LTG3'][:]
LTG3_d02_T[LTG3_d02_T == np.nanmin(LTG3_d02_T)] = 0

PR92W_d02_T = ds_d02_T['PR92W'][:]
PR92W_d02_T[PR92W_d02_T == np.nanmin(PR92W_d02_T)] = 0

CAPExP_R_d02_T = ds_d02_T['CAPExP_R'][:]
CAPExP_R_d02_T[CAPExP_R_d02_T == np.nanmin(CAPExP_R_d02_T)] = 0

CAPExP_CSI_d02_T = ds_d02_T['CAPExP_CSI'][:]
CAPExP_CSI_d02_T[CAPExP_CSI_d02_T == np.nanmin(CAPExP_CSI_d02_T)] = 0

LPI_d01_G = ds_d01_G['LPI'][:]
LPI_d01_G[LPI_d01_G == np.nanmin(LPI_d01_G)] = 0

LTG3_d01_G = ds_d01_G['LTG3'][:]
LTG3_d01_G[LTG3_d01_G == np.nanmin(LTG3_d01_G)] = 0

PR92W_d01_G = ds_d01_G['PR92W'][:]
PR92W_d01_G[PR92W_d01_G == np.nanmin(PR92W_d01_G)] = 0

CAPExP_R_d01_G = ds_d01_G['CAPExP_R'][:]
CAPExP_R_d01_G[CAPExP_R_d01_G == np.nanmin(CAPExP_R_d01_G)] = 0

CAPExP_CSI_d01_G = ds_d01_G['CAPExP_CSI'][:]
CAPExP_CSI_d01_G[CAPExP_CSI_d01_G == np.nanmin(CAPExP_CSI_d01_G)] = 0

LPI_d02_G = ds_d02_G['LPI'][:]
LPI_d02_G[LPI_d02_G == np.nanmin(LPI_d02_G)] = 0

LTG3_d02_G = ds_d02_G['LTG3'][:]
LTG3_d02_G[LTG3_d02_G == np.nanmin(LTG3_d02_G)] = 0

PR92W_d02_G = ds_d02_G['PR92W'][:]
PR92W_d02_G[PR92W_d02_G == np.nanmin(PR92W_d02_G)] = 0

CAPExP_R_d02_G = ds_d02_G['CAPExP_R'][:]
CAPExP_R_d02_G[CAPExP_R_d02_G == np.nanmin(CAPExP_R_d02_G)] = 0

CAPExP_CSI_d02_G = ds_d02_G['CAPExP_CSI'][:]
CAPExP_CSI_d02_G[CAPExP_CSI_d02_G == np.nanmin(CAPExP_CSI_d02_G)] = 0

# -------------------------------------------------------------------------------------------------------------
# 'Second, a simple linear model is built that relates the observed flash rates to the parameterization output'
# -------------------------------------------------------------------------------------------------------------
# First determine cutoff value 'c', done similarly to the sorting in 'Calibration_2.0.py'
Obs = np.sort(L[L!=0])[::-1]

# Sort model and take n highest values (with n = number of observations)
c = Obs.shape[0]
Tot_flashes = np.nansum(Obs)

# THOM - 9 km
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

CAPExP_CSI_sorted_d01_T = np.sort(CAPExP_CSI_d01_T.flatten())[::-1]
CAPExP_CSI_c_d01_T = CAPExP_CSI_sorted_d01_T[0:c]
CAPExP_CSI_d01_T_new = np.where(CAPExP_CSI_d01_T>=CAPExP_CSI_sorted_d01_T[c],CAPExP_CSI_d01_T, np.nan)
print('cutoff value CAPExP_CSI_d01_T = ', CAPExP_CSI_sorted_d01_T[c])

CAPExP_R_sorted_d01_T = np.sort(CAPExP_R_d01_T.flatten())[::-1]
CAPExP_R_c_d01_T = CAPExP_R_sorted_d01_T[0:c]
CAPExP_R_d01_T_new = np.where(CAPExP_R_d01_T>=CAPExP_R_sorted_d01_T[c],CAPExP_R_d01_T, np.nan)
print('cutoff value CAPExP_R_d01_T = ', CAPExP_R_sorted_d01_T[c])

# THOM - 3 km
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

CAPExP_CSI_sorted_d02_T = np.sort(CAPExP_CSI_d02_T.flatten())[::-1]
CAPExP_CSI_c_d02_T = CAPExP_CSI_sorted_d02_T[0:c]
CAPExP_CSI_d02_T_new = np.where(CAPExP_CSI_d02_T>=CAPExP_CSI_sorted_d02_T[c],CAPExP_CSI_d02_T, np.nan)
print('cutoff value CAPExP_CSI_d02_T = ', CAPExP_CSI_sorted_d02_T[c])

CAPExP_R_sorted_d02_T = np.sort(CAPExP_R_d02_T.flatten())[::-1]
CAPExP_R_c_d02_T = CAPExP_R_sorted_d02_T[0:c]
CAPExP_R_d02_T_new = np.where(CAPExP_R_d02_T>=CAPExP_R_sorted_d02_T[c],CAPExP_R_d02_T, np.nan)
print('cutoff value CAPExP_R_d02_T = ', CAPExP_R_sorted_d02_T[c])

# G4ICE - 9 km
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

CAPExP_CSI_sorted_d01_G = np.sort(CAPExP_CSI_d01_G.flatten())[::-1]
CAPExP_CSI_c_d01_G = CAPExP_CSI_sorted_d01_G[0:c]
CAPExP_CSI_d01_G_new = np.where(CAPExP_CSI_d01_G>=CAPExP_CSI_sorted_d01_G[c],CAPExP_CSI_d01_G, np.nan)
print('cutoff value CAPExP_CSI_d01_G = ', CAPExP_CSI_sorted_d01_G[c])

CAPExP_R_sorted_d01_G = np.sort(CAPExP_R_d01_G.flatten())[::-1]
CAPExP_R_c_d01_G = CAPExP_R_sorted_d01_G[0:c]
CAPExP_R_d01_G_new = np.where(CAPExP_R_d01_G>=CAPExP_R_sorted_d01_G[c],CAPExP_R_d01_G, np.nan)
print('cutoff value CAPExP_R_d01_G = ', CAPExP_R_sorted_d01_G[c])

# G4ICE - 3 km
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

CAPExP_CSI_sorted_d02_G = np.sort(CAPExP_CSI_d02_G.flatten())[::-1]
CAPExP_CSI_c_d02_G = CAPExP_CSI_sorted_d02_G[0:c]
CAPExP_CSI_d02_G_new = np.where(CAPExP_CSI_d02_G>=CAPExP_CSI_sorted_d02_G[c],CAPExP_CSI_d02_G, np.nan)
print('cutoff value CAPExP_CSI_d02_G = ', CAPExP_CSI_sorted_d02_G[c])

CAPExP_R_sorted_d02_G = np.sort(CAPExP_R_d02_G.flatten())[::-1]
CAPExP_R_c_d02_G = CAPExP_R_sorted_d02_G[0:c]
CAPExP_R_d02_G_new = np.where(CAPExP_R_d02_G>=CAPExP_R_sorted_d02_G[c],CAPExP_R_d02_G, np.nan)
print('cutoff value CAPExP_R_d02_G = ', CAPExP_R_sorted_d02_G[c])
# -------------------------------------------------------------------------------------------------------------
# Linear regression: X_adj = aX+b if X >= c
# -------------------------------------------------------------------------------------------------------------
# THOM - 9 km
reg_LPI_d01_T = lr().fit(LPI_c_d01_T.reshape((-1,1)), Obs)
LPI_d01_T_adj =np.add(reg_LPI_d01_T.intercept_, np.multiply(reg_LPI_d01_T.coef_,LPI_d01_T_new))
LPI_d01_T_adj[np.isnan(LPI_d01_T_adj)] = 0

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

reg_CAPExP_CSI_d01_T = lr().fit(CAPExP_CSI_c_d01_T.reshape((-1,1)), Obs)
CAPExP_CSI_d01_T_adj =np.add(reg_CAPExP_CSI_d01_T.intercept_, np.multiply(reg_CAPExP_CSI_d01_T.coef_,CAPExP_CSI_d01_T_new))
CAPExP_CSI_d01_T_adj[np.isnan(CAPExP_CSI_d01_T_adj)] = 0

reg_CAPExP_R_d01_T = lr().fit(CAPExP_R_c_d01_T.reshape((-1,1)), Obs)
CAPExP_R_d01_T_adj =np.add(reg_CAPExP_R_d01_T.intercept_, np.multiply(reg_CAPExP_R_d01_T.coef_,CAPExP_R_d01_T_new))
CAPExP_R_d01_T_adj[np.isnan(CAPExP_R_d01_T_adj)] = 0

# THOM - 3 km
reg_LPI_d02_T = lr().fit(LPI_c_d02_T.reshape((-1,1)), Obs)
LPI_d02_T_adj =np.add(reg_LPI_d02_T.intercept_, np.multiply(reg_LPI_d02_T.coef_,LPI_d02_T_new))
LPI_d02_T_adj[np.isnan(LPI_d02_T_adj)] = 0

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

reg_CAPExP_CSI_d02_T = lr().fit(CAPExP_CSI_c_d02_T.reshape((-1,1)), Obs)
CAPExP_CSI_d02_T_adj =np.add(reg_CAPExP_CSI_d02_T.intercept_, np.multiply(reg_CAPExP_CSI_d02_T.coef_,CAPExP_CSI_d02_T_new))
CAPExP_CSI_d02_T_adj[np.isnan(CAPExP_CSI_d02_T_adj)] = 0

reg_CAPExP_R_d02_T = lr().fit(CAPExP_R_c_d02_T.reshape((-1,1)), Obs)
CAPExP_R_d02_T_adj =np.add(reg_CAPExP_R_d02_T.intercept_, np.multiply(reg_CAPExP_R_d02_T.coef_,CAPExP_R_d02_T_new))
CAPExP_R_d02_T_adj[np.isnan(CAPExP_R_d02_T_adj)] = 0

# G4ICE - 9 km
reg_LPI_d01_G = lr().fit(LPI_c_d01_G.reshape((-1,1)), Obs)
LPI_d01_G_adj =np.add(reg_LPI_d01_G.intercept_, np.multiply(reg_LPI_d01_G.coef_,LPI_d01_G_new))
LPI_d01_G_adj[np.isnan(LPI_d01_G_adj)] = 0

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

reg_CAPExP_CSI_d01_G = lr().fit(CAPExP_CSI_c_d01_G.reshape((-1,1)), Obs)
CAPExP_CSI_d01_G_adj =np.add(reg_CAPExP_CSI_d01_G.intercept_, np.multiply(reg_CAPExP_CSI_d01_G.coef_,CAPExP_CSI_d01_G_new))
CAPExP_CSI_d01_G_adj[np.isnan(CAPExP_CSI_d01_G_adj)] = 0

reg_CAPExP_R_d01_G = lr().fit(CAPExP_R_c_d01_G.reshape((-1,1)), Obs)
CAPExP_R_d01_G_adj =np.add(reg_CAPExP_R_d01_G.intercept_, np.multiply(reg_CAPExP_R_d01_G.coef_,CAPExP_R_d01_G_new))
CAPExP_R_d01_G_adj[np.isnan(CAPExP_R_d01_G_adj)] = 0

# G4ICE - 3 km
reg_LPI_d02_G = lr().fit(LPI_c_d02_G.reshape((-1,1)), Obs)
LPI_d02_G_adj =np.add(reg_LPI_d02_G.intercept_, np.multiply(reg_LPI_d02_G.coef_,LPI_d02_G_new))
LPI_d02_G_adj[np.isnan(LPI_d02_G_adj)] = 0

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

reg_CAPExP_CSI_d02_G = lr().fit(CAPExP_CSI_c_d02_G.reshape((-1,1)), Obs)
CAPExP_CSI_d02_G_adj =np.add(reg_CAPExP_CSI_d02_G.intercept_, np.multiply(reg_CAPExP_CSI_d02_G.coef_,CAPExP_CSI_d02_G_new))
CAPExP_CSI_d02_G_adj[np.isnan(CAPExP_CSI_d02_G_adj)] = 0

reg_CAPExP_R_d02_G = lr().fit(CAPExP_R_c_d02_G.reshape((-1,1)), Obs)
CAPExP_R_d02_G_adj =np.add(reg_CAPExP_R_d02_G.intercept_, np.multiply(reg_CAPExP_R_d02_G.coef_,CAPExP_R_d02_G_new))
CAPExP_R_d02_G_adj[np.isnan(CAPExP_R_d02_G_adj)] = 0

# sort the newly adjusted data descending
LPI_adj_sorted_d01_T = np.sort(LPI_d01_T_adj.flatten())[::-1]
LTG3_adj_sorted_d01_T = np.sort(LTG3_d01_T_adj.flatten())[::-1]
PR92W_adj_sorted_d01_T = np.sort(PR92W_d01_T_adj.flatten())[::-1]
CAPExP_CSI_adj_sorted_d01_T = np.sort(CAPExP_CSI_d01_T_adj.flatten())[::-1]
CAPExP_R_adj_sorted_d01_T = np.sort(CAPExP_R_d01_T_adj.flatten())[::-1]
LPI_adj_sorted_d02_T = np.sort(LPI_d02_T_adj.flatten())[::-1]
LTG3_adj_sorted_d02_T = np.sort(LTG3_d02_T_adj.flatten())[::-1]
PR92W_adj_sorted_d02_T = np.sort(PR92W_d02_T_adj.flatten())[::-1]
CAPExP_CSI_adj_sorted_d02_T = np.sort(CAPExP_CSI_d02_T_adj.flatten())[::-1]
CAPExP_R_adj_sorted_d02_T = np.sort(CAPExP_R_d02_T_adj.flatten())[::-1]

LPI_adj_sorted_d01_G = np.sort(LPI_d01_G_adj.flatten())[::-1]
LTG3_adj_sorted_d01_G = np.sort(LTG3_d01_G_adj.flatten())[::-1]
PR92W_adj_sorted_d01_G = np.sort(PR92W_d01_G_adj.flatten())[::-1]
CAPExP_CSI_adj_sorted_d01_G = np.sort(CAPExP_CSI_d01_G_adj.flatten())[::-1]
CAPExP_R_adj_sorted_d01_G = np.sort(CAPExP_R_d01_G_adj.flatten())[::-1]
LPI_adj_sorted_d02_G = np.sort(LPI_d02_G_adj.flatten())[::-1]
LTG3_adj_sorted_d02_G = np.sort(LTG3_d02_G_adj.flatten())[::-1]
PR92W_adj_sorted_d02_G = np.sort(PR92W_d02_G_adj.flatten())[::-1]
CAPExP_CSI_adj_sorted_d02_G = np.sort(CAPExP_CSI_d02_G_adj.flatten())[::-1]
CAPExP_R_adj_sorted_d02_G = np.sort(CAPExP_R_d02_G_adj.flatten())[::-1]

# Make the total amount of modeled flashes equal to total amount of observed flashes by taking the N highest flash densities.
sum_lpi_d01_T = 0
sum_ltg3_d01_T = 0
sum_pr92w_d01_T = 0
sum_capexp_csi_d01_T = 0
sum_capexp_r_d01_T = 0
sum_lpi_d02_T = 0
sum_ltg3_d02_T = 0
sum_pr92w_d02_T = 0
sum_capexp_csi_d02_T = 0
sum_capexp_r_d02_T = 0
sum_lpi_d01_G = 0
sum_ltg3_d01_G = 0
sum_pr92w_d01_G = 0
sum_capexp_csi_d01_G = 0
sum_capexp_r_d01_G = 0
sum_lpi_d02_G = 0
sum_ltg3_d02_G = 0
sum_pr92w_d02_G = 0
sum_capexp_csi_d02_G = 0
sum_capexp_r_d02_G = 0

# THOM - 9 km
i = 0
while sum_lpi_d01_T < Tot_flashes:
    sum_lpi_d01_T = sum_lpi_d01_T + LPI_adj_sorted_d01_T[i]
    i += 1
else:
    LPI_adj_sorted_d01_T[i:] = 0
    LPI_cutoff_d01_T = LPI_adj_sorted_d01_T[i]
    LPI_adj_sorted_d01_T.reshape((13248,63,109))
    print(i)

i = 0
while (sum_ltg3_d01_T < Tot_flashes) & (i < len(LTG3_adj_sorted_d01_T)):
    sum_ltg3_d01_T = sum_ltg3_d01_T + LTG3_adj_sorted_d01_T[i]
    i += 1
else:
    LTG3_adj_sorted_d01_T[i:] = 0
    LTG3_cutoff_d01_T = LTG3_adj_sorted_d01_T[i]
    LTG3_adj_sorted_d01_T.reshape((13248, 63, 109))
    print(i)

i=0
while sum_pr92w_d01_T < Tot_flashes:
    sum_pr92w_d01_T = sum_pr92w_d01_T + PR92W_adj_sorted_d01_T[i]
    i += 1
else:
    PR92W_adj_sorted_d01_T[i:] = 0
    PR92W_cutoff_d01_T = PR92W_adj_sorted_d01_T[i]
    PR92W_adj_sorted_d01_T.reshape((13248, 63, 109))
    print(i)

i=0
while sum_capexp_csi_d01_T < Tot_flashes:
    sum_capexp_csi_d01_T = sum_capexp_csi_d01_T + CAPExP_CSI_adj_sorted_d01_T[i]
    i += 1
else:
    CAPExP_CSI_adj_sorted_d01_T[i:] = 0
    CAPExP_CSI_cutoff_d01_T = CAPExP_CSI_adj_sorted_d01_T[i]
    CAPExP_CSI_adj_sorted_d01_T.reshape((13248, 63, 109))
    print(i)

i=0
while sum_capexp_r_d01_T < Tot_flashes:
    sum_capexp_r_d01_T = sum_capexp_r_d01_T + CAPExP_R_adj_sorted_d01_T[i]
    i += 1
else:
    CAPExP_R_adj_sorted_d01_T[i:] = 0
    CAPExP_R_cutoff_d01_T = CAPExP_R_adj_sorted_d01_T[i]
    CAPExP_R_adj_sorted_d01_T.reshape((13248, 63, 109))
    print(i)

# THOM - 3 km
i=0
while sum_lpi_d02_T < Tot_flashes:
    sum_lpi_d02_T = sum_lpi_d02_T + LPI_adj_sorted_d02_T[i]
    i += 1
else:
    LPI_adj_sorted_d02_T[i:] = 0
    LPI_cutoff_d02_T = LPI_adj_sorted_d02_T[i]
    LPI_adj_sorted_d02_T.reshape((13248, 63, 109))
    print(i)

i=0
while (sum_ltg3_d02_T < Tot_flashes) & (i < len(LTG3_adj_sorted_d02_T)):
    sum_ltg3_d02_T = sum_ltg3_d02_T + LTG3_adj_sorted_d02_T[i]
    i += 1
else:
    LTG3_adj_sorted_d02_T[i:] = 0
    LTG3_cutoff_d02_T = LTG3_adj_sorted_d02_T[i]
    LTG3_adj_sorted_d02_T.reshape((13248, 63, 109))
    print(i)

i=0
while sum_pr92w_d02_T < Tot_flashes:
    sum_pr92w_d02_T = sum_pr92w_d02_T + PR92W_adj_sorted_d02_T[i]
    i += 1
else:
    PR92W_adj_sorted_d02_T[i:] = 0
    PR92W_cutoff_d02_T = PR92W_adj_sorted_d02_T[i]
    PR92W_adj_sorted_d02_T.reshape((13248, 63, 109))
    print(i)

i=0
while sum_capexp_csi_d02_T < Tot_flashes:
    sum_capexp_csi_d02_T = sum_capexp_csi_d02_T + CAPExP_CSI_adj_sorted_d02_T[i]
    i += 1
else:
    CAPExP_CSI_adj_sorted_d02_T[i:] = 0
    CAPExP_CSI_cutoff_d02_T = CAPExP_CSI_adj_sorted_d02_T[i]
    CAPExP_CSI_adj_sorted_d02_T.reshape((13248, 63, 109))
    print(i)

i=0
while sum_capexp_r_d02_T < Tot_flashes:
    sum_capexp_r_d02_T = sum_capexp_r_d02_T + CAPExP_R_adj_sorted_d02_T[i]
    i += 1
else:
    CAPExP_R_adj_sorted_d02_T[i:] = 0
    CAPExP_R_cutoff_d02_T = CAPExP_R_adj_sorted_d02_T[i]
    CAPExP_R_adj_sorted_d02_T.reshape((13248, 63, 109))
    print(i)

# G4ICE - 9 km
i=0
while sum_lpi_d01_G < Tot_flashes:
    sum_lpi_d01_G = sum_lpi_d01_G + LPI_adj_sorted_d01_G[i]
    i += 1
else:
    LPI_adj_sorted_d01_G[i:] = 0
    LPI_cutoff_d01_G = LPI_adj_sorted_d01_G[i]
    LPI_adj_sorted_d01_G.reshape((13248,63,109))
    print(i)

i = 0
while (sum_ltg3_d01_G < Tot_flashes) & (i < len(LTG3_adj_sorted_d01_G)):
    sum_ltg3_d01_G = sum_ltg3_d01_G + LTG3_adj_sorted_d01_G[i]
    i += 1
else:
    LTG3_adj_sorted_d01_G[i:] = 0
    LTG3_cutoff_d01_G = LTG3_adj_sorted_d01_G[i]
    LTG3_adj_sorted_d01_G.reshape((13248, 63, 109))
    print(i)

i=0
while sum_pr92w_d01_G < Tot_flashes:
    sum_pr92w_d01_G = sum_pr92w_d01_G + PR92W_adj_sorted_d01_G[i]
    i += 1
else:
    PR92W_adj_sorted_d01_G[i:] = 0
    PR92W_cutoff_d01_G = PR92W_adj_sorted_d01_G[i]
    PR92W_adj_sorted_d01_G.reshape((13248, 63, 109))
    print(i)

i=0
while sum_capexp_csi_d01_G < Tot_flashes:
    sum_capexp_csi_d01_G = sum_capexp_csi_d01_G + CAPExP_CSI_adj_sorted_d01_G[i]
    i += 1
else:
    CAPExP_CSI_adj_sorted_d01_G[i:] = 0
    CAPExP_CSI_cutoff_d01_G = CAPExP_CSI_adj_sorted_d01_G[i]
    CAPExP_CSI_adj_sorted_d01_G.reshape((13248, 63, 109))
    print(i)

i=0
while sum_capexp_r_d01_G < Tot_flashes:
    sum_capexp_r_d01_G = sum_capexp_r_d01_G + CAPExP_R_adj_sorted_d01_G[i]
    i += 1
else:
    CAPExP_R_adj_sorted_d01_G[i:] = 0
    CAPExP_R_cutoff_d01_G = CAPExP_R_adj_sorted_d01_G[i]
    CAPExP_R_adj_sorted_d01_G.reshape((13248, 63, 109))
    print(i)

# G4ICE - 3 km
i=0
while sum_lpi_d02_G < Tot_flashes:
    sum_lpi_d02_G = sum_lpi_d02_G + LPI_adj_sorted_d02_G[i]
    i += 1
else:
    LPI_adj_sorted_d02_G[i:] = 0
    LPI_cutoff_d02_G = LPI_adj_sorted_d02_G[i]
    LPI_adj_sorted_d02_G.reshape((13248, 63, 109))
    print(i)

i=0
while (sum_ltg3_d02_G < Tot_flashes) & (i < len(LTG3_adj_sorted_d02_G)):
    sum_ltg3_d02_G = sum_ltg3_d02_G + LTG3_adj_sorted_d02_G[i]
    i += 1
else:
    LTG3_adj_sorted_d02_G[i:] = 0
    # LTG3_cutoff_d02_G = LTG3_adj_sorted_d02_G[i]
    LTG3_adj_sorted_d02_G.reshape((13248, 63, 109))
    print(i)

i=0
while sum_pr92w_d02_G < Tot_flashes:
    sum_pr92w_d02_G = sum_pr92w_d02_G + PR92W_adj_sorted_d02_G[i]
    i += 1
else:
    PR92W_adj_sorted_d02_G[i:] = 0
    PR92W_cutoff_d02_G = PR92W_adj_sorted_d02_G[i]
    PR92W_adj_sorted_d02_G.reshape((13248, 63, 109))
    print(i)

i=0
while sum_capexp_csi_d02_G < Tot_flashes:
    sum_capexp_csi_d02_G = sum_capexp_csi_d02_G + CAPExP_CSI_adj_sorted_d02_G[i]
    i += 1
else:
    CAPExP_CSI_adj_sorted_d02_G[i:] = 0
    CAPExP_CSI_cutoff_d02_G = CAPExP_CSI_adj_sorted_d02_G[i]
    CAPExP_CSI_adj_sorted_d02_G.reshape((13248, 63, 109))
    print(i)

i=0
while sum_capexp_r_d02_G < Tot_flashes:
    sum_capexp_r_d02_G = sum_capexp_r_d02_G + CAPExP_R_adj_sorted_d02_G[i]
    i += 1
else:
    CAPExP_R_adj_sorted_d02_G[i:] = 0
    CAPExP_R_cutoff_d02_G = CAPExP_R_adj_sorted_d02_G[i]
    CAPExP_R_adj_sorted_d02_G.reshape((13248, 63, 109))
    print(i)

# ---------------------------------------------------------------------------------------------
# PUT CALIBRATED DATA IN NEW NETCDF FILE
# ---------------------------------------------------------------------------------------------
# Create .nc file for Thompson
ds = Dataset('/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/NU-WRF/wrfout_nc_files/data_calibrated_ax_d01_Thompson.nc', mode='w', format='NETCDF4')
ds.createDimension('time', 13248)
ds.createDimension('lat', 63)
ds.createDimension('lon', 109)
ds.createVariable('time','int', dimensions=('time',),zlib=True)
ds.variables['time'][:] = ds_d01_T['time'][:]
ds.createVariable('lat',ds_d01_T['lat'][:].dtype, dimensions=('lat','lon',),zlib=True)
ds.variables['lat'][:] = ds_d01_T['lat'][:]
ds.createVariable('lon',ds_d01_T['lon'][:].dtype, dimensions=('lat','lon',),zlib=True)
ds.variables['lon'][:] = ds_d01_T['lon'][:]
ds.createVariable('LPI_d01', 'f4', dimensions=('time','lat','lon',), zlib=True)
ds.createVariable('LTG3_d01', 'f4', dimensions=('time','lat','lon',), zlib=True)
ds.createVariable('PR92W_d01', 'f4', dimensions=('time','lat','lon',), zlib=True)
ds.createVariable('CAPExP_R_d01', 'f4', dimensions=('time','lat','lon',), zlib=True)
ds.createVariable('CAPExP_CSI_d01', 'f4', dimensions=('time','lat','lon',), zlib=True)
ds.createVariable('LPI_d02', 'f4', dimensions=('time','lat','lon',), zlib=True)
ds.createVariable('LTG3_d02', 'f4', dimensions=('time','lat','lon',), zlib=True)
ds.createVariable('PR92W_d02', 'f4', dimensions=('time','lat','lon',), zlib=True)
ds.createVariable('CAPExP_R_d02', 'f4', dimensions=('time','lat','lon',), zlib=True)
ds.createVariable('CAPExP_CSI_d02', 'f4', dimensions=('time','lat','lon',), zlib=True)
ds.createVariable('Obs', 'f4', dimensions=('time','lat','lon',), zlib=True)

# Fill it
ds['LPI_d01'][:] = LPI_d01_T_adj[:]
ds['LPI_d02'][:] = LPI_d02_T_adj[:]
ds['PR92W_d01'][:] = PR92W_d01_T_adj[:]
ds['PR92W_d02'][:] = PR92W_d02_T_adj[:]
ds['LTG3_d01'][:] = LTG3_d01_T_adj[:]
ds['LTG3_d02'][:] = LTG3_d02_T_adj[:]
ds['CAPExP_CSI_d01'][:] = CAPExP_CSI_d01_T_adj[:]
ds['CAPExP_CSI_d02'][:] = CAPExP_CSI_d02_T_adj[:]
ds['CAPExP_R_d01'][:] = CAPExP_R_d01_T_adj[:]
ds['CAPExP_R_d02'][:] = CAPExP_R_d02_T_adj[:]
ds['Obs'][:] = ds_obs['Flashdensity_CC'][:,:,:].data + ds_obs['Flashdensity_CG'][:,:,:].data
ds.close()

# Create .nc file for Goddard
ds = Dataset('/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/NU-WRF/wrfout_nc_files/data_calibrated_ax_d01.nc', mode='w', format='NETCDF4')
ds.createDimension('time', 13248)
ds.createDimension('lat', 63)
ds.createDimension('lon', 109)
ds.createVariable('time','int', dimensions=('time',),zlib=True)
ds.variables['time'][:] = ds_d01_G['time'][:]
ds.createVariable('lat',ds_d01_G['lat'][:].dtype, dimensions=('lat','lon',),zlib=True)
ds.variables['lat'][:] = ds_d01_G['lat'][:]
ds.createVariable('lon',ds_d01_G['lon'][:].dtype, dimensions=('lat','lon',),zlib=True)
ds.variables['lon'][:] = ds_d01_G['lon'][:]
ds.createVariable('LPI_d01', 'f4', dimensions=('time','lat','lon',), zlib=True)
ds.createVariable('LTG3_d01', 'f4', dimensions=('time','lat','lon',), zlib=True)
ds.createVariable('PR92W_d01', 'f4', dimensions=('time','lat','lon',), zlib=True)
ds.createVariable('CAPExP_R_d01', 'f4', dimensions=('time','lat','lon',), zlib=True)
ds.createVariable('CAPExP_CSI_d01', 'f4', dimensions=('time','lat','lon',), zlib=True)
ds.createVariable('LPI_d02', 'f4', dimensions=('time','lat','lon',), zlib=True)
ds.createVariable('LTG3_d02', 'f4', dimensions=('time','lat','lon',), zlib=True)
ds.createVariable('PR92W_d02', 'f4', dimensions=('time','lat','lon',), zlib=True)
ds.createVariable('CAPExP_R_d02', 'f4', dimensions=('time','lat','lon',), zlib=True)
ds.createVariable('CAPExP_CSI_d02', 'f4', dimensions=('time','lat','lon',), zlib=True)
ds.createVariable('Obs', 'f4', dimensions=('time','lat','lon',), zlib=True)

# Fill it
ds['LPI_d01'][:] = LPI_d01_G_adj[:]
ds['LPI_d02'][:] = LPI_d02_G_adj[:]
ds['PR92W_d01'][:] = PR92W_d01_G_adj[:]
ds['PR92W_d02'][:] = PR92W_d02_G_adj[:]
ds['LTG3_d01'][:] = LTG3_d01_G_adj[:]
ds['LTG3_d02'][:] = LTG3_d02_G_adj[:]
ds['CAPExP_CSI_d01'][:] = CAPExP_CSI_d01_G_adj[:]
ds['CAPExP_CSI_d02'][:] = CAPExP_CSI_d02_G_adj[:]
ds['CAPExP_R_d01'][:] = CAPExP_R_d01_G_adj[:]
ds['CAPExP_R_d02'][:] = CAPExP_R_d02_G_adj[:]
ds['Obs'][:] = ds_obs['Flashdensity_CC'][:,:,:].data + ds_obs['Flashdensity_CG'][:,:,:].data
ds.close()

# -------------------------------------------------------------------------------------------------------------
# FREQUENCY IFO HOURLY GRID FLASH DENSITY PLOTS
# -------------------------------------------------------------------------------------------------------------
# First set all 0 to nan
LPI_adj_sorted_d01_T[LPI_adj_sorted_d01_T == 0] = np.nan
LPI_adj_sorted_d02_T[LPI_adj_sorted_d02_T == 0] = np.nan
LTG3_adj_sorted_d01_T[LTG3_adj_sorted_d01_T == 0] = np.nan
LTG3_adj_sorted_d02_T[LTG3_adj_sorted_d02_T == 0] = np.nan
PR92W_adj_sorted_d01_T[PR92W_adj_sorted_d01_T == 0] = np.nan
PR92W_adj_sorted_d02_T[PR92W_adj_sorted_d02_T == 0] = np.nan
CAPExP_CSI_adj_sorted_d01_T[CAPExP_CSI_adj_sorted_d01_T == 0] = np.nan
CAPExP_CSI_adj_sorted_d02_T[CAPExP_CSI_adj_sorted_d02_T == 0] = np.nan
CAPExP_R_adj_sorted_d02_T[CAPExP_R_adj_sorted_d02_T == 0] = np.nan
CAPExP_R_adj_sorted_d01_T[CAPExP_R_adj_sorted_d01_T == 0] = np.nan
LPI_adj_sorted_d01_G[LPI_adj_sorted_d01_G == 0] = np.nan
LPI_adj_sorted_d02_G[LPI_adj_sorted_d02_G == 0] = np.nan
LTG3_adj_sorted_d01_G[LTG3_adj_sorted_d01_G == 0] = np.nan
LTG3_adj_sorted_d02_G[LTG3_adj_sorted_d02_G == 0] = np.nan
PR92W_adj_sorted_d01_G[PR92W_adj_sorted_d01_G == 0] = np.nan
PR92W_adj_sorted_d02_G[PR92W_adj_sorted_d02_G == 0] = np.nan
CAPExP_CSI_adj_sorted_d01_G[CAPExP_CSI_adj_sorted_d01_G == 0] = np.nan
CAPExP_CSI_adj_sorted_d02_G[CAPExP_CSI_adj_sorted_d02_G == 0] = np.nan
CAPExP_R_adj_sorted_d02_G[CAPExP_R_adj_sorted_d02_G == 0] = np.nan
CAPExP_R_adj_sorted_d01_G[CAPExP_R_adj_sorted_d01_G == 0] = np.nan

# specify bin edges
my_bins = np.linspace(-3,2,20)
# DOMAIN 1
n, bins_log10, patches = plt.hist(np.log10(L[L!=0]),bins = my_bins)
plt.close()
# THOM - 9 km
n_CAPExP_R_T, bins_CAPExP_R_log10_T, patches_CAPExP_R_T = plt.hist(np.log10(CAPExP_R_adj_sorted_d01_T[CAPExP_R_adj_sorted_d01_T!=0]),bins=my_bins)
plt.close()
n_CAPExP_CSI_T, bins_CAPExP_CSI_log10_T, patches_CAPExP_CSI_T = plt.hist(np.log10(CAPExP_CSI_adj_sorted_d01_T[CAPExP_CSI_adj_sorted_d01_T!=0]),bins=my_bins)
plt.close()
n_lpi_T, bins_lpi_log10_T, patches_lpi_T = plt.hist(np.log10(LPI_adj_sorted_d01_T[LPI_adj_sorted_d01_T!=0]),bins=my_bins)
plt.close()
n_LTG3_T, bins_LTG3_log10_T, patches_LTG3_T = plt.hist(np.log10(LTG3_adj_sorted_d01_T[LTG3_adj_sorted_d01_T!=0]),bins=my_bins)
plt.close()
n_PR92W_T, bins_PR92W_log10_T, patches_PR92W_T = plt.hist(np.log10(PR92W_adj_sorted_d01_T[PR92W_adj_sorted_d01_T!=0]),bins=my_bins)
plt.close()
# G4ICE - 9 km
n_CAPExP_R_G, bins_CAPExP_R_log10_G, patches_CAPExP_R_G = plt.hist(np.log10(CAPExP_R_adj_sorted_d01_G[CAPExP_R_adj_sorted_d01_G!=0]),bins=my_bins)
plt.close()
n_CAPExP_CSI_G, bins_CAPExP_CSI_log10_G, patches_CAPExP_CSI_G = plt.hist(np.log10(CAPExP_CSI_adj_sorted_d01_G[CAPExP_CSI_adj_sorted_d01_G!=0]),bins=my_bins)
plt.close()
n_lpi_G, bins_lpi_log10_G, patches_lpi_G = plt.hist(np.log10(LPI_adj_sorted_d01_G[LPI_adj_sorted_d01_G!=0]),bins=my_bins)
plt.close()
n_LTG3_G, bins_LTG3_log10_G, patches_LTG3_G = plt.hist(np.log10(LTG3_adj_sorted_d01_G[LTG3_adj_sorted_d01_G!=0]),bins=my_bins)
plt.close()
n_PR92W_G, bins_PR92W_log10_G, patches_PR92W_G = plt.hist(np.log10(PR92W_adj_sorted_d01_G[PR92W_adj_sorted_d01_G!=0]),bins=my_bins)
plt.close()

# set bins with 0 to the left to np.nan for a line plot, lines should just end to the left when there are no values
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
    if n_CAPExP_CSI_T[i] == 0:
        n_CAPExP_CSI_T[i] = np.nan
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
    if n_CAPExP_CSI_G[i] == 0:
        n_CAPExP_CSI_G[i] = np.nan
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

# Calculate center of bins to create line instead of bars
bins_log10_centered = bins_log10[0:-1]+0.5*(bins_log10[1:]-bins_log10[0:-1])
bins_centered = 10**bins_log10_centered

# THOM - 9 km
bins_CAPExP_R_log10_centered_T = bins_CAPExP_R_log10_T[0:-1]+0.5*(bins_CAPExP_R_log10_T[1:]-bins_CAPExP_R_log10_T[0:-1])
bins_CAPExP_R_centered_T = 10**bins_CAPExP_R_log10_centered_T

bins_lpi_log10_centered_T = bins_lpi_log10_T[0:-1]+0.5*(bins_lpi_log10_T[1:]-bins_lpi_log10_T[0:-1])
bins_lpi_centered_T = 10**bins_lpi_log10_centered_T

bins_LTG3_log10_centered_T = bins_LTG3_log10_T[0:-1]+0.5*(bins_LTG3_log10_T[1:]-bins_LTG3_log10_T[0:-1])
bins_LTG3_centered_T = 10**bins_LTG3_log10_centered_T

bins_CAPExP_CSI_log10_centered_T = bins_CAPExP_CSI_log10_T[0:-1]+0.5*(bins_CAPExP_CSI_log10_T[1:]-bins_CAPExP_CSI_log10_T[0:-1])
bins_CAPExP_CSI_centered_T = 10**bins_CAPExP_CSI_log10_centered_T

bins_PR92W_log10_centered_T = bins_PR92W_log10_T[0:-1]+0.5*(bins_PR92W_log10_T[1:]-bins_PR92W_log10_T[0:-1])
bins_PR92W_centered_T = 10**bins_PR92W_log10_centered_T

# G4ICE - 9 km
bins_CAPExP_R_log10_centered_G = bins_CAPExP_R_log10_G[0:-1]+0.5*(bins_CAPExP_R_log10_G[1:]-bins_CAPExP_R_log10_G[0:-1])
bins_CAPExP_R_centered_G = 10**bins_CAPExP_R_log10_centered_G

bins_lpi_log10_centered_G = bins_lpi_log10_G[0:-1]+0.5*(bins_lpi_log10_G[1:]-bins_lpi_log10_G[0:-1])
bins_lpi_centered_G = 10**bins_lpi_log10_centered_G

bins_LTG3_log10_centered_G = bins_LTG3_log10_G[0:-1]+0.5*(bins_LTG3_log10_G[1:]-bins_LTG3_log10_G[0:-1])
bins_LTG3_centered_G = 10**bins_LTG3_log10_centered_G

bins_CAPExP_CSI_log10_centered_G = bins_CAPExP_CSI_log10_G[0:-1]+0.5*(bins_CAPExP_CSI_log10_G[1:]-bins_CAPExP_CSI_log10_G[0:-1])
bins_CAPExP_CSI_centered_G = 10**bins_CAPExP_CSI_log10_centered_G

bins_PR92W_log10_centered_G = bins_PR92W_log10_G[0:-1]+0.5*(bins_PR92W_log10_G[1:]-bins_PR92W_log10_G[0:-1])
bins_PR92W_centered_G = 10**bins_PR92W_log10_centered_G

# DOMAIN 2
n, bins_log10, patches = plt.hist(np.log10(L[L!=0]),bins = my_bins)
plt.close()
# THOM - 3 km
n_CAPExP_R_d02_T, bins_CAPExP_R_log10_d02_T, patches_CAPExP_R_d02_T = plt.hist(np.log10(CAPExP_R_adj_sorted_d02_T[CAPExP_R_adj_sorted_d02_T!=0]),bins=my_bins)
plt.close()
n_CAPExP_CSI_d02_T, bins_CAPExP_CSI_log10_d02_T, patches_CAPExP_CSI_d02_T = plt.hist(np.log10(CAPExP_CSI_adj_sorted_d02_T[CAPExP_CSI_adj_sorted_d02_T!=0]),bins=my_bins)
plt.close()
n_lpi_d02_T, bins_lpi_log10_d02_T, patches_lpi_d02_T = plt.hist(np.log10(LPI_adj_sorted_d02_T[LPI_adj_sorted_d02_T!=0]),bins=my_bins)
plt.close()
n_LTG3_d02_T, bins_LTG3_log10_d02_T, patches_LTG3_d02_T = plt.hist(np.log10(LTG3_adj_sorted_d02_T[LTG3_adj_sorted_d02_T!=0]),bins=my_bins)
plt.close()
n_PR92W_d02_T, bins_PR92W_log10_d02_T, patches_PR92W_d02_T = plt.hist(np.log10(PR92W_adj_sorted_d02_T[PR92W_adj_sorted_d02_T!=0]),bins=my_bins)
plt.close()

# G4ICE - 3 km
n_CAPExP_R_d02_G, bins_CAPExP_R_log10_d02_G, patches_CAPExP_R_d02_G = plt.hist(np.log10(CAPExP_R_adj_sorted_d02_G[CAPExP_R_adj_sorted_d02_G!=0]),bins=my_bins)
plt.close()
n_CAPExP_CSI_d02_G, bins_CAPExP_CSI_log10_d02_G, patches_CAPExP_CSI_d02_G = plt.hist(np.log10(CAPExP_CSI_adj_sorted_d02_G[CAPExP_CSI_adj_sorted_d02_G!=0]),bins=my_bins)
plt.close()
n_lpi_d02_G, bins_lpi_log10_d02_G, patches_lpi_d02_G = plt.hist(np.log10(LPI_adj_sorted_d02_G[LPI_adj_sorted_d02_G!=0]),bins=my_bins)
plt.close()
n_LTG3_d02_G, bins_LTG3_log10_d02_G, patches_LTG3_d02_G = plt.hist(np.log10(LTG3_adj_sorted_d02_G[LTG3_adj_sorted_d02_G!=0]),bins=my_bins)
plt.close()
n_PR92W_d02_G, bins_PR92W_log10_d02_G, patches_PR92W_d02_G = plt.hist(np.log10(PR92W_adj_sorted_d02_G[PR92W_adj_sorted_d02_G!=0]),bins=my_bins)
plt.close()

# set bins with 0 to the left to np.nan for a line plot, lines should just end to the left when there are no values
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
    if n_CAPExP_CSI_d02_T[i] == 0:
        n_CAPExP_CSI_d02_T[i] = np.nan
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
    if n_CAPExP_CSI_d02_G[i] == 0:
        n_CAPExP_CSI_d02_G[i] = np.nan
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

# Calculate center of bins to create line instead of bars
# THOM - 3 km
bins_CAPExP_R_log10_centered_d02_T = bins_CAPExP_R_log10_d02_T[0:-1]+0.5*(bins_CAPExP_R_log10_d02_T[1:]-bins_CAPExP_R_log10_d02_T[0:-1])
bins_CAPExP_R_centered_d02_T = 10**bins_CAPExP_R_log10_centered_d02_T

bins_lpi_log10_centered_d02_T = bins_lpi_log10_d02_T[0:-1]+0.5*(bins_lpi_log10_d02_T[1:]-bins_lpi_log10_d02_T[0:-1])
bins_lpi_centered_d02_T = 10**bins_lpi_log10_centered_d02_T

bins_LTG3_log10_centered_d02_T = bins_LTG3_log10_d02_T[0:-1]+0.5*(bins_LTG3_log10_d02_T[1:]-bins_LTG3_log10_d02_T[0:-1])
bins_LTG3_centered_d02_T = 10**bins_LTG3_log10_centered_d02_T

bins_CAPExP_CSI_log10_centered_d02_T = bins_CAPExP_CSI_log10_d02_T[0:-1]+0.5*(bins_CAPExP_CSI_log10_d02_T[1:]-bins_CAPExP_CSI_log10_d02_T[0:-1])
bins_CAPExP_CSI_centered_d02_T = 10**bins_CAPExP_CSI_log10_centered_d02_T

bins_PR92W_log10_centered_d02_T = bins_PR92W_log10_d02_T[0:-1]+0.5*(bins_PR92W_log10_d02_T[1:]-bins_PR92W_log10_d02_T[0:-1])
bins_PR92W_centered_d02_T = 10**bins_PR92W_log10_centered_d02_T

# G4ICE - 3 km
bins_CAPExP_R_log10_centered_d02_G = bins_CAPExP_R_log10_d02_G[0:-1]+0.5*(bins_CAPExP_R_log10_d02_G[1:]-bins_CAPExP_R_log10_d02_G[0:-1])
bins_CAPExP_R_centered_d02_G = 10**bins_CAPExP_R_log10_centered_d02_G

bins_lpi_log10_centered_d02_G = bins_lpi_log10_d02_G[0:-1]+0.5*(bins_lpi_log10_d02_G[1:]-bins_lpi_log10_d02_G[0:-1])
bins_lpi_centered_d02_G = 10**bins_lpi_log10_centered_d02_G

bins_LTG3_log10_centered_d02_G = bins_LTG3_log10_d02_G[0:-1]+0.5*(bins_LTG3_log10_d02_G[1:]-bins_LTG3_log10_d02_G[0:-1])
bins_LTG3_centered_d02_G = 10**bins_LTG3_log10_centered_d02_G

bins_CAPExP_CSI_log10_centered_d02_G = bins_CAPExP_CSI_log10_d02_G[0:-1]+0.5*(bins_CAPExP_CSI_log10_d02_G[1:]-bins_CAPExP_CSI_log10_d02_G[0:-1])
bins_CAPExP_CSI_centered_d02_G = 10**bins_CAPExP_CSI_log10_centered_d02_G

bins_PR92W_log10_centered_d02_G = bins_PR92W_log10_d02_G[0:-1]+0.5*(bins_PR92W_log10_d02_G[1:]-bins_PR92W_log10_d02_G[0:-1])
bins_PR92W_centered_d02_G = 10**bins_PR92W_log10_centered_d02_G

# ---------------------------------------------------------------------------------------------
# PLOTS
# ---------------------------------------------------------------------------------------------
# Set all bins with frequency 0 to nan to get rid of vertical lines in the plots
n[n==0] = np.nan
n_lpi_T[n_lpi_T==0] = np.nan
n_LTG3_T[n_LTG3_T==0] = np.nan
n_PR92W_T[n_PR92W_T==0] = np.nan
n_CAPExP_R_T[n_CAPExP_R_T == 0] = np.nan
n_CAPExP_CSI_T[n_CAPExP_CSI_T==0] = np.nan

n_lpi_d02_T[n_lpi_d02_T==0] = np.nan
n_LTG3_d02_T[n_LTG3_d02_T==0] = np.nan
n_PR92W_d02_T[n_PR92W_d02_T==0] = np.nan
n_CAPExP_R_d02_T[n_CAPExP_R_d02_T == 0] = np.nan
n_CAPExP_CSI_d02_T[n_CAPExP_CSI_d02_T==0] = np.nan

n_lpi_G[n_lpi_G==0] = np.nan
n_LTG3_G[n_LTG3_G==0] = np.nan
n_PR92W_G[n_PR92W_G==0] = np.nan
n_CAPExP_R_G[n_CAPExP_R_G == 0] = np.nan
n_CAPExP_CSI_G[n_CAPExP_CSI_G==0] = np.nan

n_lpi_d02_G[n_lpi_d02_G==0] = np.nan
n_LTG3_d02_G[n_LTG3_d02_G==0] = np.nan
n_PR92W_d02_G[n_PR92W_d02_G==0] = np.nan
n_CAPExP_R_d02_G[n_CAPExP_R_d02_G == 0] = np.nan
n_CAPExP_CSI_d02_G[n_CAPExP_CSI_d02_G==0] = np.nan

# Select a colorblind palette to help the colorblind
palette = sns.color_palette('colorblind')

fig, axes = plt.subplots(2,2, figsize=(15.3, 8.27), dpi=150)

axes[0][0].loglog(bins_centered,n,'lightgrey', linewidth = 6)
axes[0][0].plot(bins_lpi_centered_G,n_lpi_G,marker='d', ms=7, color=palette[0])
axes[0][0].plot(bins_LTG3_centered_G,n_LTG3_G,marker='s', ms=7, color=palette[1])
axes[0][0].plot(bins_PR92W_centered_G,n_PR92W_G,marker='v', ms=7, color=palette[4])
axes[0][0].plot(bins_CAPExP_R_centered_G,n_CAPExP_R_G,marker='x', ms=7, color=palette[3])
axes[0][0].grid(which='major', axis='both', color='lightgray')
axes[0][0].set_title('G4ICE - 9 km', fontsize=14)
axes[0][0].set_yticks(ticks=[10, 1000, 100000])
axes[0][0].set_xlim(0.007, 30)
axes[0][0].tick_params(axis='both', which='major', labelsize=12)
axes[0][0].set_ylabel('Frequency', fontsize = 12)

axes[0][1].loglog(bins_centered,n,'lightgrey', linewidth = 6)
axes[0][1].plot(bins_lpi_centered_T,n_lpi_T,marker='d', ms=7, color=palette[0])
axes[0][1].plot(bins_LTG3_centered_T,n_LTG3_T,marker='s', ms=7, color=palette[1])
axes[0][1].plot(bins_PR92W_centered_T,n_PR92W_T,marker='v', ms=7, color=palette[4])
axes[0][1].plot(bins_CAPExP_R_centered_T,n_CAPExP_R_T,marker='x', ms=7, color=palette[3])
axes[0][1].grid(which='major', axis='both', color='lightgray')
axes[0][1].set_title('THOM - 9 km', fontsize=14)
axes[0][1].set_yticks(ticks=[10, 1000, 100000])
axes[0][1].set_xlim(0.007, 30)
axes[0][1].tick_params(axis='both', which='major', labelsize=12)

axes[1][0].loglog(bins_centered,n,'lightgrey', linewidth = 6)
axes[1][0].plot(bins_lpi_centered_d02_G,n_lpi_d02_G,marker='d', ms=7, color=palette[0])
axes[1][0].plot(bins_LTG3_centered_d02_G,n_LTG3_d02_G,marker='s', ms=7, color=palette[1])
axes[1][0].plot(bins_PR92W_centered_d02_G,n_PR92W_d02_G,marker='v', ms=7, color=palette[4])
axes[1][0].plot(bins_CAPExP_R_centered_d02_G,n_CAPExP_R_d02_G,marker='x', ms=7, color=palette[3])
axes[1][0].grid(which='major', axis='both', color='lightgray')
axes[1][0].set_title('G4ICE - 3 km', fontsize=14)
axes[1][0].set_yticks(ticks=[10, 1000, 100000])
axes[1][0].tick_params(axis='both', which='major', labelsize=12)
axes[1][0].set_ylabel('Frequency', fontsize = 12)
axes[1][0].set_xlim(0.007, 30)
axes[1][0].set_xlabel('Flash density (# hr$^{-1}$ km$^{-2}$)', fontsize = 12)

axes[1][1].loglog(bins_centered,n,'lightgrey', linewidth = 6, label='CLDN')
axes[1][1].plot(bins_lpi_centered_d02_T,n_lpi_d02_T,marker='d', ms=7, label='LPI', color=palette[0])
axes[1][1].plot(bins_LTG3_centered_d02_T,n_LTG3_d02_T,marker='s', ms=7, label='LT', color=palette[1])
axes[1][1].plot(bins_PR92W_centered_d02_T,n_PR92W_d02_T,marker='v', ms=7,  label='PR92W', color=palette[4])
axes[1][1].plot(bins_CAPExP_R_centered_d02_T,n_CAPExP_R_d02_T,marker='x', ms=7, label='CAPEXP', color=palette[3])
axes[1][1].grid(which='major', axis='both', color='lightgray')
axes[1][1].set_title('THOM - 3 km', fontsize=14)
axes[1][1].set_yticks(ticks=[10, 1000, 100000])
axes[1][1].tick_params(axis='both', which='major', labelsize=12)
axes[1][1].set_xlabel('Flash density (# hr$^{-1}$ km$^{-2}$)', fontsize = 12)
axes[1][1].legend(bbox_to_anchor=(0.5, -0.35), ncol=5, fontsize=12)
axes[1][1].set_xlim(0.007, 30)

plt.show()

# fig, axes = plt.subplots(2,2, figsize=(11.69, 8.27), dpi=600)
#
# axes[0][0].loglog(bins_centered,n,'lightgrey', linewidth = 6)
# axes[0][0].plot(bins_lpi_centered_G,n_lpi_G,marker='d', ms=7, color=palette[0])
# axes[0][0].plot(bins_LTG3_centered_G,n_LTG3_G,marker='s', ms=7, color=palette[1])
# axes[0][0].plot(bins_PR92W_centered_G,n_PR92W_G,marker='v', ms=7, color=palette[4])
# axes[0][0].plot(bins_CAPExP_R_centered_G,n_CAPExP_R_G,marker='x', ms=7, color=palette[3])
# # axes[0][0].plot(bins_CAPExP_CSI_centered_G,n_CAPExP_CSI_G,marker='P', ms=7, color=palette[4])
# axes[0][0].grid(which='major', axis='both', color='lightgray')
# axes[0][0].set_title('G4ICE - 9 km', fontsize=30)
# axes[0][0].set_yticks(ticks=[10, 1000, 100000])
# axes[0][0].set_xlim(0.007, 30)
# axes[0][0].tick_params(axis='both', which='major', labelsize=20)
# axes[0][0].set_ylabel('Frequency', fontsize = 28)
#
# axes[0][1].loglog(bins_centered,n,'lightgrey', linewidth = 6)
# axes[0][1].plot(bins_lpi_centered_T,n_lpi_T,marker='d', ms=7, color=palette[0])
# axes[0][1].plot(bins_LTG3_centered_T,n_LTG3_T,marker='s', ms=7, color=palette[1])
# axes[0][1].plot(bins_PR92W_centered_T,n_PR92W_T,marker='v', ms=7, color=palette[4])
# axes[0][1].plot(bins_CAPExP_R_centered_T,n_CAPExP_R_T,marker='x', ms=7, color=palette[3])
# # axes[0][1].plot(bins_CAPExP_CSI_centered_T,n_CAPExP_CSI_T,marker='P', ms=7, color=palette[4])
# axes[0][1].grid(which='major', axis='both', color='lightgray')
# axes[0][1].set_title('THOM - 9 km', fontsize=30)
# axes[0][1].set_yticks(ticks=[10, 1000, 100000])
# axes[0][1].set_xlim(0.007, 30)
# axes[0][1].tick_params(axis='both', which='major', labelsize=20)
#
# axes[1][0].loglog(bins_centered,n,'lightgrey', linewidth = 6)
# axes[1][0].plot(bins_lpi_centered_d02_G,n_lpi_d02_G,marker='d', ms=7, color=palette[0])
# axes[1][0].plot(bins_LTG3_centered_d02_G,n_LTG3_d02_G,marker='s', ms=7, color=palette[1])
# axes[1][0].plot(bins_PR92W_centered_d02_G,n_PR92W_d02_G,marker='v', ms=7, color=palette[4])
# axes[1][0].plot(bins_CAPExP_R_centered_d02_G,n_CAPExP_R_d02_G,marker='x', ms=7, color=palette[3])
# # axes[1][0].plot(bins_CAPExP_CSI_centered_d02_G,n_CAPExP_CSI_d02_G,marker='P', ms=7, color=palette[4])
# axes[1][0].grid(which='major', axis='both', color='lightgray')
# axes[1][0].set_title('G4ICE - 3 km', fontsize=30)
# axes[1][0].set_yticks(ticks=[10, 1000, 100000])
# axes[1][0].tick_params(axis='both', which='major', labelsize=20)
# axes[1][0].set_ylabel('Frequency', fontsize = 28)
# axes[1][0].set_xlim(0.007, 30)
# axes[1][0].set_xlabel('Flash density (# hr$^{-1}$ km$^{-2}$)', fontsize = 26)
#
# axes[1][1].loglog(bins_centered,n,'lightgrey', linewidth = 6, label='CLDN')
# axes[1][1].plot(bins_lpi_centered_d02_T,n_lpi_d02_T,marker='d', ms=7, label='LPI', color=palette[0])
# axes[1][1].plot(bins_LTG3_centered_d02_T,n_LTG3_d02_T,marker='s', ms=7, label='LT3', color=palette[1])
# axes[1][1].plot(bins_PR92W_centered_d02_T,n_PR92W_d02_T,marker='v', ms=7,  label='PR92W', color=palette[4])
# axes[1][1].plot(bins_CAPExP_R_centered_d02_T,n_CAPExP_R_d02_T,marker='x', ms=7, label='CAPEXP', color=palette[3])
# # axes[1][1].plot(bins_CAPExP_CSI_centered_d02_T,n_CAPExP_CSI_d02_T,marker='P', ms=7, label='CAPEXP_CSI', color=palette[4])
# axes[1][1].grid(which='major', axis='both', color='lightgray')
# axes[1][1].set_title('THOM - 3 km', fontsize=30)
# axes[1][1].set_yticks(ticks=[10, 1000, 100000])
# axes[1][1].tick_params(axis='both', which='major', labelsize=20)
# axes[1][1].set_xlabel('Flash density (# hr$^{-1}$ km$^{-2}$)', fontsize = 26)
# axes[1][1].legend(bbox_to_anchor=(0.7, -0.45), ncol=5,fontsize=20)
# axes[1][1].set_xlim(0.007, 30)
#
# plt.show()
