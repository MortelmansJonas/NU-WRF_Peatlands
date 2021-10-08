# ---------------------------------------------------------------------------------------------
# MODULES
# ---------------------------------------------------------------------------------------------
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from sklearn.linear_model import LinearRegression as lr
# ---------------------------------------------------------------------------------------------
# LOAD DATA
# ---------------------------------------------------------------------------------------------
ds_d01 = Dataset('/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/wrfout_nc_files/d01_all.nc','r')
ds_d02 = Dataset('/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/wrfout_nc_files/domain2_at_domain1_all.nc','r')
ds_obs = Dataset('/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/wrfout_nc_files/CLDN_at_domain1_all_v4.nc')

L = ds_obs['Flashdensity_CC'][:,:,:].data.flatten() + ds_obs['Flashdensity_CG'][:,:,:].data.flatten()
# ---------------------------------------------------------------------------------------------
# bins
# ---------------------------------------------------------------------------------------------
# specify bin edges
my_bins = np.linspace(-3,0,20)
# DOMAIN 1
n, bins_log10, patches = plt.hist(np.log10(L[L!=0]),bins = my_bins)
plt.close()
n_CAPExP_R, bins_CAPExP_R_log10, patches_CAPExP_R = plt.hist(np.log10(ds_d01['CAPExP_R'][:].flatten()),bins=my_bins)
plt.close()
n_CAPExP_CSI, bins_CAPExP_CSI_log10, patches_CAPExP_CSI = plt.hist(np.log10(ds_d01['CAPExP_CSI'][:].flatten()),bins=my_bins)
plt.close()
n_lpi, bins_lpi_log10, patches_lpi = plt.hist(np.log10(ds_d01['LPI'][:].flatten()),bins=my_bins)
plt.close()
n_LTG3, bins_LTG3_log10, patches_LTG3 = plt.hist(np.log10(ds_d01['LTG3'][:].flatten()),bins=my_bins)
plt.close()
n_PR92W, bins_PR92W_log10, patches_PR92W = plt.hist(np.log10(ds_d01['PR92W'][:].flatten()),bins=my_bins)
plt.close()

# reg = lr().fit(n_lpi, n)
# print(reg.score(n_lpi, n))
# print(reg.coef_)
# print(reg.intercept_)
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

plt.scatter(n_PR92W, n, label='PR92W')
plt.scatter(n_lpi, n, label='lpi')
plt.scatter(n_CAPExP_R, n, label='CAPExP_R')
plt.scatter(n_CAPExP_CSI, n, label='CAPExP_CSI')
plt.scatter(n_LTG3, n, label='LT3')
plt.legend()
plt.show()

n[np.isnan(n)] = 0
transformer = pf(degree=1)
x_ = transformer.fit_transform(n_lpi.reshape((-1,1)))
reg = lr().fit(x_, n)
print(reg.score(x_, n))
print(reg.coef_)
print(reg.intercept_)
