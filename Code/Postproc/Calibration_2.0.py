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
ds_d01 = Dataset('/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/wrfout_nc_files/d01_all.nc','r')
ds_d02 = Dataset('/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/wrfout_nc_files/domain2_at_domain1_all.nc','r')
ds_obs = Dataset('/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/wrfout_nc_files/CLDN_at_domain1_all_v4.nc')

lat = ds_d01['lat'][:]
lon = ds_d01['lon'][:]

L = ds_obs['Flashdensity_CC'][:,:,:].data.flatten() + ds_obs['Flashdensity_CG'][:,:,:].data.flatten()

# ---------------------------------------------------------------------------------------------
# SORTING
# ---------------------------------------------------------------------------------------------
# First sort non-zero observations
Obs = np.sort(L[L!=0])

# Sort model and take n highest values (with n = number of observations)
n = Obs.shape[0]

LPI_d01 = np.sort(ds_d01['LPI'][:].flatten())[::-1]
LPI_d01 = LPI_d01[0:n]
print(LPI_d01)

LTG3_d01 = np.sort(ds_d01['LTG3'][:].flatten())[::-1]
LTG3_d01 = LTG3_d01[0:n]
print(LTG3_d01)

PR92W_d01 = np.sort(ds_d01['PR92W'][:].flatten())[::-1]
PR92W_d01 = PR92W_d01[0:n]
print(PR92W_d01)

CAPExP_CSI_d01 = np.sort(ds_d01['CAPExP_CSI'][:].flatten())[::-1]
CAPExP_CSI_d01 = CAPExP_CSI_d01[0:n]
print(CAPExP_CSI_d01)

CAPExP_R_d01 = np.sort(ds_d01['CAPExP_R'][:].flatten())[::-1]
CAPExP_R_d01 = CAPExP_R_d01[0:n]
print(CAPExP_R_d01)

LPI_d02 = np.sort(ds_d02['LPI'][:].flatten())[::-1]
LPI_d02 = LPI_d02[0:n]
print(LPI_d02)

LTG3_d02 = np.sort(ds_d02['LTG3'][:].flatten())[::-1]
LTG3_d02 = LTG3_d02[0:n]
print(LTG3_d02)

PR92W_d02 = np.sort(ds_d02['PR92W'][:].flatten())[::-1]
PR92W_d02 = PR92W_d02[0:n]
print(PR92W_d02)

CAPExP_CSI_d02 = np.sort(ds_d02['CAPExP_CSI'][:].flatten())[::-1]
CAPExP_CSI_d02 = CAPExP_CSI_d02[0:n]
print(CAPExP_CSI_d02)

CAPExP_R_d02 = np.sort(ds_d02['CAPExP_R'][:].flatten())[::-1]
CAPExP_R_d02 = CAPExP_R_d02[0:n]
print(CAPExP_R_d02)

# ---------------------------------------------------------------------------------------------
# LINEAR REGRESSION
# ---------------------------------------------------------------------------------------------
# DOMAIN 1
reg_LPI_d01 = lr().fit(LPI_d01.reshape((-1,1)), Obs)
print('LPI_d01')
print(reg_LPI_d01.coef_)
print(reg_LPI_d01.intercept_)
print(reg_LPI_d01.score(LPI_d01.reshape((-1,1)), Obs))
linear_LPI_d01 =np.add(reg_LPI_d01.intercept_, np.multiply(reg_LPI_d01.coef_,ds_d01['LPI'][:]))

reg_LTG3_d01 = lr().fit(LTG3_d01.reshape((-1,1)), Obs)
print('LTG3_d01')
print(reg_LTG3_d01.coef_)
print(reg_LTG3_d01.intercept_)
print(reg_LTG3_d01.score(LTG3_d01.reshape((-1,1)), Obs))
# linear_LTG3_d01 =np.multiply(reg_LTG3_d01.coef_,ds_d01['LTG3'][:])
linear_LTG3_d01 =np.add(reg_LTG3_d01.intercept_, np.multiply(reg_LTG3_d01.coef_,ds_d01['LTG3'][:]))

reg_PR92W_d01 = lr().fit(PR92W_d01.reshape((-1,1)), Obs)
print('PR92W_d01')
print(reg_PR92W_d01.coef_)
print(reg_PR92W_d01.intercept_)
print(reg_PR92W_d01.score(PR92W_d01.reshape((-1,1)), Obs))
# linear_PR92W_d01 =np.multiply(reg_PR92W_d01.coef_,ds_d01['PR92W'][:])
linear_PR92W_d01 =np.add(reg_PR92W_d01.intercept_, np.multiply(reg_PR92W_d01.coef_,ds_d01['PR92W'][:]))

reg_CAPExP_CSI_d01 = lr().fit(CAPExP_CSI_d01.reshape((-1,1)), Obs)
print('CAPExP_CSI_d01')
print(reg_CAPExP_CSI_d01.coef_)
print(reg_CAPExP_CSI_d01.intercept_)
print(reg_CAPExP_CSI_d01.score(CAPExP_CSI_d01.reshape((-1,1)), Obs))
# linear_CAPExP_CSI_d01 =np.multiply(reg_CAPExP_CSI_d01.coef_,ds_d01['CAPExP_CSI'][:])
linear_CAPExP_CSI_d01 =np.add(reg_CAPExP_CSI_d01.intercept_, np.multiply(reg_CAPExP_CSI_d01.coef_,ds_d01['CAPExP_CSI'][:]))

reg_CAPExP_R_d01 = lr().fit(CAPExP_R_d01.reshape((-1,1)), Obs)
print('CAPExP_R_d01')
print(reg_CAPExP_R_d01.coef_)
print(reg_CAPExP_R_d01.intercept_)
print(reg_CAPExP_R_d01.score(CAPExP_R_d01.reshape((-1,1)), Obs))
# linear_CAPExP_R_d01 =np.multiply(reg_CAPExP_R_d01.coef_,ds_d01['CAPExP_R'][:])
linear_CAPExP_R_d01 =np.add(reg_CAPExP_R_d01.intercept_, np.multiply(reg_CAPExP_R_d01.coef_,ds_d01['CAPExP_R'][:]))

# DOMAIN 2
reg_LPI_d02 = lr().fit(LPI_d02.reshape((-1,1)), Obs)
print('LPI_d02')
print(reg_LPI_d02.coef_)
print(reg_LPI_d02.intercept_)
print(reg_LPI_d02.score(LPI_d02.reshape((-1,1)), Obs))
linear_LPI_d02 =np.add(reg_LPI_d02.intercept_, np.multiply(reg_LPI_d02.coef_,ds_d02['LPI'][:]))

reg_LTG3_d02 = lr().fit(LTG3_d02.reshape((-1,1)), Obs)
print('LTG3_d02')
print(reg_LTG3_d02.coef_)
print(reg_LTG3_d02.intercept_)
print(reg_LTG3_d02.score(LTG3_d02.reshape((-1,1)), Obs))
# linear_LTG3_d02 =np.multiply(reg_LTG3_d02.coef_,ds_d02['LTG3'][:])
linear_LTG3_d02 =np.add(reg_LTG3_d02.intercept_, np.multiply(reg_LTG3_d02.coef_,ds_d02['LTG3'][:]))

reg_PR92W_d02 = lr().fit(PR92W_d02.reshape((-1,1)), Obs)
print('PR92W_d02')
print(reg_PR92W_d02.coef_)
print(reg_PR92W_d02.intercept_)
print(reg_PR92W_d02.score(PR92W_d02.reshape((-1,1)), Obs))
# linear_PR92W_d02 =np.multiply(reg_PR92W_d02.coef_,ds_d02['PR92W'][:])
linear_PR92W_d02 =np.add(reg_PR92W_d02.intercept_, np.multiply(reg_PR92W_d02.coef_,ds_d02['PR92W'][:]))

reg_CAPExP_CSI_d02 = lr().fit(CAPExP_CSI_d02.reshape((-1,1)), Obs)
print('CAPExP_CSI_d02')
print(reg_CAPExP_CSI_d02.coef_)
print(reg_CAPExP_CSI_d02.intercept_)
print(reg_CAPExP_CSI_d02.score(CAPExP_CSI_d02.reshape((-1,1)), Obs))
# linear_CAPExP_CSI_d02 =np.multiply(reg_CAPExP_CSI_d02.coef_,ds_d02['CAPExP_CSI'][:])
linear_CAPExP_CSI_d02 =np.add(reg_CAPExP_CSI_d02.intercept_, np.multiply(reg_CAPExP_CSI_d02.coef_,ds_d02['CAPExP_CSI'][:]))

reg_CAPExP_R_d02 = lr().fit(CAPExP_R_d02.reshape((-1,1)), Obs)
print('CAPExP_R_d02')
print(reg_CAPExP_R_d02.coef_)
print(reg_CAPExP_R_d02.intercept_)
print(reg_CAPExP_R_d02.score(CAPExP_R_d02.reshape((-1,1)), Obs))
# linear_CAPExP_R_d02 =np.multiply(reg_CAPExP_R_d02.coef_,ds_d02['CAPExP_R'][:])
linear_CAPExP_R_d02 =np.add(reg_CAPExP_R_d02.intercept_, np.multiply(reg_CAPExP_R_d02.coef_,ds_d02['CAPExP_R'][:]))

# ---------------------------------------------------------------------------------------------
# PLOTS TO CHECK LINEAR REGRESSION
# ---------------------------------------------------------------------------------------------
plt.scatter(linear_PR92W_d01, L)
plt.ylabel('Observations')
plt.xlabel('Calibrated PR92W')
plt.title('d01')
plt.show()
