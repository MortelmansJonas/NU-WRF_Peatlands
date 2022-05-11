## import modules
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset

# load data
ds_G = Dataset('/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/NU-WRF/wrfout_nc_files/domain1_6h_72x72.nc', 'r')
ds_T = Dataset('/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/NU-WRF/wrfout_nc_files/domain1_6h_72x72_Thompson.nc', 'r')
# DOMAIN 1 GODDARD
# Create A, B, C, and D of contingency table
A_LPI_d01_G = np.count_nonzero(np.where((ds_G['LPI_d01'][:]> 0) & (ds_G['Obs'][:]>0), 1,0))
B_LPI_d01_G = np.count_nonzero(np.where((ds_G['LPI_d01'][:]> 0) & (ds_G['Obs'][:]==0), 1,0))
C_LPI_d01_G = np.count_nonzero(np.where((ds_G['LPI_d01'][:]== 0) & (ds_G['Obs'][:]>0), 1,0))
D_LPI_d01_G = np.count_nonzero(np.where((ds_G['LPI_d01'][:]== 0) & (ds_G['Obs'][:]==0), 1,0))

POD_LPI_d01_G = A_LPI_d01_G/(A_LPI_d01_G + C_LPI_d01_G)
FAR_LPI_d01_G = B_LPI_d01_G/(A_LPI_d01_G + B_LPI_d01_G)
Bias_LPI_d01_G = (A_LPI_d01_G + B_LPI_d01_G)/(A_LPI_d01_G + C_LPI_d01_G)
CSI_LPI_d01_G = A_LPI_d01_G/(A_LPI_d01_G + B_LPI_d01_G + C_LPI_d01_G)
SR_LPI_d01_G = 1-FAR_LPI_d01_G

A_LTG3_d01_G = np.count_nonzero(np.where((ds_G['LTG3_d01'][:]> 0) & (ds_G['Obs'][:]>0), 1,0))
B_LTG3_d01_G = np.count_nonzero(np.where((ds_G['LTG3_d01'][:]> 0) & (ds_G['Obs'][:]==0), 1,0))
C_LTG3_d01_G = np.count_nonzero(np.where((ds_G['LTG3_d01'][:]== 0) & (ds_G['Obs'][:]>0), 1,0))
D_LTG3_d01_G = np.count_nonzero(np.where((ds_G['LTG3_d01'][:]== 0) & (ds_G['Obs'][:]==0), 1,0))

POD_LTG3_d01_G = A_LTG3_d01_G/(A_LTG3_d01_G + C_LTG3_d01_G)
FAR_LTG3_d01_G = B_LTG3_d01_G/(A_LTG3_d01_G + B_LTG3_d01_G)
Bias_LTG3_d01_G = (A_LTG3_d01_G + B_LTG3_d01_G)/(A_LTG3_d01_G + C_LTG3_d01_G)
CSI_LTG3_d01_G = A_LTG3_d01_G/(A_LTG3_d01_G + B_LTG3_d01_G + C_LTG3_d01_G)
SR_LTG3_d01_G = 1-FAR_LTG3_d01_G

A_PR92W_d01_G = np.count_nonzero(np.where((ds_G['PR92W_d01'][:]> 0) & (ds_G['Obs'][:]>0), 1,0))
B_PR92W_d01_G = np.count_nonzero(np.where((ds_G['PR92W_d01'][:]> 0) & (ds_G['Obs'][:]==0), 1,0))
C_PR92W_d01_G = np.count_nonzero(np.where((ds_G['PR92W_d01'][:]== 0) & (ds_G['Obs'][:]>0), 1,0))
D_PR92W_d01_G = np.count_nonzero(np.where((ds_G['PR92W_d01'][:]== 0) & (ds_G['Obs'][:]==0), 1,0))

POD_PR92W_d01_G = A_PR92W_d01_G/(A_PR92W_d01_G + C_PR92W_d01_G)
FAR_PR92W_d01_G = B_PR92W_d01_G/(A_PR92W_d01_G + B_PR92W_d01_G)
Bias_PR92W_d01_G = (A_PR92W_d01_G + B_PR92W_d01_G)/(A_PR92W_d01_G + C_PR92W_d01_G)
CSI_PR92W_d01_G = A_PR92W_d01_G/(A_PR92W_d01_G + B_PR92W_d01_G + C_PR92W_d01_G)
SR_PR92W_d01_G = 1-FAR_PR92W_d01_G

A_CAPExP_R_d01_G = np.count_nonzero(np.where((ds_G['CAPExP_R_d01'][:]> 0) & (ds_G['Obs'][:]>0), 1,0))
B_CAPExP_R_d01_G = np.count_nonzero(np.where((ds_G['CAPExP_R_d01'][:]> 0) & (ds_G['Obs'][:]==0), 1,0))
C_CAPExP_R_d01_G = np.count_nonzero(np.where((ds_G['CAPExP_R_d01'][:]== 0) & (ds_G['Obs'][:]>0), 1,0))
D_CAPExP_R_d01_G = np.count_nonzero(np.where((ds_G['CAPExP_R_d01'][:]== 0) & (ds_G['Obs'][:]==0), 1,0))

POD_CAPExP_R_d01_G = A_CAPExP_R_d01_G/(A_CAPExP_R_d01_G + C_CAPExP_R_d01_G)
FAR_CAPExP_R_d01_G = B_CAPExP_R_d01_G/(A_CAPExP_R_d01_G + B_CAPExP_R_d01_G)
Bias_CAPExP_R_d01_G = (A_CAPExP_R_d01_G + B_CAPExP_R_d01_G)/(A_CAPExP_R_d01_G + C_CAPExP_R_d01_G)
CSI_CAPExP_R_d01_G = A_CAPExP_R_d01_G/(A_CAPExP_R_d01_G + B_CAPExP_R_d01_G + C_CAPExP_R_d01_G)
SR_CAPExP_R_d01_G = 1-FAR_CAPExP_R_d01_G

A_CAPExP_CSI_d01_G = np.count_nonzero(np.where((ds_G['CAPExP_CSI_d01'][:]> 0) & (ds_G['Obs'][:]>0), 1,0))
B_CAPExP_CSI_d01_G = np.count_nonzero(np.where((ds_G['CAPExP_CSI_d01'][:]> 0) & (ds_G['Obs'][:]==0), 1,0))
C_CAPExP_CSI_d01_G = np.count_nonzero(np.where((ds_G['CAPExP_CSI_d01'][:]== 0) & (ds_G['Obs'][:]>0), 1,0))
D_CAPExP_CSI_d01_G = np.count_nonzero(np.where((ds_G['CAPExP_CSI_d01'][:]== 0) & (ds_G['Obs'][:]==0), 1,0))

POD_CAPExP_CSI_d01_G = A_CAPExP_CSI_d01_G/(A_CAPExP_CSI_d01_G + C_CAPExP_CSI_d01_G)
FAR_CAPExP_CSI_d01_G = B_CAPExP_CSI_d01_G/(A_CAPExP_CSI_d01_G + B_CAPExP_CSI_d01_G)
Bias_CAPExP_CSI_d01_G = (A_CAPExP_CSI_d01_G + B_CAPExP_CSI_d01_G)/(A_CAPExP_CSI_d01_G + C_CAPExP_CSI_d01_G)
CSI_CAPExP_CSI_d01_G = A_CAPExP_CSI_d01_G/(A_CAPExP_CSI_d01_G + B_CAPExP_CSI_d01_G + C_CAPExP_CSI_d01_G)
SR_CAPExP_CSI_d01_G = 1-FAR_CAPExP_CSI_d01_G

# DOMAIN 2 GODDARD
# Create A, B, C, and D of contingency table
A_LPI_d02_G = np.count_nonzero(np.where((ds_G['LPI_d02'][:]> 0) & (ds_G['Obs'][:]>0), 1,0))
B_LPI_d02_G = np.count_nonzero(np.where((ds_G['LPI_d02'][:]> 0) & (ds_G['Obs'][:]==0), 1,0))
C_LPI_d02_G = np.count_nonzero(np.where((ds_G['LPI_d02'][:]== 0) & (ds_G['Obs'][:]>0), 1,0))
D_LPI_d02_G = np.count_nonzero(np.where((ds_G['LPI_d02'][:]== 0) & (ds_G['Obs'][:]==0), 1,0))

POD_LPI_d02_G = A_LPI_d02_G/(A_LPI_d02_G + C_LPI_d02_G)
FAR_LPI_d02_G = B_LPI_d02_G/(A_LPI_d02_G + B_LPI_d02_G)
Bias_LPI_d02_G = (A_LPI_d02_G + B_LPI_d02_G)/(A_LPI_d02_G + C_LPI_d02_G)
CSI_LPI_d02_G = A_LPI_d02_G/(A_LPI_d02_G + B_LPI_d02_G + C_LPI_d02_G)
SR_LPI_d02_G = 1-FAR_LPI_d02_G

A_LTG3_d02_G = np.count_nonzero(np.where((ds_G['LTG3_d02'][:]> 0) & (ds_G['Obs'][:]>0), 1,0))
B_LTG3_d02_G = np.count_nonzero(np.where((ds_G['LTG3_d02'][:]> 0) & (ds_G['Obs'][:]==0), 1,0))
C_LTG3_d02_G = np.count_nonzero(np.where((ds_G['LTG3_d02'][:]== 0) & (ds_G['Obs'][:]>0), 1,0))
D_LTG3_d02_G = np.count_nonzero(np.where((ds_G['LTG3_d02'][:]== 0) & (ds_G['Obs'][:]==0), 1,0))

POD_LTG3_d02_G = A_LTG3_d02_G/(A_LTG3_d02_G + C_LTG3_d02_G)
FAR_LTG3_d02_G = B_LTG3_d02_G/(A_LTG3_d02_G + B_LTG3_d02_G)
Bias_LTG3_d02_G = (A_LTG3_d02_G + B_LTG3_d02_G)/(A_LTG3_d02_G + C_LTG3_d02_G)
CSI_LTG3_d02_G = A_LTG3_d02_G/(A_LTG3_d02_G + B_LTG3_d02_G + C_LTG3_d02_G)
SR_LTG3_d02_G = 1-FAR_LTG3_d02_G

A_PR92W_d02_G = np.count_nonzero(np.where((ds_G['PR92W_d02'][:]> 0) & (ds_G['Obs'][:]>0), 1,0))
B_PR92W_d02_G = np.count_nonzero(np.where((ds_G['PR92W_d02'][:]> 0) & (ds_G['Obs'][:]==0), 1,0))
C_PR92W_d02_G = np.count_nonzero(np.where((ds_G['PR92W_d02'][:]== 0) & (ds_G['Obs'][:]>0), 1,0))
D_PR92W_d02_G = np.count_nonzero(np.where((ds_G['PR92W_d02'][:]== 0) & (ds_G['Obs'][:]==0), 1,0))

POD_PR92W_d02_G = A_PR92W_d02_G/(A_PR92W_d02_G + C_PR92W_d02_G)
FAR_PR92W_d02_G = B_PR92W_d02_G/(A_PR92W_d02_G + B_PR92W_d02_G)
Bias_PR92W_d02_G = (A_PR92W_d02_G + B_PR92W_d02_G)/(A_PR92W_d02_G + C_PR92W_d02_G)
CSI_PR92W_d02_G = A_PR92W_d02_G/(A_PR92W_d02_G + B_PR92W_d02_G + C_PR92W_d02_G)
SR_PR92W_d02_G = 1-FAR_PR92W_d02_G

A_CAPExP_R_d02_G = np.count_nonzero(np.where((ds_G['CAPExP_R_d02'][:]> 0) & (ds_G['Obs'][:]>0), 1,0))
B_CAPExP_R_d02_G = np.count_nonzero(np.where((ds_G['CAPExP_R_d02'][:]> 0) & (ds_G['Obs'][:]==0), 1,0))
C_CAPExP_R_d02_G = np.count_nonzero(np.where((ds_G['CAPExP_R_d02'][:]== 0) & (ds_G['Obs'][:]>0), 1,0))
D_CAPExP_R_d02_G = np.count_nonzero(np.where((ds_G['CAPExP_R_d02'][:]== 0) & (ds_G['Obs'][:]==0), 1,0))

POD_CAPExP_R_d02_G = A_CAPExP_R_d02_G/(A_CAPExP_R_d02_G + C_CAPExP_R_d02_G)
FAR_CAPExP_R_d02_G = B_CAPExP_R_d02_G/(A_CAPExP_R_d02_G + B_CAPExP_R_d02_G)
Bias_CAPExP_R_d02_G = (A_CAPExP_R_d02_G + B_CAPExP_R_d02_G)/(A_CAPExP_R_d02_G + C_CAPExP_R_d02_G)
CSI_CAPExP_R_d02_G = A_CAPExP_R_d02_G/(A_CAPExP_R_d02_G + B_CAPExP_R_d02_G + C_CAPExP_R_d02_G)
SR_CAPExP_R_d02_G = 1-FAR_CAPExP_R_d02_G

A_CAPExP_CSI_d02_G = np.count_nonzero(np.where((ds_G['CAPExP_CSI_d02'][:]> 0) & (ds_G['Obs'][:]>0), 1,0))
B_CAPExP_CSI_d02_G = np.count_nonzero(np.where((ds_G['CAPExP_CSI_d02'][:]> 0) & (ds_G['Obs'][:]==0), 1,0))
C_CAPExP_CSI_d02_G = np.count_nonzero(np.where((ds_G['CAPExP_CSI_d02'][:]== 0) & (ds_G['Obs'][:]>0), 1,0))
D_CAPExP_CSI_d02_G = np.count_nonzero(np.where((ds_G['CAPExP_CSI_d02'][:]== 0) & (ds_G['Obs'][:]==0), 1,0))

POD_CAPExP_CSI_d02_G = A_CAPExP_CSI_d02_G/(A_CAPExP_CSI_d02_G + C_CAPExP_CSI_d02_G)
FAR_CAPExP_CSI_d02_G = B_CAPExP_CSI_d02_G/(A_CAPExP_CSI_d02_G + B_CAPExP_CSI_d02_G)
Bias_CAPExP_CSI_d02_G = (A_CAPExP_CSI_d02_G + B_CAPExP_CSI_d02_G)/(A_CAPExP_CSI_d02_G + C_CAPExP_CSI_d02_G)
CSI_CAPExP_CSI_d02_G = A_CAPExP_CSI_d02_G/(A_CAPExP_CSI_d02_G + B_CAPExP_CSI_d02_G + C_CAPExP_CSI_d02_G)
SR_CAPExP_CSI_d02_G = 1-FAR_CAPExP_CSI_d02_G

# DOMAIN 1 THOMPSON
# Create A, B, C, and D of contingency table
A_LPI_d01_T = np.count_nonzero(np.where((ds_T['LPI_d01'][:]> 0) & (ds_T['Obs'][:]>0), 1,0))
B_LPI_d01_T = np.count_nonzero(np.where((ds_T['LPI_d01'][:]> 0) & (ds_T['Obs'][:]==0), 1,0))
C_LPI_d01_T = np.count_nonzero(np.where((ds_T['LPI_d01'][:]== 0) & (ds_T['Obs'][:]>0), 1,0))
D_LPI_d01_T = np.count_nonzero(np.where((ds_T['LPI_d01'][:]== 0) & (ds_T['Obs'][:]==0), 1,0))

POD_LPI_d01_T = A_LPI_d01_T/(A_LPI_d01_T + C_LPI_d01_T)
FAR_LPI_d01_T = B_LPI_d01_T/(A_LPI_d01_T + B_LPI_d01_T)
Bias_LPI_d01_T = (A_LPI_d01_T + B_LPI_d01_T)/(A_LPI_d01_T + C_LPI_d01_T)
CSI_LPI_d01_T = A_LPI_d01_T/(A_LPI_d01_T + B_LPI_d01_T + C_LPI_d01_T)
SR_LPI_d01_T = 1-FAR_LPI_d01_T

A_LTG3_d01_T = np.count_nonzero(np.where((ds_T['LTG3_d01'][:]> 0) & (ds_T['Obs'][:]>0), 1,0))
B_LTG3_d01_T = np.count_nonzero(np.where((ds_T['LTG3_d01'][:]> 0) & (ds_T['Obs'][:]==0), 1,0))
C_LTG3_d01_T = np.count_nonzero(np.where((ds_T['LTG3_d01'][:]== 0) & (ds_T['Obs'][:]>0), 1,0))
D_LTG3_d01_T = np.count_nonzero(np.where((ds_T['LTG3_d01'][:]== 0) & (ds_T['Obs'][:]==0), 1,0))

POD_LTG3_d01_T = A_LTG3_d01_T/(A_LTG3_d01_T + C_LTG3_d01_T)
FAR_LTG3_d01_T = B_LTG3_d01_T/(A_LTG3_d01_T + B_LTG3_d01_T)
Bias_LTG3_d01_T = (A_LTG3_d01_T + B_LTG3_d01_T)/(A_LTG3_d01_T + C_LTG3_d01_T)
CSI_LTG3_d01_T = A_LTG3_d01_T/(A_LTG3_d01_T + B_LTG3_d01_T + C_LTG3_d01_T)
SR_LTG3_d01_T = 1-FAR_LTG3_d01_T

A_PR92W_d01_T = np.count_nonzero(np.where((ds_T['PR92W_d01'][:]> 0) & (ds_T['Obs'][:]>0), 1,0))
B_PR92W_d01_T = np.count_nonzero(np.where((ds_T['PR92W_d01'][:]> 0) & (ds_T['Obs'][:]==0), 1,0))
C_PR92W_d01_T = np.count_nonzero(np.where((ds_T['PR92W_d01'][:]== 0) & (ds_T['Obs'][:]>0), 1,0))
D_PR92W_d01_T = np.count_nonzero(np.where((ds_T['PR92W_d01'][:]== 0) & (ds_T['Obs'][:]==0), 1,0))

POD_PR92W_d01_T = A_PR92W_d01_T/(A_PR92W_d01_T + C_PR92W_d01_T)
FAR_PR92W_d01_T = B_PR92W_d01_T/(A_PR92W_d01_T + B_PR92W_d01_T)
Bias_PR92W_d01_T = (A_PR92W_d01_T + B_PR92W_d01_T)/(A_PR92W_d01_T + C_PR92W_d01_T)
CSI_PR92W_d01_T = A_PR92W_d01_T/(A_PR92W_d01_T + B_PR92W_d01_T + C_PR92W_d01_T)
SR_PR92W_d01_T = 1-FAR_PR92W_d01_T

A_CAPExP_R_d01_T = np.count_nonzero(np.where((ds_T['CAPExP_R_d01'][:]> 0) & (ds_T['Obs'][:]>0), 1,0))
B_CAPExP_R_d01_T = np.count_nonzero(np.where((ds_T['CAPExP_R_d01'][:]> 0) & (ds_T['Obs'][:]==0), 1,0))
C_CAPExP_R_d01_T = np.count_nonzero(np.where((ds_T['CAPExP_R_d01'][:]== 0) & (ds_T['Obs'][:]>0), 1,0))
D_CAPExP_R_d01_T = np.count_nonzero(np.where((ds_T['CAPExP_R_d01'][:]== 0) & (ds_T['Obs'][:]==0), 1,0))

POD_CAPExP_R_d01_T = A_CAPExP_R_d01_T/(A_CAPExP_R_d01_T + C_CAPExP_R_d01_T)
FAR_CAPExP_R_d01_T = B_CAPExP_R_d01_T/(A_CAPExP_R_d01_T + B_CAPExP_R_d01_T)
Bias_CAPExP_R_d01_T = (A_CAPExP_R_d01_T + B_CAPExP_R_d01_T)/(A_CAPExP_R_d01_T + C_CAPExP_R_d01_T)
CSI_CAPExP_R_d01_T = A_CAPExP_R_d01_T/(A_CAPExP_R_d01_T + B_CAPExP_R_d01_T + C_CAPExP_R_d01_T)
SR_CAPExP_R_d01_T = 1-FAR_CAPExP_R_d01_T

A_CAPExP_CSI_d01_T = np.count_nonzero(np.where((ds_T['CAPExP_CSI_d01'][:]> 0) & (ds_T['Obs'][:]>0), 1,0))
B_CAPExP_CSI_d01_T = np.count_nonzero(np.where((ds_T['CAPExP_CSI_d01'][:]> 0) & (ds_T['Obs'][:]==0), 1,0))
C_CAPExP_CSI_d01_T = np.count_nonzero(np.where((ds_T['CAPExP_CSI_d01'][:]== 0) & (ds_T['Obs'][:]>0), 1,0))
D_CAPExP_CSI_d01_T = np.count_nonzero(np.where((ds_T['CAPExP_CSI_d01'][:]== 0) & (ds_T['Obs'][:]==0), 1,0))

POD_CAPExP_CSI_d01_T = A_CAPExP_CSI_d01_T/(A_CAPExP_CSI_d01_T + C_CAPExP_CSI_d01_T)
FAR_CAPExP_CSI_d01_T = B_CAPExP_CSI_d01_T/(A_CAPExP_CSI_d01_T + B_CAPExP_CSI_d01_T)
Bias_CAPExP_CSI_d01_T = (A_CAPExP_CSI_d01_T + B_CAPExP_CSI_d01_T)/(A_CAPExP_CSI_d01_T + C_CAPExP_CSI_d01_T)
CSI_CAPExP_CSI_d01_T = A_CAPExP_CSI_d01_T/(A_CAPExP_CSI_d01_T + B_CAPExP_CSI_d01_T + C_CAPExP_CSI_d01_T)
SR_CAPExP_CSI_d01_T = 1-FAR_CAPExP_CSI_d01_T

# DOMAIN 2 THOMPSON
# Create A, B, C, and D of contingency table
A_LPI_d02_T = np.count_nonzero(np.where((ds_T['LPI_d02'][:]> 0) & (ds_T['Obs'][:]>0), 1,0))
B_LPI_d02_T = np.count_nonzero(np.where((ds_T['LPI_d02'][:]> 0) & (ds_T['Obs'][:]==0), 1,0))
C_LPI_d02_T = np.count_nonzero(np.where((ds_T['LPI_d02'][:]== 0) & (ds_T['Obs'][:]>0), 1,0))
D_LPI_d02_T = np.count_nonzero(np.where((ds_T['LPI_d02'][:]== 0) & (ds_T['Obs'][:]==0), 1,0))

POD_LPI_d02_T = A_LPI_d02_T/(A_LPI_d02_T + C_LPI_d02_T)
FAR_LPI_d02_T = B_LPI_d02_T/(A_LPI_d02_T + B_LPI_d02_T)
Bias_LPI_d02_T = (A_LPI_d02_T + B_LPI_d02_T)/(A_LPI_d02_T + C_LPI_d02_T)
CSI_LPI_d02_T = A_LPI_d02_T/(A_LPI_d02_T + B_LPI_d02_T + C_LPI_d02_T)
SR_LPI_d02_T = 1-FAR_LPI_d02_T

A_LTG3_d02_T = np.count_nonzero(np.where((ds_T['LTG3_d02'][:]> 0) & (ds_T['Obs'][:]>0), 1,0))
B_LTG3_d02_T = np.count_nonzero(np.where((ds_T['LTG3_d02'][:]> 0) & (ds_T['Obs'][:]==0), 1,0))
C_LTG3_d02_T = np.count_nonzero(np.where((ds_T['LTG3_d02'][:]== 0) & (ds_T['Obs'][:]>0), 1,0))
D_LTG3_d02_T = np.count_nonzero(np.where((ds_T['LTG3_d02'][:]== 0) & (ds_T['Obs'][:]==0), 1,0))

POD_LTG3_d02_T = A_LTG3_d02_T/(A_LTG3_d02_T + C_LTG3_d02_T)
FAR_LTG3_d02_T = B_LTG3_d02_T/(A_LTG3_d02_T + B_LTG3_d02_T)
Bias_LTG3_d02_T = (A_LTG3_d02_T + B_LTG3_d02_T)/(A_LTG3_d02_T + C_LTG3_d02_T)
CSI_LTG3_d02_T = A_LTG3_d02_T/(A_LTG3_d02_T + B_LTG3_d02_T + C_LTG3_d02_T)
SR_LTG3_d02_T = 1-FAR_LTG3_d02_T

A_PR92W_d02_T = np.count_nonzero(np.where((ds_T['PR92W_d02'][:]> 0) & (ds_T['Obs'][:]>0), 1,0))
B_PR92W_d02_T = np.count_nonzero(np.where((ds_T['PR92W_d02'][:]> 0) & (ds_T['Obs'][:]==0), 1,0))
C_PR92W_d02_T = np.count_nonzero(np.where((ds_T['PR92W_d02'][:]== 0) & (ds_T['Obs'][:]>0), 1,0))
D_PR92W_d02_T = np.count_nonzero(np.where((ds_T['PR92W_d02'][:]== 0) & (ds_T['Obs'][:]==0), 1,0))

POD_PR92W_d02_T = A_PR92W_d02_T/(A_PR92W_d02_T + C_PR92W_d02_T)
FAR_PR92W_d02_T = B_PR92W_d02_T/(A_PR92W_d02_T + B_PR92W_d02_T)
Bias_PR92W_d02_T = (A_PR92W_d02_T + B_PR92W_d02_T)/(A_PR92W_d02_T + C_PR92W_d02_T)
CSI_PR92W_d02_T = A_PR92W_d02_T/(A_PR92W_d02_T + B_PR92W_d02_T + C_PR92W_d02_T)
SR_PR92W_d02_T = 1-FAR_PR92W_d02_T

A_CAPExP_R_d02_T = np.count_nonzero(np.where((ds_T['CAPExP_R_d02'][:]> 0) & (ds_T['Obs'][:]>0), 1,0))
B_CAPExP_R_d02_T = np.count_nonzero(np.where((ds_T['CAPExP_R_d02'][:]> 0) & (ds_T['Obs'][:]==0), 1,0))
C_CAPExP_R_d02_T = np.count_nonzero(np.where((ds_T['CAPExP_R_d02'][:]== 0) & (ds_T['Obs'][:]>0), 1,0))
D_CAPExP_R_d02_T = np.count_nonzero(np.where((ds_T['CAPExP_R_d02'][:]== 0) & (ds_T['Obs'][:]==0), 1,0))

POD_CAPExP_R_d02_T = A_CAPExP_R_d02_T/(A_CAPExP_R_d02_T + C_CAPExP_R_d02_T)
FAR_CAPExP_R_d02_T = B_CAPExP_R_d02_T/(A_CAPExP_R_d02_T + B_CAPExP_R_d02_T)
Bias_CAPExP_R_d02_T = (A_CAPExP_R_d02_T + B_CAPExP_R_d02_T)/(A_CAPExP_R_d02_T + C_CAPExP_R_d02_T)
CSI_CAPExP_R_d02_T = A_CAPExP_R_d02_T/(A_CAPExP_R_d02_T + B_CAPExP_R_d02_T + C_CAPExP_R_d02_T)
SR_CAPExP_R_d02_T = 1-FAR_CAPExP_R_d02_T

A_CAPExP_CSI_d02_T = np.count_nonzero(np.where((ds_T['CAPExP_CSI_d02'][:]> 0) & (ds_T['Obs'][:]>0), 1,0))
B_CAPExP_CSI_d02_T = np.count_nonzero(np.where((ds_T['CAPExP_CSI_d02'][:]> 0) & (ds_T['Obs'][:]==0), 1,0))
C_CAPExP_CSI_d02_T = np.count_nonzero(np.where((ds_T['CAPExP_CSI_d02'][:]== 0) & (ds_T['Obs'][:]>0), 1,0))
D_CAPExP_CSI_d02_T = np.count_nonzero(np.where((ds_T['CAPExP_CSI_d02'][:]== 0) & (ds_T['Obs'][:]==0), 1,0))

POD_CAPExP_CSI_d02_T = A_CAPExP_CSI_d02_T/(A_CAPExP_CSI_d02_T + C_CAPExP_CSI_d02_T)
FAR_CAPExP_CSI_d02_T = B_CAPExP_CSI_d02_T/(A_CAPExP_CSI_d02_T + B_CAPExP_CSI_d02_T)
Bias_CAPExP_CSI_d02_T = (A_CAPExP_CSI_d02_T + B_CAPExP_CSI_d02_T)/(A_CAPExP_CSI_d02_T + C_CAPExP_CSI_d02_T)
CSI_CAPExP_CSI_d02_T = A_CAPExP_CSI_d02_T/(A_CAPExP_CSI_d02_T + B_CAPExP_CSI_d02_T + C_CAPExP_CSI_d02_T)
SR_CAPExP_CSI_d02_T = 1-FAR_CAPExP_CSI_d02_T

## PLOT

fig = plt.figure(figsize=(15.3, 8.27), dpi=150)

grid_ticks = np.arange(0, 1.01, 0.01)
sr_g, pod_g = np.meshgrid(grid_ticks, grid_ticks)
bias = pod_g / sr_g
csi = 1.0 / (1.0 / sr_g + 1.0 / pod_g - 1.0)
csi_contour = plt.contourf(sr_g, pod_g, csi, np.arange(0.1, 1.1, 0.1), cmap="Blues")
b_contour = plt.contour(sr_g, pod_g, bias, [0.3, 0.5, 0.7, 1, 1.5, 2, 4], colors="k", linestyles="dashed")
plt.clabel(b_contour, fmt="%1.1f", fontsize = 12)
plt.scatter(SR_LPI_d01_G, POD_LPI_d01_G, color='tab:blue', marker='o', s=30)
plt.scatter(SR_LTG3_d01_G, POD_LTG3_d01_G, color='tab:blue', marker='v', s=30)
plt.scatter(SR_PR92W_d01_G, POD_PR92W_d01_G, color='tab:blue', marker='d', s=30)
# plt.scatter(SR_CAPExP_CSI_d01_G, POD_CAPExP_CSI_d01_G, color='tab:blue', marker='s')
plt.scatter(SR_CAPExP_R_d01_G, POD_CAPExP_R_d01_G, color='tab:blue', marker='x', s=30)
plt.scatter(SR_LPI_d02_G, POD_LPI_d02_G, color='tab:brown', marker='o', s=30)
plt.scatter(SR_LTG3_d02_G, POD_LTG3_d02_G, color='tab:brown', marker='v', s=30)
plt.scatter(SR_PR92W_d02_G, POD_PR92W_d02_G, color='tab:brown', marker='d', s=30)
# plt.scatter(SR_CAPExP_CSI_d02_G, POD_CAPExP_CSI_d02_G, color='tab:brown', marker='s')
plt.scatter(SR_CAPExP_R_d02_G, POD_CAPExP_R_d02_G, color='tab:brown', marker='x', s=30)
plt.scatter(SR_LPI_d01_T, POD_LPI_d01_T, color='tab:pink', marker='o', s=30)
plt.scatter(SR_LTG3_d01_T, POD_LTG3_d01_T, color='tab:pink', marker='v', s=30)
plt.scatter(SR_PR92W_d01_T, POD_PR92W_d01_T, color='tab:pink', marker='d', s=30)
# plt.scatter(SR_CAPExP_CSI_d01_T, POD_CAPExP_CSI_d01_T, color='tab:pink', marker='s')
plt.scatter(SR_CAPExP_R_d01_T, POD_CAPExP_R_d01_T, color='tab:pink', marker='x', s=30)
plt.scatter(SR_LPI_d02_T, POD_LPI_d02_T, color='tab:orange', marker='o', s=30)
plt.scatter(SR_LTG3_d02_T, POD_LTG3_d02_T, color='tab:orange', marker='v', s=30)
plt.scatter(SR_PR92W_d02_T, POD_PR92W_d02_T, color='tab:orange', marker='d', s=30)
# plt.scatter(SR_CAPExP_CSI_d02_T, POD_CAPExP_CSI_d02_T, color='tab:orange', marker='s')
plt.scatter(SR_CAPExP_R_d02_T, POD_CAPExP_R_d02_T, color='tab:orange', marker='x', s=30)

# Create legend
plt.scatter([],[], color='tab:blue', marker='s', label='G4ICE - 9 km')
plt.scatter([],[], color='tab:brown', marker='s', label='G4ICE - 3 km')
plt.scatter([],[], color='tab:pink', marker='s', label='THOM - 9 km')
plt.scatter([],[], color='tab:orange', marker='s', label='THOM - 3 km')
plt.scatter([],[], color='k', marker='o', label='LPI')
plt.scatter([],[], color='k', marker='v', label='LT')
plt.scatter([],[], color='k', marker='d', label='PR92W')
plt.scatter([],[], color='k', marker='x', label='CAPExP')

cbar = plt.colorbar(csi_contour)
cbar.set_label("Critical Success Index", fontsize=12)
cbar.ax.tick_params(labelsize=12)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.ylabel("Probability of Detection", fontsize=12)
plt.xlabel("Success Ratio", fontsize=12)
plt.legend(fontsize=12, loc='lower right')
plt.show()
