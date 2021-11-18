# ---------------------------------------------------------------------------------------------
# MODULES
# ---------------------------------------------------------------------------------------------
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset

# ---------------------------------------------------------------------------------------------
# LOAD DATA
# ---------------------------------------------------------------------------------------------
ds_d01 = Dataset('/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/wrfout_nc_files/d01_diurnal.nc', 'r')
ds_d02 = Dataset('/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/wrfout_nc_files/d02_diurnal.nc', 'r')

Hours = np.linspace(0,24, num=24, endpoint=False)

# ---------------------------------------------------------------------------------------------
# NORMALIZE DATA
# ---------------------------------------------------------------------------------------------
LPI_d01 = (ds_d01['LPI'][:])
LPI_d02 = (ds_d02['LPI'][:])
LTG3_d01 = (ds_d01['LTG3'][:])
LTG3_d02 = (ds_d02['LTG3'][:])
PR92W_d01 = (ds_d01['PR92W'][:])
PR92W_d02 = (ds_d02['PR92W'][:])
CAPExP_R_d01 = (ds_d01['CAPExP_R'][:])
CAPExP_R_d02 = (ds_d02['CAPExP_R'][:])
CAPExP_CSI_d01 = (ds_d01['CAPExP_CSI'][:])
CAPExP_CSI_d02 = (ds_d02['CAPExP_CSI'][:])
OBS_d01 = (ds_d01['Obs'][:])
OBS_d02 = (ds_d02['Obs'][:])

# ---------------------------------------------------------------------------------------------
# PLOT D01
# ---------------------------------------------------------------------------------------------
x1 = 0.0005

fig = plt.figure()
ax1 = plt.subplot2grid((3,2),(0,0))
ax2 = plt.subplot2grid((3,2),(1,0))
ax3 = plt.subplot2grid((3,2),(0,1))
ax4 = plt.subplot2grid((3,2),(1,1))
ax5 = plt.subplot2grid((3,2),(2,0))
ax6 = plt.subplot2grid((3,2),(2,1))

ax1.plot(Hours, OBS_d01, 'grey', label='CLDN')
ax1.set_title('CLDN')
ax1.annotate('(a) \n', xy=(x1, ax1.get_ylim()[1]),annotation_clip=False)
ax1.set_ylabel('Hourly flash density \n(# hr$^{-1}$ km$^{-2}$)')
ax1.grid(which='major', axis='both', color='lightgray')

# ax2.plot(Hours, OBS_d01, 'grey', label='CLDN')
ax2.plot(Hours, LTG3_d01, 'k', label='LT3')
ax2.set_title('LT3')
ax2.annotate('(c) \n', xy=(x1, ax2.get_ylim()[1]),annotation_clip=False)
ax2.set_ylabel('Hourly flash density \n(# hr$^{-1}$ km$^{-2}$)')
ax2.grid(which='major', axis='both', color='lightgray')

ax3.plot(Hours, OBS_d01, 'grey', label='CLDN')
ax3.plot(Hours, LPI_d01, 'k', label='LPI')
ax3.set_title('LPI')
ax3.annotate('(b) \n', xy=(x1, ax3.get_ylim()[1]),annotation_clip=False)
ax3.grid(which='major', axis='both', color='lightgray')

ax4.plot(Hours, OBS_d01, 'grey', label='CLDN')
ax4.plot(Hours, PR92W_d01, 'k', label='PR92W')
ax4.set_title('PR92W')
ax4.annotate('(d) \n', xy=(x1, ax4.get_ylim()[1]),annotation_clip=False)
ax4.grid(which='major', axis='both', color='lightgray')

ax5.plot(Hours, OBS_d01, 'grey', label='CLDN')
ax5.plot(Hours, CAPExP_R_d01, 'k', label='CAPExP_R')
ax5.set_title('CAPExP_R')
ax5.annotate('(e) \n', xy=(x1, ax5.get_ylim()[1]),annotation_clip=False)
ax5.set_ylabel('Hourly flash density \n(# hr$^{-1}$ km$^{-2}$)')
ax5.grid(which='major', axis='both', color='lightgray')

ax6.plot(Hours, OBS_d01, 'grey', label='CLDN')
ax6.plot(Hours, CAPExP_CSI_d01, 'k', label='CAPExP_CSI')
ax6.set_title('CAPExP_CSI')
ax6.annotate('(f) \n', xy=(x1, ax6.get_ylim()[1]),annotation_clip=False)
ax6.grid(which='major', axis='both', color='lightgray')

fig.suptitle('Convection-parameterized (9 km)')
fig.tight_layout()
plt.show()

# ---------------------------------------------------------------------------------------------
# PLOT D01
# ---------------------------------------------------------------------------------------------
fig = plt.figure()
ax1 = plt.subplot2grid((3,2),(0,0))
ax2 = plt.subplot2grid((3,2),(1,0))
ax3 = plt.subplot2grid((3,2),(0,1))
ax4 = plt.subplot2grid((3,2),(1,1))
ax5 = plt.subplot2grid((3,2),(2,0))
ax6 = plt.subplot2grid((3,2),(2,1))

ax1.plot(Hours, OBS_d02, 'grey', label='CLDN')
ax1.set_title('CLDN')
ax1.annotate('(a) \n', xy=(x1, ax1.get_ylim()[1]),annotation_clip=False)
ax1.set_ylabel('Hourly flash density \n(# hr$^{-1}$ km$^{-2}$)')
ax1.grid(which='major', axis='both', color='lightgray')

# ax2.plot(Hours, OBS_d02, 'grey', label='CLDN')
ax2.plot(Hours, LTG3_d02, 'k', label='LT3')
ax2.set_title('LT3')
ax2.annotate('(c) \n', xy=(x1, ax2.get_ylim()[1]),annotation_clip=False)
ax2.set_ylabel('Hourly flash density \n(# hr$^{-1}$ km$^{-2}$)')
ax2.grid(which='major', axis='both', color='lightgray')

ax3.plot(Hours, OBS_d02, 'grey', label='CLDN')
ax3.plot(Hours, LPI_d02, 'k', label='LPI')
ax3.set_title('LPI')
ax3.annotate('(b) \n', xy=(x1, ax3.get_ylim()[1]),annotation_clip=False)
ax3.grid(which='major', axis='both', color='lightgray')

ax4.plot(Hours, OBS_d02, 'grey', label='CLDN')
ax4.plot(Hours, PR92W_d02, 'k', label='PR92W')
ax4.set_title('PR92W')
ax4.annotate('(d) \n', xy=(x1, ax4.get_ylim()[1]),annotation_clip=False)
ax4.grid(which='major', axis='both', color='lightgray')

ax5.plot(Hours, OBS_d02, 'grey', label='CLDN')
ax5.plot(Hours, CAPExP_R_d02, 'k', label='CAPExP_R')
ax5.set_title('CAPExP_R')
ax5.annotate('(e) \n', xy=(x1, ax5.get_ylim()[1]),annotation_clip=False)
ax5.set_ylabel('Hourly flash density \n(# hr$^{-1}$ km$^{-2}$)')
ax5.grid(which='major', axis='both', color='lightgray')

ax6.plot(Hours, OBS_d02, 'grey', label='CLDN')
ax6.plot(Hours, CAPExP_CSI_d02, 'k', label='CAPExP_CSI')
ax6.set_title('CAPExP_CSI')
ax6.annotate('(f) \n', xy=(x1, ax6.get_ylim()[1]),annotation_clip=False)
ax6.grid(which='major', axis='both', color='lightgray')

fig.suptitle('Convection-permitting (3 km)')
fig.tight_layout()
plt.show()
