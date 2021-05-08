# ---------------------------------------------------------------------------------------------
# MODULES
# ---------------------------------------------------------------------------------------------
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset

# ---------------------------------------------------------------------------------------------
# LOAD DATA
# ---------------------------------------------------------------------------------------------
ds_d01 = Dataset('/scratch/leuven/336/vsc33651/nu-wrf-dev/domain1_diurnal.nc', 'r')
ds_d02 = Dataset('/scratch/leuven/336/vsc33651/nu-wrf-dev/domain2_diurnal.nc', 'r')

Hours = np.linspace(0,24, num=24, endpoint=False)

# ---------------------------------------------------------------------------------------------
# NORMALIZE DATA
# ---------------------------------------------------------------------------------------------
LPI_d01 = (ds_d01['LPI'][:])/(np.nanmax(ds_d01['LPI'][:]))
CAPExP_d01 = (ds_d01['CAPExP'][:])/(np.nanmax(ds_d01['CAPExP'][:]))
PR92W_d01 = (ds_d01['PR92W'][:])/(np.nanmax(ds_d01['PR92W'][:]))
LTG3_d01 = (ds_d01['LTG3'][:])/(np.nanmax(ds_d01['LTG3'][:]))
Obs_d01 = (ds_d01['Obs'][:])/(np.nanmax(ds_d01['Obs'][:]))

LPI_d02 = (ds_d02['LPI'][:])/(np.nanmax(ds_d02['LPI'][:]))
PR92W_d02 = (ds_d02['PR92W'][:])/(np.nanmax(ds_d02['PR92W'][:]))
LTG3_d02 = (ds_d02['LTG3'][:])/(np.nanmax(ds_d02['LTG3'][:]))
Obs_d02 = (ds_d02['Obs'][:])/(np.nanmax(ds_d02['Obs'][:]))

# ---------------------------------------------------------------------------------------------
# PLOT D01
# ---------------------------------------------------------------------------------------------
x1 = 0.0005

fig = plt.figure()
ax1 = plt.subplot2grid((3,2),(0,0))
ax2 = plt.subplot2grid((3,2),(2,0))
ax3 = plt.subplot2grid((3,2),(0,1))
ax4 = plt.subplot2grid((3,2),(1,0))
ax5 = plt.subplot2grid((3,2),(1,1))

ax1.plot(Hours, Obs_d01, 'grey', label='CLDN')
ax1.set_title('CLDN')
ax1.set_ylim(0,1.1)
ax1.annotate('(a) \n ', xy=(x1, ax1.get_ylim()[1]),annotation_clip=False)
ax1.set_ylabel('Max-normalized\nhourly flash density (-)\n*Max-normalized\nhourly LPI (-)')
ax1.grid(which='major', axis='both', color='lightgray')

ax2.plot(Hours, Obs_d01, 'grey', label='CLDN')
ax2.plot(Hours, CAPExP_d01, 'k', label='CAPExP')
ax2.set_title('CAPExP')
ax2.set_ylim(0,1.1)
ax2.annotate('(e) \n ', xy=(x1, ax2.get_ylim()[1]),annotation_clip=False)
ax2.set_ylabel('Max-normalized\nhourly flash density (-)')
ax2.grid(which='major', axis='both', color='lightgray')
ax2.set_xlabel('Hour of the day')

ax3.plot(Hours, Obs_d01, 'grey', label='CLDN')
ax3.plot(Hours, LPI_d01, 'k', label='LPI')
ax3.set_title('LPI*')
ax3.set_ylim(0,1.1)
ax3.annotate('(b) \n ', xy=(x1, ax3.get_ylim()[1]),annotation_clip=False)
ax3.grid(which='major', axis='both', color='lightgray')

ax4.plot(Hours, Obs_d01, 'grey', label='CLDN')
ax4.plot(Hours, LTG3_d01, 'k', label='LT3')
ax4.set_title('LT3')
ax4.set_ylim(0,1.1)
ax4.annotate('(c) \n ', xy=(x1, ax4.get_ylim()[1]),annotation_clip=False)
ax4.grid(which='major', axis='both', color='lightgray')
ax4.set_ylabel('Max-normalized\nhourly flash density (-)')

ax5.plot(Hours, Obs_d01, 'grey', label='CLDN')
ax5.plot(Hours, PR92W_d01, 'k', label='PR92W')
ax5.set_title('PR92W')
ax5.set_ylim(0,1.1)
ax5.set_xlabel('Hour of the day')
ax5.annotate('(d) \n ', xy=(x1, ax6.get_ylim()[1]),annotation_clip=False)
ax5.grid(which='major', axis='both', color='lightgray')

fig.suptitle('Convection-parameterized (9 km)')
fig.tight_layout()
plt.show()

# ---------------------------------------------------------------------------------------------
# PLOT D02
# ---------------------------------------------------------------------------------------------
fig = plt.figure()
ax1 = plt.subplot2grid((3,2),(0,0))
ax2 = plt.subplot2grid((3,2),(0,1))
ax3 = plt.subplot2grid((3,2),(1,0))
ax4 = plt.subplot2grid((3,2),(1,1))

ax1.plot(Hours, Obs_d02, 'grey', label='CLDN')
ax1.set_title('CLDN')
ax1.set_ylim(0,1.1)
ax1.annotate('(a) \n ', xy=(x1, ax1.get_ylim()[1]),annotation_clip=False)
ax1.set_ylabel('Max-normalized\nhourly flash density (-)  \n *Max-normalized\nhourly LPI (-)')
ax1.grid(which='major', axis='both', color='lightgray')

ax2.plot(Hours, Obs_d02, 'grey', label='CLDN')
ax2.plot(Hours, LPI_d02, 'k', label='LPI')
ax2.set_title('LPI*')
ax2.set_ylim(0,1.1)
ax2.annotate('(b) \n ', xy=(x1, ax2.get_ylim()[1]),annotation_clip=False)
ax2.grid(which='major', axis='both', color='lightgray')

ax3.plot(Hours, Obs_d02, 'grey', label='CLDN')
ax3.plot(Hours, LTG3_d02, 'k', label='LT3')
ax3.set_title('LT3')
ax3.set_ylim(0,1.1)
ax3.annotate('(c) \n ', xy=(x1, ax3.get_ylim()[1]),annotation_clip=False)
ax3.set_ylabel('Max-normalized\nhourly flash density (-)')
ax3.grid(which='major', axis='both', color='lightgray')

ax4.plot(Hours, Obs_d02, 'grey', label='CLDN')
ax4.plot(Hours, PR92W_d02, 'k', label='PR92W')
ax4.set_title('PR92W')
ax4.set_ylim(0,1.1)
ax4.set_xlabel('Hour of the day')
ax4.annotate('(d) \n ', xy=(x1, ax4.get_ylim()[1]),annotation_clip=False)
ax4.grid(which='major', axis='both', color='lightgray')

fig.suptitle('Convection-permitting (3 km)')
fig.tight_layout()
plt.show()
