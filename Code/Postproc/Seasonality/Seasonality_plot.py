# ---------------------------------------------------------------------------------------------
# MODULES
# ---------------------------------------------------------------------------------------------
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from netCDF4 import Dataset

# ---------------------------------------------------------------------------------------------
# LOAD DATA
# ---------------------------------------------------------------------------------------------

ds_d01 = Dataset('/scratch/leuven/336/vsc33651/nu-wrf-dev/domain1_seasonality.nc', 'r')
ds_d02 = Dataset('/scratch/leuven/336/vsc33651/nu-wrf-dev/domain2_seasonality.nc', 'r')

# ---------------------------------------------------------------------------------------------
# SEASONALITY 2.0
# ---------------------------------------------------------------------------------------------
weeks = pd.date_range('2015-06-01','2015-08-30',freq='2W')

biweekly_LPI_d02 = np.zeros((7))
biweekly_LTG3_d02 = np.zeros((7))
biweekly_PR92W_d02 = np.zeros((7))
biweekly_Obs_d02 = np.zeros((7))
biweekly_LPI_d01 = np.zeros((7))
biweekly_LTG3_d01 = np.zeros((7))
biweekly_PR92W_d01 = np.zeros((7))
biweekly_CAPExP_d01 = np.zeros((7))
biweekly_Obs_d01 = np.zeros((7))

for i in range(0,7):
    j = (i+1)*2
    k = i*2
    biweekly_LPI_d02[i] = np.nanmean(ds_d02['LPI'][k:j])
    biweekly_LTG3_d02[i] = np.nanmean(ds_d02['LTG3'][k:j])
    biweekly_PR92W_d02[i] = np.nanmean(ds_d02['PR92W'][k:j])
    biweekly_Obs_d02[i] = np.nanmean(ds_d02['Obs'][k:j])
    biweekly_LPI_d01[i] = np.nanmean(ds_d01['LPI'][k:j])
    biweekly_LTG3_d01[i] = np.nanmean(ds_d01['LTG3'][k:j])
    biweekly_PR92W_d01[i] = np.nanmean(ds_d01['PR92W'][k:j])
    biweekly_CAPExP_d01[i] = np.nanmean(ds_d01['CAPExP'][k:j])
    biweekly_Obs_d01[i] = np.nanmean(ds_d01['Obs'][k:j])

# NORMALIZE DATA
LPI_d01 = (biweekly_LPI_d01[:])/(np.nanmax(biweekly_LPI_d01[:]))
CAPExP_d01 = (biweekly_CAPExP_d01[:])/(np.nanmax(biweekly_CAPExP_d01[:]))
PR92W_d01 = (biweekly_PR92W_d01[:])/(np.nanmax(biweekly_PR92W_d01[:]))
LTG3_d01 = (biweekly_LTG3_d01[:])/(np.nanmax(biweekly_LTG3_d01[:]))
Obs_d01 = (biweekly_Obs_d01[:])/(np.nanmax(biweekly_Obs_d01[:]))

LPI_d02 = (biweekly_LPI_d02[:])/(np.nanmax(biweekly_LPI_d02[:]))
PR92W_d02 = (biweekly_PR92W_d02[:])/(np.nanmax(biweekly_PR92W_d02[:]))
LTG3_d02 = (biweekly_LTG3_d02[:])/(np.nanmax(biweekly_LTG3_d02[:]))
Obs_d02 = (biweekly_Obs_d02[:])/(np.nanmax(biweekly_Obs_d02[:]))

# DOMAIN 1
x1 = 0.002

fig = plt.figure()
ax1 = plt.subplot2grid((3,2),(0,0))
ax2 = plt.subplot2grid((3,2),(2,0))
ax3 = plt.subplot2grid((3,2),(0,1))
ax4 = plt.subplot2grid((3,2),(1,0))
ax6 = plt.subplot2grid((3,2),(1,1))

ax1.plot(weeks, Obs_d01, 'grey', label='CLDN')
ax1.set_title('CLDN')
ax1.set_ylim(0,1.1)
ax1.grid(which='major', axis='both', color='lightgray')
ax1.set_xticks(ticks=[pd.to_datetime('2015-06', format='%Y-%m'), pd.to_datetime('2015-07', format='%Y-%m'),
                     pd.to_datetime('2015-08', format='%Y-%m'), pd.to_datetime('2015-09', format='%Y-%m')])
ax1.set_xticklabels(['Jun', 'Jul', 'Aug', 'Sep'])
ax1.annotate('(a) \n ', xy=(x1, ax1.get_ylim()[1]),annotation_clip=False)
ax1.set_ylabel('Max-normalized\nbiweekly flash density (-) \n *Max-normalized\nbiweekly LPI (-)')

ax2.plot(weeks, Obs_d01, 'grey',label='CLDN')
ax2.plot(weeks, CAPExP_d01, 'k', label='CAPExP')
ax2.set_title('CAPExP')
ax2.set_ylim(0,1.1)
ax2.grid(which='major', axis='both', color='lightgray')
ax2.set_xticks(ticks=[pd.to_datetime('2015-06', format='%Y-%m'), pd.to_datetime('2015-07', format='%Y-%m'),
                     pd.to_datetime('2015-08', format='%Y-%m'), pd.to_datetime('2015-09', format='%Y-%m')])
ax2.set_xticklabels(['Jun', 'Jul', 'Aug', 'Sep'])
ax2.annotate('(e) \n ', xy=(x1, ax2.get_ylim()[1]),annotation_clip=False)
ax2.set_ylabel('Max-normalized\nbiweekly flash density (-)')

ax3.plot(weeks, Obs_d01, 'grey',label='CLDN')
ax3.plot(weeks, LPI_d01, 'k', label='LPI')
ax3.set_title('LPI*')
ax3.set_ylim(0,1.1)
ax3.grid(which='major', axis='both', color='lightgray')
ax3.set_xticks(ticks=[pd.to_datetime('2015-06', format='%Y-%m'), pd.to_datetime('2015-07', format='%Y-%m'),
                     pd.to_datetime('2015-08', format='%Y-%m'), pd.to_datetime('2015-09', format='%Y-%m')])
ax3.set_xticklabels(['Jun', 'Jul', 'Aug', 'Sep'])
ax3.annotate('(b) \n ', xy=(x1, ax3.get_ylim()[1]),annotation_clip=False)

ax4.plot(weeks, Obs_d01, 'grey',label='CLDN')
ax4.plot(weeks, LTG3_d01, 'k', label='LT3')
ax4.set_title('LT3')
ax4.set_ylim(0,1.1)
ax4.grid(which='major', axis='both', color='lightgray')
ax4.set_xticks(ticks=[pd.to_datetime('2015-06', format='%Y-%m'), pd.to_datetime('2015-07', format='%Y-%m'),
                     pd.to_datetime('2015-08', format='%Y-%m'), pd.to_datetime('2015-09', format='%Y-%m')])
ax4.set_xticklabels(['Jun', 'Jul', 'Aug', 'Sep'])
ax4.annotate('(c) \n ', xy=(x1, ax4.get_ylim()[1]),annotation_clip=False)
ax4.set_ylabel('Max-normalized\nbiweekly flash density (-)')

ax6.plot(weeks, Obs_d01, 'grey',label='CLDN')
ax6.plot(weeks, PR92W_d01, 'k', label='PR92W')
ax6.set_title('PR92W')
ax6.set_ylim(0,1.1)
ax6.grid(which='major', axis='both', color='lightgray')
ax6.set_xticks(ticks=[pd.to_datetime('2015-06', format='%Y-%m'), pd.to_datetime('2015-07', format='%Y-%m'),
                     pd.to_datetime('2015-08', format='%Y-%m'), pd.to_datetime('2015-09', format='%Y-%m')])
ax6.set_xticklabels(['Jun', 'Jul', 'Aug', 'Sep'])
ax6.annotate('(d) \n ', xy=(x1, ax6.get_ylim()[1]),annotation_clip=False)

fig.suptitle('Convection-parameterized (9 km)')
fig.tight_layout()
plt.show()

# DOMAIN 2
fig = plt.figure()
ax1 = plt.subplot2grid((2,2),(0,0))
ax3 = plt.subplot2grid((2,2),(0,1))
ax4 = plt.subplot2grid((2,2),(1,0))
ax6 = plt.subplot2grid((2,2),(1,1))

ax1.plot(weeks, Obs_d02, 'grey', label='CLDN')
ax1.set_title('CLDN')
ax1.set_ylim(0,1.1)
ax1.grid(which='major', axis='both', color='lightgray')
ax1.set_xticks(ticks=[pd.to_datetime('2015-06', format='%Y-%m'), pd.to_datetime('2015-07', format='%Y-%m'),
                     pd.to_datetime('2015-08', format='%Y-%m'), pd.to_datetime('2015-09', format='%Y-%m')])
ax1.set_xticklabels(['Jun', 'Jul', 'Aug', 'Sep'])
ax1.annotate('(a) \n ', xy=(x1, ax1.get_ylim()[1]),annotation_clip=False)
ax1.set_ylabel('Max-normalized\nbiweekly flash density (-) \n *Max-normalized\nbiweekly LPI (-)')

ax3.plot(weeks, Obs_d02, 'grey',label='CLDN')
ax3.plot(weeks, LPI_d02, 'k', label='LPI')
ax3.set_title('LPI*')
ax3.set_ylim(0,1.1)
ax3.grid(which='major', axis='both', color='lightgray')
ax3.set_xticks(ticks=[pd.to_datetime('2015-06', format='%Y-%m'), pd.to_datetime('2015-07', format='%Y-%m'),
                     pd.to_datetime('2015-08', format='%Y-%m'), pd.to_datetime('2015-09', format='%Y-%m')])
ax3.set_xticklabels(['Jun', 'Jul', 'Aug', 'Sep'])
ax3.annotate('(b) \n ', xy=(x1, ax3.get_ylim()[1]),annotation_clip=False)

ax4.plot(weeks, Obs_d02, 'grey',label='CLDN')
ax4.plot(weeks, LTG3_d02, 'k', label='LT3')
ax4.set_title('LT3')
ax4.set_ylim(0,1.1)
ax4.grid(which='major', axis='both', color='lightgray')
ax4.set_xticks(ticks=[pd.to_datetime('2015-06', format='%Y-%m'), pd.to_datetime('2015-07', format='%Y-%m'),
                     pd.to_datetime('2015-08', format='%Y-%m'), pd.to_datetime('2015-09', format='%Y-%m')])
ax4.set_xticklabels(['Jun', 'Jul', 'Aug', 'Sep'])
ax4.annotate('(c) \n ', xy=(x1, ax4.get_ylim()[1]),annotation_clip=False)
ax4.set_ylabel('Max-normalized\nbiweekly flash density (-)')

ax6.plot(weeks, Obs_d02, 'grey',label='CLDN')
ax6.plot(weeks, PR92W_d02, 'k', label='PR92W')
ax6.set_title('PR92W')
ax6.set_ylim(0,1.1)
ax6.grid(which='major', axis='both', color='lightgray')
ax6.set_xticks(ticks=[pd.to_datetime('2015-06', format='%Y-%m'), pd.to_datetime('2015-07', format='%Y-%m'),
                     pd.to_datetime('2015-08', format='%Y-%m'), pd.to_datetime('2015-09', format='%Y-%m')])
ax6.set_xticklabels(['Jun', 'Jul', 'Aug', 'Sep'])
ax6.annotate('(d) \n ', xy=(x1, ax6.get_ylim()[1]),annotation_clip=False)

fig.suptitle('Convection-permitting (3 km)')
fig.tight_layout()
plt.show()
