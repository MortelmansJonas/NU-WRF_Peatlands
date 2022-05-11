#!/usr/bin/env python
## Load Modules
from netCDF4 import Dataset
import netCDF4 as nc


def create_file_from_source(src_file, trg_file):
    src = nc.Dataset(src_file)
    trg = nc.Dataset(trg_file, mode='w')

    # Create the dimensions of the file
    for name, dim in src.dimensions.items():
        trg.createDimension(name, len(dim) if not dim.isunlimited() else None)

    # Copy the global attributes
    trg.setncatts({a:src.getncattr(a) for a in src.ncattrs()})

    # Create the variables in the file
    for name, var in src.variables.items():
        trg.createVariable(name, var.dtype, var.dimensions)

        # Copy the variable attributes
        trg.variables[name].setncatts({a:var.getncattr(a) for a in var.ncattrs()})

        # Copy the variables values (as 'f4' eventually)
        trg.variables[name][:] = src.variables[name][:]

    # Save the file
    trg.close()


infile = '/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/NU-WRF/wrfout_nc_files/d01_all.nc'
ds01 = Dataset(infile, 'r')
outfile = '/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/NU-WRF/wrfout_nc_files/domain2_at_domain1_all.nc'
create_file_from_source(infile,outfile)
ds_out = nc.Dataset(outfile,"a")

LI = 'LTG3'
infile = '/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/NU-WRF/wrfout_nc_files/domain2_at_domain1_'+LI+'.nc'
ds_in = Dataset(infile, 'r')
ds_out[LI][:,:,:] = ds_in[LI][:,:,:]

LI = 'CAPExP_R'
infile = '/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/NU-WRF/wrfout_nc_files/domain2_at_domain1_'+LI+'.nc'
ds_in = Dataset(infile, 'r')
ds_out[LI][:,:,:] = ds_in[LI][:,:,:]

LI = 'PR92W'
infile = '/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/NU-WRF/wrfout_nc_files/domain2_at_domain1_'+LI+'.nc'
ds_in = Dataset(infile, 'r')
ds_out[LI][:,:,:] = ds_in[LI][:,:,:]

LI = 'LPI'
infile = '/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/NU-WRF/wrfout_nc_files/domain2_at_domain1_'+LI+'.nc'
ds_in = Dataset(infile, 'r')
ds_out[LI][:,:,:] = ds_in[LI][:,:,:]

ds_out.close()

