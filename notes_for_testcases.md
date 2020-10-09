#Notes on running testcases of nu-wrf on HPC of KU Leuven\

#Src code of nu-wrf-dev-noelec (group has read access):\
ls /data/leuven/317/vsc31786/src_code/nu-wrf-dev-noelec\

#Path to nu-wrf input, mix of things from staging (linked),\
#files from Isis (linked), and newly added files (e.g. topo_30s).\ 
#Use this path to add new data (group with read/write access)\
#download files from ucar wrf page (e.g. done for topo_30s):\
#https://www2.mmm.ucar.edu/wrf/users/download/get_sources_wps_geog.html\
ls /staging/leuven/stg_00024/OUTPUT/michelb/l_data/model_param/nu-wrf

#Project folder (read access only, copy to your data to work on testcases)\
ls /data/leuven/317/vsc31786/projects/nu-wrf-dev-noelec/testcases/wrflis/noahmp36_modis_merra2

#copy with rsync. This preserves file attributes like permissions and executable yes/no, etc. :\
rsync -av /data/leuven/317/vsc31786/projects/nu-wrf-dev-noelec $VSC_DATA/projects/.

#To run code, use your folder on scratch (I made my one readable for you\
#rsync from my data to your scratch\ 
#(in the variable $VSC_SCRATCH is already your scratch path defined, convenient to use)\
rsync -av /data/leuven/317/vsc31786/projects/nu-wrf-dev-noelec $VSC_SCRATCH/.

#go to this folder\
cd $VSC_SCRATCH/nu-wrf-dev-noelec/testcases/wrflis/noahmp36_modis_merra2

#load modules and path definitions\
source common.reg 

#The region check shows with the next command\
ncl $NUWRFDIR/WPS/util/plotgrids_new.ncl

#run geogrid.reg, it uses input from the namelist.wps\
./geogrid.reg

#the current error message asks for the modis_landuse (see below).\ 
#From the error message, it's clear which files you need to download from the above webpage.\ 
#As said before, put it into my /staging/leuven/stg_00024/OUTPUT/michelb/l_data/model_param/nu-wrf\
#the input files are stored in two folders, LS_PARAMETERS (for LIS) and geog (needed for geogrid)\
#While the error message is clear, it's not fully clear to me a priori which files are needed.\
#It is somehow contained in this file geogrid/GEOGRID.TBL\
#Don't know yet this format...\
vsc31786@tier2-p-login-3:/scratch/leuven/317/vsc31786/nu-wrf-dev-noelec/testcases/wrflis/noahmp36_modis_merra2$ ./geogrid.reg\ 
ERROR: Could not open /staging/leuven/stg_00024/OUTPUT/michelb/l_data/model_param/nu-wrf/geog/modis_landuse_21class_30s/index

#Group storage containing forcing data, static input data, satellite data, etc.\
#Worth to check what is available already\
/staging/leuven/stg_00024/

NOTE: scratch folder are removed after 1 month of inactivity. So relevant output should be saved at some point on your data disk.
