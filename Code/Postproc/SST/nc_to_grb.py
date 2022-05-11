# ---------------------------------------------------------------------------------------------
# MODULES (CHANGE INTERPRETER TO py3_iris:
# ---------------------------------------------------------------------------------------------
import iris
import iris_grib
from iris import coord_systems
from iris.coords import DimCoord, AuxCoord
import os
import glob
import pandas as pd
import numpy as np

os.chdir('/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/SST/')
for file in glob.glob('*.nc'):
    file_open = os.path.join('/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/SST/'+file)
    filename = os.path.join('/scratch/leuven/projects/lt1_2020_es_pilot/project_output/rsda/vsc33651/SST/','SSTRSS:' + file[0:4] + '-' + file[4:6] + '-' + file[6:8] + '_' + file[8:10] + '.grb2')
    data1 = iris.load(file_open, 'sea_surface_foundation_temperature')[0]

    time1_dims = data1.coord_dims('time')[0]
    data1.Originating_or_generating_Center = 7
    data1.level_type = 1
    data1.remove_coord('time')
    date = ((pd.to_datetime(file[0:4] + '-' + file[4:6] + '-' + file[6:8])) - pd.to_datetime('1981-01-01')).days*24
    new_t_coord = iris.coords.DimCoord(date, standard_name='time', units= 'hours since 1981-01-01')
    data1.add_dim_coord(new_t_coord, time1_dims)

    lon_coord1 = data1.coords("longitude")[0]
    lat_coord1 = data1.coords("latitude")[0]
    lon_coord1.coord_system = iris.coord_systems.GeogCS(iris.fileformats.pp.EARTH_RADIUS)
    lat_coord1.coord_system = iris.coord_systems.GeogCS(iris.fileformats.pp.EARTH_RADIUS)

    data1.add_aux_coord(iris.coords.DimCoord(0, standard_name='forecast_period', units='hours'))
    data1.add_aux_coord(iris.coords.DimCoord(0, "height", units="m"))

    data1.var_name = "SST"
    data1.long_name = "Sea Surface Temperature"
    data1.standard_name = 'sea_surface_temperature'
    data1.units = 'K'
    data1.discipline = 0
    data1.category = 0
    data1.typeOfFirstFixedSurface = 1
    time = file[0:4] + '-' + file[4:6] + '-' + file[6:8]
    time1 = file[0:4] + '-' + file[4:6] + '-' + file[6:8]

    print(data1)
    iris_grib.save_grib2(data1, filename)