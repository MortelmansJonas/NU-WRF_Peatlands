------------------
Postproc
------------------

- The folder 'Atmospheric_variables' contains all scripts that are used to compare the modeled precipitation and heat fluxes with the MERRA-2 input data. 

- The folder 'UPP' contains the scripts used to run UPP and example files to combine the UPP output for all years for both domains (only for 2015). The same scripts are used for all the other years by simply changing the dates in the script to the corresponding year.
  
- The files 'WRF_OUTPUT_2015_d01_nc.py' and 'WRF_OUTPUT_2015_d02_nc.py' are used to extract the wanted variables of the WRF hourly output files and combine all the data per year and per domain. To get this for other years, the dates in these scripts need to be changed to the corresponding years.

- The files 'd01_nc.py' and 'd02_nc.py' are used to calculate the indices for each year and store them all in one netcdf file.

- The files 'Daily_nc.py' and 'Daily_nc_regridded.py' contain the script to calculate the daily averages for each index and use the output of 'd01_nc.py' and 'd02_nc.py' as input. The '...regridded' file is specifically used for the regridded d02 data, while 'Daily_nc.py' is used for the data at its original resolution.

- The files 'Weekly_nc.py' and 'Weekly_nc_regridded.py' contain the script to calculate the weekly averages based on the daily averages. The regridded file is again specifically used for the regridded data of d02.
  
- The file 'Spatial_pattern.py' contains the script to map the spatial distribution patterns for each index. To do so, it uses the daily averaged data from 'Daily_nc.py'.

- The file 'Timeseries.py' is used to calculate and plot the 6-year timeseries of the daily flash density for both domains. It used the data from 'Daily_nc.py' and 'Weekly_nc.py'.

- The folder 'Diurnal' contains two scripts. The first one, 'Diurnal_Cycle.py', is used to calculate the diurnal cycle based on the hourly, regridded data and stores this data in a netcdf file. The other file, 'Diurnal_plot.py' then uses this file to make the plots.
  
- The folder 'Seasonality' also contains two scripts. 'Seasonality_nc.py' calculates the seasonality, based on the data from 'weekly_nc_regridded.py' and stores it in a netcdf file. The script 'Seasonality_plot.py' then plots the data stored in this file.
  
- Lastly, the file 'Frequency' uses the regridded hourly data to calculate and plot the frequency distribution per index and per domain.
