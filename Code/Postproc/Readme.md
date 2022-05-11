------------------
POSTPROC
------------------

This folder contains all the files used in the postprocessing of the WRF output to calculate and evaluate the lightning indices as shown in the paper. For additional files (e.g. comparison of seasonality or other analysis that are not shown in the paper, see the Pycharm environment).

11/05/2022:
-----------

- SST
	This folder contains two scripts used as a workaround for the erroneous SST data provided by the original WRF scripts. It contains two files to open and read GRIB data.

- WRF_output
	This folder contains the first files in the postprocessing chain.
	* LPI_nc.py: This script combines all LPI hourly files into one file for a specific NU-WRF run (i.e. one year, one experiment, one domain).
	* wrfout_script.py: This file takes all wrf output files for a specific experiment. domain, and year and combines all hours into one .nc file. This file needs to be run for all experiments, all years, and both domains separately.
	* d01_nc.py and d02_nc.py: These files combine all output files of the previous scripts (LPI and wrfout) into one file. They combines all years for a specific experiment for domain 1 and domain 2, respectively.
	* MERRA2_maps.py: This script is the first script to really produce output as it creates the maps to evaluate the latent and sensible heat fluxes.

- Regridding:
	* Regridding_model_all.py and Regridding_model_merge.py: These files are used to first regrid all data from d02 to the resolution of d01 (best run for each variable separately). The merge file then merges all regridded files into one file that can then be used for further analysis.

- Calibration
	This folder contains two scripts to linearly rescale the model output to the observations and plot the Frequency distributions.

- Input_maps
	* LIS_HIST_comb.py: This file combines LIS_HIST files of all timesteps and extracts the albedo and greenness fraction.
	* Mapping_input.py: This file then maps the input (albedo, greenness fraction, elevation, and LUC)

- Spatial_pattern, Diurnal Cycle, and Performance_diagram:
	These folders contain all files necessary to calculate and evaluate the spatial patterns, diurnal cycle and overall event-based skill assessment, respectively



10/01/2022:
-----------

- This file contains all files used in the postprocessing of the WRF output to calculate and evaluate the lightning parameterization schemes:

	- 'WRF_output_2015_d01_Goddard.py' is the first file in the postprocessing chain. It takes all wrf output files of a specific simulation run (1 specific year, resolution, microphysics scheme) and puts all desired variables in a single file.

	- These files are created for each year, resolution and microphysics scheme and are then combined for all years per microphysics scheme and resolution 

	- The folder 'Regridding' contains the files used to regrid the model output at convection-permitting resolution to convection-parameterized resolution. This is needed before calculating the diurnal cycle, seasonality and performance diagram.

	- The scripts in 'Calibration' are then used to adjust the model output to the observations and to calculate the PSS values/evaluate the adjustment.

	- The folders 'Spatial_pattern', 'Diurnal_cycle', 'Seasonality', and 'Performance_diagram' contain all scripts used in the calculation of the spatial pattern, diurnal cycle, seasonality, and the performance diagram.

	- Lastly, the 'Input_maps' folder contains two files used to extract and plot the input data.
