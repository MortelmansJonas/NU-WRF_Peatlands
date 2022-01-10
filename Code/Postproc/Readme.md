------------------
POSTPROC
------------------

Notes (10/01/2022):

- This file contains all files used in the postprocessing of the WRF output to calculate and evaluate the lightning parameterization schemes:

	- 'WRF_output_2015_d01_Goddard.py' pis the first file in the postprocessing chain. It takes all wrf output files of a specific simulation run (1 specific year, resolution, microphysics scheme) and puts all desired variables in a single file.

	- These files are created for each year, resolution and microphysics scheme and are then combined for all years per microphysics scheme and resolution 

	- The folder 'Regridding' contains the files used to regrid the model output at convection-permitting resolution to convection-parameterized resolution. This is needed before calculating the diurnal cycle, seasonality and performance diagram.

	- The scripts in 'Calibration' are then used to adjust the model output to the observations and to calculate the PSS values/evaluate the adjustment.

	- The folders 'Spatial_pattern', 'Diurnal_cycle', 'Seasonality', and 'Performance_diagram' contain all scripts used in the calculation of the spatial pattern, diurnal cycle, seasonality, and the performance diagram.

	- Lastly, the 'Input_maps' folder contains two files used to extract and plot the input data.
