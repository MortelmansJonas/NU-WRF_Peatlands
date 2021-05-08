------------------
Observational Data
------------------

Notes (08/05/2021):

- 'Lightning_data_0.1d.py' does multiple things:

	- It plots the amount of flashes per day and per type (CC and CG)
	- It plots the amount of flashes per day, per type and per location (latitudes and longitudes are rounded to 1 decimal places.)
	- It creates a netCDF file with, for each type separately, the flashedensity and the average Estimated Peak Current (EPC).
	
- 'Lightning_data_0.01d.py' does the same thing as 'Lightning_data_0.1d.py', but rounds the latitudes and longitudes to 2 decimal places.
