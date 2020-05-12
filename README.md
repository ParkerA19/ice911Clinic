# ice911Clinic

Each of the 3 main files are standalone files that were used throughout the year to calculate variables and generate animations with data. They are all python files and can be run in an interactive python environment such as ipython. In general these files demonstrate how to load in netCDF files to access the climate data, how to loop through this data to calculate new variables, and how to write this new data to netCDF files.

# Spi.py
This file was used initally for calculating SPI given precipitation data. It was modeled after the Climate Indices implementation in python, however after some time it was decided to just use the pre-built command line tool that Climate Indices provides. Setup information and documentation on this tool is available on the [website](https://climate-indices.readthedocs.io/en/latest/). There are some other functions in this file that allow yearly averages to be taken of data which should be useful and the generalizable to many different data variables. 


# Fireweather.py

# Gif.py
