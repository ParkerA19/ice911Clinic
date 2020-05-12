# ice911Clinic

Each of the 3 main files are standalone files that were used throughout the year to calculate variables and generate animations with data. They are all python files and can be run in an interactive python environment such as ipython. In general these files demonstrate how to load in netCDF files to access the climate data, how to loop through this data to calculate new variables, and how to write this new data to netCDF files. Each file has functions which are well documented and explain the inputs and outputs.

# Spi.py
This file was used initally for calculating SPI given precipitation data. It was modeled after the Climate Indices implementation in python, however after some time it was decided to just use the pre-built command line tool that Climate Indices provides. Setup information and documentation on this tool is available on the [website](https://climate-indices.readthedocs.io/en/latest/). There are some other functions in this file that allow yearly averages to be taken of data which should be useful and the generalizable to many different data variables. 


# Fireweather.py
This file was used to calculate the [Fosberg Fire Weather Index](https://a.atmos.washington.edu/wrfrt/descript/definitions/fosbergindex.html). This index requires temperature, relative humidity, and wind speed data. In order to store this data, we used pickle files and one of the functions allowed us to store a dictionary of each of these files. However, it does require a specific folder structure as shown here ![folders]{/folders.png}.

# Gif.py
This is just a simple python file that concatenates all the image files in a given folder into a gif. This was mainly used in accordance with Panoply as we were able to export a series of png images, and then pass these images into our gif function and animate the data.
