# SPEI Calculator for netCDF files

import datetime as dt
from dateutil.relativedelta import relativedelta
import numpy as np
import os
from netCDF4 import Dataset
import pickle
import sharppy
from sharppy.sharptab import profile, fire, thermo, utils

# method to generate time bounds given a starting month and how many time steps
def timeBounds(months, startingMonth=1):
    """ Takes in amount of months we need timebounds for
        also starting month, Jan=1, Feb=2 etc.
        returns a 2d numpy array of the first and last day in the month, with the starting month having 0 as the first term
    """
    # amount of days in each month starting in January
    steps = [31.0, 28.0, 31.0, 30.0, 31.0, 30.0, 31.0, 31.0, 30.0, 31.0, 30.0, 31.0]

    # initialize the timeBounds array
    bounds = np.zeros((months, 2))
    start = 0.0

    for m in range(months):
        currentMonthIndex = ((m + startingMonth) % 12) - 1
        end = start + steps[currentMonthIndex]
        mArray = np.array([start, end])
        bounds[m] = mArray
        start = end
    
    return bounds

# encoding = 'latin1' for precl
def fromPickle(fname, encoding=None):
    """ Gets a file from pickle
    """
    with open(fname, 'rb') as f:
        # The protocol version used is detected automatically, so we do not
        # have to specify it.
        if encoding == None:
            data = pickle.load(f)
        else:
            print(f)
            data = pickle.load(f, encoding=encoding)
    return data

# Do 'Nonintervention' for nonintervention data
# For this to work you have to be in the correct directory
# I also renamed the files just by their variable name so T, RELHUM, PRECC, etc..
# NOTE: The global and fram intervention data both start in may of 2000, control data starts in Jan of 2000
def pickleToArrays(interv='Control'):
    d = {}
    if interv.lower() == 'fram':
        intervType = 'Fram'
    elif interv.lower() == 'global':
        intervType = 'Global'
    elif interv.lower() == 'control':
        intervType = 'Nonintervention'
    else:
        print('Need to input a valid intervention type: fram, global, or control')
        return

    baseDir = 'Fire Weather Data/' + intervType
    for fname in os.listdir('Fire Weather Data/' + intervType):
        if not fname == '.DS_Store':
            dataList = fromPickle(baseDir + '/' +  fname, encoding='latin1')
            d[fname[:-7]] = np.array(dataList)

    # PRECT is the sum of PRECC and PRECL
    prect = d['PRECC'] + d['PRECL']
    d['PRECT'] = prect

    return d

def calcFWI(temp, rh, wspd):
    """ 
    Params:
        temp: temperature - in K
        rh: relative humidity - in percent
        wspd: wind speed - in m/s

    Returns:
        fwi value based on equations found at https://a.atmos.washington.edu/wrfrt/descript/definitions/fosbergindex.html
    """
    temp = thermo.ktof(temp)
    wspd = utils.MS2MPH(wspd)

    # rh = thermo.relh(prof.pres[prof.sfc], prof.tmpc[prof.sfc], prof.dwpc[prof.sfc])
    if (rh <= 10):
        em = 0.03229 + 0.281073*rh - 0.000578*rh*temp
    elif (10 < rh <= 50):
        em = 2.22749 + 0.160107*rh - 0.014784*temp
    else:
        em = 21.0606 + 0.005565*rh*rh - 0.00035*rh*temp - 0.483199*rh

    em30 = em/30
    w_sq = wspd**2
    mdc = 1 - 2*em30 + 1.5*((em30)**2) - 0.5*((em30)**3)

    fwi = (mdc*(np.sqrt(1+w_sq)))/0.3002

    return fwi

def fwi(fname, interv='Control'):
    """ Calclates Fosberg FWI for each (time, lat, lon) data point and creates a new np array.
        Writes this array to the specified NC File.
        
        The Fosberg Fire Weather Index (FFWI), also known as the Fire Weather Index (FWI), is a fire weather index created to measure the potential influence of weather on a wildfire based on model output of temperature, wind and relative humidity.

        The index represents expected flame length and fuel drying based on model output fields of temperature, wind and humidity. Large values of the FWI imply high flame lengths and rapid drying. The values were designed so that a index rating of 100 is equal to a moisture content of 0 and a wind speed of 30 mph. Larger combinations of these values still result in an index of 100.
    """
    # Global Intervention, Fram Intervention, or Nonintervention data

    # Getting the np array data Temp, RelH, and Wind
    Arrs = pickleToArrays(interv)

    # Get the necessary arrays
    T = Arrs['T']           # Temperature
    RH = Arrs['RELHUM']     # Relative Humidity
    WS = Arrs['U10']        # Wind Speed
    
    time,lat,lon = Arrs['T'].shape   # This assumes that they all have the same shape

    FWI = np.zeros((time, lat, lon))

    for t in range(time):
        if t%100 == 0:
            print(t)
        for i in range(lat):
            for j in range(lon):
                # Get the temp, relh, and wspd values
                temp = T[t][i][j]
                relh = RH[t][i][j]
                wspd = WS[t][i][j]

                # Calculate FWI given the values
                val = calcFWI(temp, relh, wspd)
                
                # put these values into the final arrays
                FWI[t][i][j] = val

    #Store this data in an nc file given by the input filenames
    data = Dataset(fname, 'r+')

    variables = data.variables

    dtype = np.dtype('float32')
    dims = ('time', 'lat', 'lon')
    description = 'Fosberg Fire Weather Index'

    # Initialize the FWI variable if necessary
    initialize = False
    if 'FWI' not in variables.keys():
        data.createVariable('FWI', dtype, dims)
        initialize = True

    fwiVar = variables['FWI']
    
    # if this is the first time creating it then 
    if initialize:
        fwiVar.long_name = description

    for t in range(time):
        fwiVar[t] = FWI[t]

    return FWI


def yearlyAvgJan(monthly, startingMonth):
    """ Converts from a monthly data to a nparray of m-month averaged data (so yearly data when m=12)
        Adjust the data to start the years in January
    """
    time,lat,lon = monthly.shape

    # How many months you need to ignore at the begging to get the first year starting in January
    startIndex = ((12-startingMonth) + 1) % 12

    # How many years after we get rid of the excess original months
    years = (time-startIndex) // 12

    # How many extra months at the end we need to ignore
    correction = (time-startIndex) % 12

    print('years:', years)
    yearly = np.zeros((years, lat, lon))

    count = 0
    for i in range(lat):
        for j in range(lon):
            if count%1000==0:
                print(count)
            count += 1
            series = np.copy(monthly[:,i,j])[startIndex:time-correction]

            # splits the values up into years , so 12 months per array
            series = series.reshape((years, 12))

            average = np.nanmean(series, axis=1)

            yearly[:,i,j] = average

    return yearly


def yearlyFWI(
        fname="framFWI.nc",
        newVar='FWIYearlyAvg',
        description='Yearly Average of FWI (2001 - 2081)',
        startingMonth=5
    ):
    """
    Paramaters:
        fname: filename path
        newVar: variable name
        description: short description of the variable
        startingMonth: int - month your data starts

    Return:
        3d array of yearly fwi data - averages fwi over every year
        saves this data to the specified nc file (have not implemented this yet)
    """
     # reading netCDF dataset into variable and able to overwrite the data as well with 'r+'
    data = Dataset(fname, 'r+')

    # Get the SPI variable
    variables = data.variables
    fwi = variables['FWI']

    # Call the helper function to calculate the yearly average
    fwiYearly = yearlyAvgJan(fwi, startingMonth)

    years,_,_= fwiYearly.shape

    return fwiYearly
