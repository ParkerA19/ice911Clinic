# SPEI Calculator for netCDF files
# Just use spi()

import numpy as np
import os
from netCDF4 import Dataset
import pickle
import math
import multiprocessing as mp
import scipy.stats
from scipy.stats import gamma

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
            data = pickle.load(f, encoding=encoding)
    return data

def gamma_parameters(
        values: np.ndarray,
) -> (np.ndarray, np.ndarray):
    """
    Computes the gamma distribution parameters alpha and beta.

    :param values: 2-D array of values, with each row typically representing a year
                   containing twelve columns representing the respective calendar
                   months, or 366 days per column as if all years were leap years
    :return: two 2-D arrays of gamma fitting parameter values, corresponding in size
        and shape of the input array
    :rtype: tuple of two 2-D numpy.ndarrays of floats, alphas and betas
    """

    # replace zeros with NaNs
    values[values == 0] = np.NaN

    calibration_values = values

    # compute the gamma distribution's shape and scale parameters, alpha and beta
    # TODO explain this better
    means = np.nanmean(calibration_values, axis=0)
    log_means = np.log(means)
    logs = np.log(calibration_values)
    mean_logs = np.nanmean(logs, axis=0)
    a = log_means - mean_logs
    alphas = (1 + np.sqrt(1 + 4 * a / 3)) / (4 * a)
    betas = means / alphas

    # print('alphas: ', alphas, 'beta:', betas)
    return alphas, betas

def computeSigmas(values, hist_values):

    zeros = (values == 0).sum(axis=0)
    # print('zeros:', zeros)
    probabilities_of_zero = zeros / values.shape[0]
    # print('z prob:', probabilities_of_zero)
    alphas, betas = gamma_parameters(hist_values)
    gprobs = gamma.cdf(values, a=alphas, scale=betas)
    # print('gprobs: ', gprobs[0])

    # (normalize including the probability of zero, putting into the range [0..1]?)
    probabilities = probabilities_of_zero + ((1 - probabilities_of_zero) * gprobs)
    # print('probs:', prob   abilities[0])

    return scipy.stats.norm.ppf(probabilities)

def calculate_over_full_series_with_distribution(data, hist_data, starting_month, calibration_years, option=True):
    """ Takes in 3d np array of data and hist_data
        Fits a gamma distribution based on hist data
        compares the data to this distribution to get the standard deviations for each point
        returns a 3d array of spi values
    """

    # Number of months in data
    n_months = np.shape(data)[0]
    nh_months = np.shape(hist_data)[0]

    # Pre-allocate SPI
    spi = np.zeros(np.shape(data))*np.nan

    # Single date series
    if data.ndim == 1:
        # print('data ndim is 1')
        data = data.reshape(len(data), 1)
        spi = spi.reshape(len(spi), 1)

    if hist_data.ndim == 1:
        hist_data = hist_data.reshape(len(hist_data),1)


    count = 0
    for i in range(np.shape(data)[1]):          # Loop over all latitude
        for k in range(np.shape(data)[2]):      # Loop over all longitude
            if count%1000==100:
                print(count)
            count += 1

            # Make a series of numbers (at a single lat long and each is at a different time step)
            data_one_series = np.copy(data[:,i,k])
            hist_data_one_series = np.copy(hist_data[:,i,k])

            completeYears = n_months // 12
            extraMonths = n_months % 12
            if extraMonths == 0:
                values = data_one_series.reshape(completeYears, 12)
            else:
                extras = np.ones(12-extraMonths)    # so it doesnt affect the zero probability
                values = np.concatenate((data_one_series, extras)).reshape(completeYears+1, 12)

            hist_values = hist_data_one_series.reshape(calibration_years, 12)
            sigmas = computeSigmas(values, hist_values)
            # print('dists: ', dists)

            sigmas = sigmas.flatten()      # gets them into 1 dimension
            totalvals = sigmas.shape[0]

            # Get rid of the extra months
            sigmas = sigmas[:totalvals-(12-extraMonths)]

            # set all values > 3.1 to 3.1 and all values < -3.1 to -3.1
            sigmas[sigmas > 3.1] = 3.1
            sigmas[sigmas < -3.1] = -3.1

            spi[:,i,k] = sigmas

    return spi

def spi(
        fname="Interns2019/f09_g16.B.cobalt.FRAM.MAY.PRECT.200005-208106.nc",
        hfname="controlPRECT.nc",
        var='PRECT',
        histYears=30,
        newVar='SPInew',
        description='Standard Precipitation Index'
    ):
    """
    Inputs:
        fname: nc file path with precipitation data
        hfname: nc file path with historical precipitation data
        var: name of precipitation variable in the nc files
        histYears: int representing the amount of years of historical data to use to fit the gamma distribution
        newVar: name of new spi variable
        description: short description/full name of spi variable: ex: 'Standard Preipitation Index'
    Output:
        Calculates SPI into a 3d array
        Copies this data to the nc file
        returns the spi array
    """
    # reading netCDF dataset into variable and able to overwrite the data as well with 'r+'
    data = Dataset(fname, 'r+')
    hData = Dataset(hfname, 'r+')

    # get variable metadata
    variables = data.variables
    varData = variables[var]

    # get control data variables
    hVariables = hData.variables
    hVarData = hVariables[var]

    # get units, dims, and dtype
    # units = varData.units
    dims = varData.dimensions
    dtype = varData.dtype
    timeS, latS, lonS = varData.shape

    # get the time steps
    time = variables['time']
    steps = time.shape[0]

    initialize = False
    if newVar not in variables.keys():
        data.createVariable(newVar, dtype, dims)
        initialize = True

    spiv = variables[newVar]

    # if this is the first time creating it then 
    if initialize:
        spiv.long_name = description
    
    spiArr = calculate_over_full_series_with_distribution(varData, hVarData[0:histYears*12], starting_month=5, calibration_years=30)

    for t in range(timeS):
        spiv[t] = spiArr[t]

    return spiv

def yearlyAvg(monthly):
    """ Converts from a monthly data to a nparray of yearly data
        Does not adjust the years so if you start in June, then it will be a yearly average from Jun to the next May
    """
    time,lat,lon = monthly.shape
    
    # How many years after we get rid of the excess original months
    years = time // 12

    # How many extra months at the end we need to ignore
    correction = time % 12

    yearly = np.zeros((years, lat, lon))

    count = 0
    for i in range(lat):
        for j in range(lon):
            if count%1000==0:
                print(count)
            count += 1
            series = np.copy(monthly[:,i,j])[time-correction]

            series = series.reshape((years, 12))

            average = np.average(series, axis=1)

            yearly[:,i,j] = average

    return yearly

def seasonalAvg(monthly, startingMonth):
    """ Converts from a monthly data to a nparray of seasonal averaged data
        Seasons are Dec-Feb, March-May, Jun-August, Sep-Nov
    """
    time,lat,lon = monthly.shape

    # How many month you need to ignore at the begging to get the first season
    # Start month must me Dec, March, Jun, or Sep
    startIndex = (3-startingMonth)%3
    print('startIndex:', startIndex)

    # How many years after we get rid of the excess original months
    seasons = (time-startIndex) // 3
    print('seasons:', seasons)

    # How many extra months at the end we need to ignore
    correction = (time-startIndex) % 3
    print('correction:', correction)

    average = np.zeros((seasons, lat, lon))

    count = 0
    for i in range(lat):
        for j in range(lon):
            if count%1000==0:
                print(count)
            count += 1
            series = np.copy(monthly[:,i,j])[startIndex: time-correction]

            series = series.reshape((seasons, 3))

            seriesAverage = np.average(series, axis=1)

            average[:,i,j] = seriesAverage

    return average

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

            # Replaces all the inf and -inf values with nan values
            newSeries = np.where(series != math.inf, series, float('nan'))
            newSeries = np.where(newSeries != -math.inf, newSeries, float('nan'))

            # get the max and min arrays for each 
            maxSeries = np.nanmax(newSeries, axis=1)
            minSeries = np.nanmin(newSeries, axis=1)

            for y in range(years):
                newSeries[y] = np.where(newSeries[y] != math.inf, newSeries[y], maxSeries[y])
                newSeries[y] = np.where(newSeries[y] != -math.inf, newSeries[y], minSeries[y])


            average = np.nanmean(newSeries, axis=1)

            yearly[:,i,j] = average

    return yearly

def spiYears(
        fname="Interns2019/f09_g16.B.cobalt.FRAM.MAY.PRECT.200005-208106.nc",
        startingMonth=5
    ):
    """ Calculates yearly SPI from montly SPI data
    """
    # reading netCDF dataset into variable and able to overwrite the data as well with 'r+'
    data = Dataset(fname, 'r+')

    # Get the SPI variable
    variables = data.variables
    spi = variables['SPI']

    spiYearly = yearlyAvgJan(spi, startingMonth)

    return spiYearly  
