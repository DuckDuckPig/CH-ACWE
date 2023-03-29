#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Description:
    Collection of functions and tools for managing a dataset of .fits AIA 
    images.
Created on Tue Jun  2 12:12:43 2020

@author: jgra
"""

# In[1]:
# Import Libraries and Tools
import datetime
import os
import csv
import numpy as np
from astropy.io import fits

# In[2]:
def timeFormat(date):
    '''
    Convert Datatime object into the format that 
    Jsoc expects, ex: "2010-12-06T10:00:01"
    '''
    
    # Get elements
    year   = date.strftime('%Y')
    month  = date.strftime('%m')
    day    = date.strftime('%d')
    hour   = date.strftime('%H')
    minute = date.strftime('%M')
    sec    = date.strftime('%S')
    
    # Return formatted Date
    return year + '-' + month + '-' + day + 'T' + hour + ':' + minute + ':' + sec

def timeUnformat(date):
    '''
    Convert Jsoc time string (ex: "2010-12-06T10:00:01") 
    to Datetime object
    '''
    
    # get elements
    year   = int(float(date[ 0:4 ]))
    month  = int(float(date[ 5:7 ]))
    day    = int(float(date[ 8:10]))
    hour   = int(float(date[11:13]))
    minute = int(float(date[14:16]))
    sec    = int(float(date[17:19]))
    
    # Return dattime object
    return datetime.datetime(year,month,day,hour,minute,sec)

def timeFromFilenameMag(fname):
    '''
    Convert file name (ex: "hmi.m_720s.20100713_090000_TAI.1.magnetogram.fits") 
    to Datetime object
    '''
    
    date = os.path.basename(fname).split('.')[2]
    
    # get elements
    year   = int(float(date[ 0:4 ]))
    month  = int(float(date[ 4:6 ]))
    day    = int(float(date[ 6:8 ]))
    hour   = int(float(date[ 9:11]))
    minute = int(float(date[11:13]))
    sec    = int(float(date[13:15]))
    
    # Return dattime object
    return datetime.datetime(year,month,day,hour,minute,sec)

def timeFromFilenameAIA(fname):
    '''
    Convert file name (ex: "aia.lev1_euv_12s.2010-07-13T090008Z.193.image_lev1.fits") 
    to Datetime object
    '''
    
    date = os.path.basename(fname).split('.')[2]
    
    # get elements
    year   = int(float(date[ 0:4 ]))
    month  = int(float(date[ 5:7 ]))
    day    = int(float(date[ 8:10]))
    hour   = int(float(date[11:13]))
    minute = int(float(date[13:15]))
    sec    = int(float(date[15:17]))
    
    # Return dattime object
    return datetime.datetime(year,month,day,hour,minute,sec)


# In[3]:
def timeSearch(startDate,endDate,datasetFile,keys,
                         delimiter = ',',verbose=False):
    '''
    Given a dataset file that holds all the paths to all files returns all 
    observation times (inclusive) in which the desired images are available. 
    
    Inputs:
        startDate : datetime
            date and time used for start of search
        endDate : datetime
            date and time used for end of search
        datasetFile : str
            Path to dataset file. File must be .csv where the first row is the 
            keys and all other rows are the files. The first entry in each row
            must be the record time in the format that Jsoc expects, 
            ex: "2010-12-06T10:00:01"
        keys : [str]
            names of the images that an entry must have to be counted. Names
            must match how they appear in the top row of the dataset file
    Outputs:
        [str]
            Rows of the dataset file that meet the desired criteria
        [int]
            Position(s) within each entry that corresponds to the identified 
            keys
    '''
    
    # Inform User
    if verbose:
        print('Opening', os.path.basename(datasetFile))
    
    # Open Dataset File
    with open(datasetFile, 'r+') as f:
        csvData = csv.reader(f,delimiter=delimiter)
        fullDataset = []
        for datum in csvData:
            fullDataset.append(datum)
            
    # Separate dataset from keys
    headerRow    = fullDataset[0]
    fullDataset  = fullDataset[1:]
    
    # Determine key positions
    keyPositions = []
    for key in keys:
        for i in range(len(headerRow)):
            if key == headerRow[i]:
                keyPositions.append(i*1)
                
    # Inform User
    if verbose:
        print('Searching for valid entries')
    
    # Convert observation times to datetime
    times = []
    for datum in fullDataset:
        times.append(timeUnformat(datum[0]))
    times = np.array(times)
    
    # Find all entries in desired time range
    timeIndex = np.hstack(np.where(np.logical_and(times>=startDate, times <=endDate))) # I don't know why I need to hstack
    
    # Within Index find all entries that have the requested images
    keepIndex = []
    for index in timeIndex:
        valid = False
        try: # Check for Each Key
            v = []
            for key in keyPositions:
                v.append(fullDataset[index][key] != '')
            valid = all(v)
        except:
            pass
        if valid:
            keepIndex.append(index)
    
    # Reduce dataset to valid value
    keepDataset = []
    for index in keepIndex:
        keepDataset.append(fullDataset[index])
    return keepDataset,keyPositions

# In[4]:
def quickCheck(filePath):
    
    # Assume that image is corrupted
    valid = False
    
    # Check Image
    try:
        hdulist = fits.open(filePath) # Open File
        hdulist.verify('silentfix')   # Check if Standard Fix Works
        h = hdulist[1].header         # Check for Header
        I = hdulist[1].data           # Check for Data
        I = I*2                       # Check if data can be manipulated
        
        if h['QUALITY'] == 0 or h['QUALITY'] == '0':
            valid = True
    except:
        pass
    
    return valid
