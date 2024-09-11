#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Description:
    Download aia.lev1_euv_12s and hmi.M_720s images, at a 1 hour cadence for 
    the specified Carrington rotations. Data are downloaded if, and only if, 
    all 8 files exist and have a 'QUALITY' key of '0' 
Created on Thu Jul  1 15:45:22 2021

@author: jgra
"""

# In[1]:
# Import Libraries and Tools
import pandas as pd
import glob
import numpy as np
import os
import datetime
import DataManagmentTools as dmt
import copy
import drms
import time as sleep

# Setup DRMS Client
client = drms.Client()

# In[2]:
# Key Variables

# Carrington Rotation Data
CarringtonFile = 'Carrington Rotation Start Dates.csv'
StartRotation  =  2099 # Download from StartRotation to EndRotation, inclusive
EndRotation    =  2101
RotCadence     =  1    # Download every nth rotation

# Dataset folders
dataFolder  = '/home/jgra/Coronal Holes/newDataset/' # This is the folder where you want the dataset to be placed.
traceFolder = 'DownloadLists/' # The "traceFolder" is assumed to be inside of the folder where this script is 
                               # located. It holds the .csv files that provide a summary of the data

# Request Data
AIAserise  = 'aia.lev1_euv_12s' # Type of images to request from AIA
HMIserise  = 'hmi.M_720s'       # Type of images to request from HMI
notify     = ''                 # visit http://jsoc.stanford.edu/ajax/register_email.html to register email
AIAsegment = 'image'            # We want the image data, not the "sharps"
HMIsegment = 'magnetogram'      # This is unnecessary, but included for consistency
reckeys    = 'T_REC, QUALITY'   # Data acquired in search
protocol   = 'fits'             # include all keyword data in file headers
ttlFiles   = 8                  # There are 7 EUVs and 1 Magnetogram per set/hour

# Pause between request
sleepTime = 30 # seconds

# Inform User
verbose = True

# Redundant Error Checking for Timeout Issues
cycle = 0
cycleLimit = 5

# In[3]:
# Read File
data = pd.read_csv(CarringtonFile,header=0)
keys = data.keys()

# Starting Process
startIndex = np.where(data[keys[0]]==StartRotation)[0][0]

# In[4]:
# Cycle Through Data
for i in range(startIndex,len(data)-1,RotCadence):
    
    # Find Carrington Rotation Name/Number
    RotNum = data[keys[0]][i]
    
    # Check if Within Specified Download Range (Exit Otherwise)
    if RotNum > EndRotation:
        break
    
    # Inform User
    if verbose:
        print('Rotation:',RotNum)
    
    # Carrington Start Datetime (rounded via floor operation to nearest hour)
    year   = data[keys[1]][i]
    month  = data[keys[2]][i]
    day    = int(data[keys[3]][i])
    hour   = int((data[keys[3]][i]-day) * 24)
    minute = 0
    second = 0
    startDate = datetime.datetime(year, month, day, hour, minute, second)
    
    # Carrington End Datetime (rounded via floor operation to nearest second)
    year   = data[keys[1]][i+1]
    month  = data[keys[2]][i+1]
    day    = int(data[keys[3]][i+1])
    hour   = int((data[keys[3]][i+1]-day) * 24)
    minute = int((((data[keys[3]][i+1]-day) * 24)-hour) * 60)
    second = int((((((data[keys[3]][i+1]-day) * 24)-hour) * 60)-minute)*60)
    endDate = datetime.datetime(year, month, day, hour, minute, second)
    
    # Ensure We Do Not Pass Data End
    if endDate < datetime.datetime.now():
        
        # Create Document Name
        traceName = 'CR' + str(RotNum) + '.csv'
        traceName = traceFolder + traceName
        
        # Create Document if Needed
        if not os.path.exists(traceName):
            with open(traceName,'w+') as f:
                f.write('T_intended,T_REC,94,131,171,193,211,304,335,Magnetogram\n')
                
        # Check Download Progress Otherwise
        else:
            RotData = pd.read_csv(traceName)
            rotKeys = RotData.keys()
            # Find Last Time
            if len(RotData) != 0:
                # Find where we Left Off
                LeftOff = RotData[rotKeys[0]][len(RotData)-1]
                startDate = dmt.timeUnformat(LeftOff)
                startDate += datetime.timedelta(hours=1)
                
        # Check for and, if Needed Create Folder for Download
        crFolder = dataFolder + 'CR' + str(RotNum) + '/'
        if not os.path.exists(crFolder):
            os.mkdir(crFolder)
        
        # Define Date for Search
        date = copy.deepcopy(startDate)
        
        # Cycle Through Dates with One Hour Cadence
        while date <= endDate:
            
            try: # Catch Errors from Missing Entries - and Thus Missing Keywords
                
                # Format Data Request
                T_intended = dmt.timeFormat(date) + 'Z'
                AIArequest = AIAserise + '[' + T_intended + ']'
                HMIrequest = HMIserise + '[' + T_intended + ']'
                
                # Inform User
                if verbose:
                    print('   ',T_intended)
                    print('        Checking File Quality')
                
                # Generate Query to get Record Time and QUALITY information
                AIAquery = client.query(AIArequest,seg=AIAsegment,key=reckeys)
                HMIquery = client.query(HMIrequest,seg=HMIsegment,key=reckeys)
                T_REC = AIAquery[0][reckeys.split(', ')[0]][0]
                QUALITY = AIAquery[0][reckeys.split(', ')[1]]
                QUALITY = np.hstack([QUALITY,HMIquery[0][reckeys.split(', ')[1]]])
                
                # Download if Data are of Expected Quality and Exist
                if len(QUALITY) == ttlFiles and np.max(QUALITY) == 0:
                
                    # Generate Save Folder for Individual Hour
                    outDir = T_REC.replace(':','')
                    if not os.path.exists(crFolder + outDir):
                        os.mkdir(crFolder + outDir)
                        
                    # Inform User
                    if verbose:
                        print('        Requesting AIA DATA')
                    
                    # Perform AIA Request
                    AIArequest = AIArequest + '{' + AIAsegment + '}'
                    result = client.export(AIArequest, method='url', 
                                           protocol=protocol, email=notify)
                        
                    # Download AIA Data
                    result.download( crFolder + outDir)
                    
                    # Power Nap, so as not to overload JSOC with requests
                    sleep.sleep(int(sleepTime/2))
                    
                    # Inform User
                    if verbose:
                        print('        Requesting HMI DATA')
                    
                    # Perform HMI Request
                    HMIrequest = HMIrequest + '{' + HMIsegment + '}'
                    result = client.export(HMIrequest, method='url', 
                                           protocol=protocol, email=notify)
                    
                    # Download HMI Data
                    result.download( crFolder + outDir)
                    
                    # Inform User
                    if verbose:
                        print('        Updating',os.path.basename(traceName))
                        
                    # Generate List of Files
                    files = glob.glob(crFolder + outDir + '/*.fits')
                    
                    # Sort List
                    missingFiles = ttlFiles * 1
                    filesOrdered = ['EMPTY','EMPTY','EMPTY','EMPTY',
                                    'EMPTY','EMPTY','EMPTY','EMPTY']
                    for file in files:
                        if 'Z.94.i' in os.path.basename(file):
                            filesOrdered[0] = outDir + '/' + os.path.basename(file)
                            missingFiles = missingFiles - 1
                        elif 'Z.131.i' in os.path.basename(file):
                            filesOrdered[1] = outDir + '/' + os.path.basename(file)
                            missingFiles = missingFiles - 1
                        elif 'Z.171.i' in os.path.basename(file):
                            filesOrdered[2] = outDir + '/' + os.path.basename(file)
                            missingFiles = missingFiles - 1
                        elif 'Z.193.i' in os.path.basename(file):
                            filesOrdered[3] = outDir + '/' + os.path.basename(file)
                            missingFiles = missingFiles - 1
                        elif 'Z.211.i' in os.path.basename(file):
                            filesOrdered[4] = outDir + '/' + os.path.basename(file)
                            missingFiles = missingFiles - 1
                        elif 'Z.304.i' in os.path.basename(file):
                            filesOrdered[5] = outDir + '/' + os.path.basename(file)
                            missingFiles = missingFiles - 1
                        elif 'Z.335.i' in os.path.basename(file):
                            filesOrdered[6] = outDir + '/' + os.path.basename(file)
                            missingFiles = missingFiles - 1
                        elif 'hmi.m' in os.path.basename(file):
                            filesOrdered[7] = outDir + '/' + os.path.basename(file)
                            missingFiles = missingFiles - 1
                    
                    # Create Row in Document
                    row = T_intended + ',' + T_REC + ','
                    for j in range(len(filesOrdered)):
                        row = row + filesOrdered[j]
                        if  j!=7:
                            row = row + ','
                        else:
                            row = row + '\n'
                    
                    # Update Document if all files present
                    if missingFiles == 0:
                        with open(traceName,'a+') as f:
                                        f.write(row)
                    
                    # Update Date for Next Request
                    date += datetime.timedelta(hours=1)
                    
                    # Inform User
                    if verbose:
                        print('    Sleeping',sleepTime,'Seconds...\n')
                    
                    # Wait, so as not to Overload JSOC with Requests
                    sleep.sleep(sleepTime)
                    
                # Quality Keyword Error
                else:
                    
                    # Inform User
                    if verbose:
                        print('        Quality Error, skipping\n')
                    
                    # Update Date for Next Request
                    date += datetime.timedelta(hours=1)
            
            except: # Catch errors from missing entries - and thus missing keywords
                
                # We will attempt to download the data multiple times to ensure
                # that the error was not not a timeout issue.
                if cycle > cycleLimit: 
                    
                    # Update counter to ensure that we check the appropriate 
                    # number of times during the next try/except cycle
                    cycle = 0 
                    
                    # Inform User
                    if verbose:
                        print('        Data Error, skipping')
                        
                    # Update Date for next request (move on to next hour)
                    date += datetime.timedelta(hours=1)
                    
                    # Inform user
                    if verbose:
                        print('    Sleeping',sleepTime,'Seconds...\n')
                    
                    # Wait, so as not to overload JSOC with requests
                    sleep.sleep(sleepTime)
                    
                else:
                    # Update try/except Counter
                    cycle += 1
                    
                    # Inform User
                    if verbose:
                        print('    Data Error, sleeping',sleepTime,'seconds...\n')
                    
                    # Wait, so as not to overload JSOC with requests before
                    # trying same request again
                    sleep.sleep(sleepTime)

# In[5]:
# End Process

print('**Process Complete**')
