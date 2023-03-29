#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Description:
    This script uses a preexisting file (i.e. the ones in the 'DownloadLists' 
    folder) to attempt to redownload/rebuild the specified dataset in the 
    specified folderpath.

Created on Wed Mar 29 10:08:24 2023

@author: jgra
"""

# In[1]
# Import Libraries and Tools
import pandas as pd
import glob
import os
import drms
import time as sleep

# Setup DRMS Client
client = drms.Client()

# In[2]
# Key Variables

# Dataset folders
dataFolder  = '/home/jgra/Coronal Holes/newDataset/' # Specify folder path here.

# Dataset File
traceFolder = 'DownloadLists/'   # Change to match 'CR*.csv' location.
CR          = 'CR2099'           # Change to match 'CR*.csv' name.

# Data Request
notify     = ''                 # visit http://jsoc.stanford.edu/ajax/register_email.html to register email
protocol   = 'fits'             # include all keyword data in file headers

# Inform User
verbose = True

# Pause between request
sleepTime = 15 # seconds

# In[3]
# Read Dataset File

traceFile = glob.glob(traceFolder + CR + '*')[0]

data = pd.read_csv(traceFile,header=0)
keys = data.keys()

# In[4]
# Find Full File Path

# Format to match folder setup
crFolder = os.path.basename(traceFile).split('.')[0] # CR21001.csv contains
crFolder = dataFolder + crFolder + '/'               # data on folder CR2101

# In[5]
# Check For Each File

# Cycle Through Folders
for i in range(len(data[keys[1]])):
    
    # Cycle Through Files
    for key in keys[2:]:
        
        # Get Filename
        filename = data[key][i]
        
        # Generate Full File Path
        filepath   = crFolder + filename
        
        # Download if missing
        if not os.path.exists(filepath):
            
            # Inform User
            if verbose:
                print(os.path.basename(filename), 'Not Found')
                print('    Generating Request')
            
            # Generte Request Data
            outDir  = filename.split('/')[0] + '/'
            serise  = os.path.basename(filename).split('.')[0]
            serise  = serise + '.' + os.path.basename(filename).split('.')[1]
            segment = os.path.basename(filename).split('.')[-2].split('_')[0]
            
            if 'aia' in serise:
                
                # Time
                time = data[keys[1]][i]
                
                # Wavelength
                wavelength = os.path.basename(filename).split('.')[-3]
                
                # Build request
                request = serise + '[' + time + '][' + wavelength + ']'
            
            elif 'hmi' in serise:
                
                # Time
                time =  data[keys[0]][i]
                
                # Build request
                request = serise + '[' + time + ']'
            
            # Full request Name
            request = request + '{' + segment + '}'
            
            # Verify folder exists:
            if not os.path.exists(crFolder + outDir):
                        os.makedirs(crFolder + outDir)
            
            # Inform User
            if verbose:
                print('    Making Request')
            
            # Perform Request
            result = client.export(request, method='url', 
                                   protocol=protocol, 
                                   email=notify)
            # Download Data
            result.download( crFolder + outDir)
            
            # Inform User
            if verbose:
                print('    Sleeping',sleepTime,'Seconds...\n')
                    
            # Wait, so as not to Overload JSOC with Requests
            sleep.sleep(sleepTime)
            
            
# In[6]
# End Process

print('**Process Complete**')
