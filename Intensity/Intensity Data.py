#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  9 15:07:08 2023

@author: jgra
"""

# In[1]
# Import Libraries and Tools
import os
import pandas as pd
import numpy as np
from astropy.io import fits
import sunpy.map
from aiapy.calibrate import register, update_pointing
import warnings
warnings.filterwarnings("ignore")

# Root directory of the project
ROOT_DIR = os.path.abspath("../")

# In[2]
# Key Varibles

# Dataset
# Dataset folders
dataFolder  = '/home/jgra/Coronal Holes/newDataset/' # Update to reflect Dataset Location
traceFolder = os.path.join(ROOT_DIR, 'DatasetTools/DownloadLists/')
CR          = 'CR2133' # Update to reflect chosen Carrington Rotation

# ACWE Parameters 
acweChoice = '193'

# Inform user
verbose = True

# In[3]
# Open file and get list of images

# Open Carrington Rotation Document
CarringtonFile = traceFolder + CR + '.csv'
data = pd.read_csv(CarringtonFile,header=0)
keys = data.keys()

# Find Image group we want 
acweChoice = np.where(keys == acweChoice)[0][0]

# In[4]
# Prepare to determin intensity range

MIN = float(np.finfo('d').max)
MAX = float(np.finfo('d').min)

# In[5]
# find min and max of CR

# Inform User
if verbose:
    print(CR)

# Cycle Through Dataset
for file in data[keys[acweChoice]]:
    
    # Inform User
    if verbose:
        print('   ', os.path.basename(file))
        print('        Opening EUV Image')
    
    # Attempt to open and update file
    success = False
    while not success:
        try:
            # Extract Image and Header Data
            hdulist = fits.open(dataFolder+str(CR) +'/'+file)
            hdulist.verify('silentfix') # no clue why this is needed for successful data read
            h = hdulist[1].header
            J = hdulist[1].data
            hdulist.close()
            
            # Update to Level 1.5 Data Product
            if h['LVL_NUM'] < 1.5:
                m = sunpy.map.Map((J,h))    # Create Sunpy Map
                m = update_pointing(m)      # Update Header based on Latest Information
                m_registrered = register(m) # Recenter and rotate to Solar North
                I = m_registrered.data
                # Undo Keword Renaming
                H = dict()
                for k in m_registrered.meta.keys(): 
                    H[k.upper()] = m_registrered.meta[k] 
            # Skip if already Level 1.5
            else:
                # Convert header to dictionary
                m = sunpy.map.Map((J,h)) # Create Map
                H = dict()
                for k in m.meta.keys():
                    h[k.upper()] = m.meta[k]
                I = J*1 # Copy image
            success = True
        except:
            pass
        
    # Inform User
    if verbose:
        print('        min:',np.min(I))
        print('        max:',np.max(I))
    
    # Get Min and Max
    if np.min(I)<MIN:
        MIN = np.min(I)
    if np.max(I)>MAX:
        MAX = np.max(I)
        
# In[6]
# End Process

print()
print(CR)
print('min:',MIN)
print('max:',MAX)
print()

print('**Process Complete**')