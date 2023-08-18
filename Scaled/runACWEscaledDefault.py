#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Description:
    Run ACWE on a specified folder, and place the results in a mirror folder
    structure. In particular, this version generates ACWE outputs, with the 
    traditional thresholding methods, and a rolling alpha.

Created on Fri Jul  9 12:57:16 2021
Updated on Tue Feb 21 12:14:04 2023

@author: jgra
"""

# In[1]:
# Import Libraries and Tools
import os 
import sys
import pandas as pd
import numpy as np
from astropy.io import fits
import sunpy.map
from aiapy.calibrate import register, update_pointing
import warnings
warnings.filterwarnings("ignore")

# Root directory of the project
ROOT_DIR = os.path.abspath("../")

# Import ACWE Tools
sys.path.append(ROOT_DIR)
from ACWE_python_spring_2023 import acweFunctions_v6, acweSaveSeg_v5

# import time

# In[2]:
# Key Variables

# Choose Resize Parameter
resize_param = 4 # 1, 2, and 4 were tested

# Dataset
# Dataset folders
dataFolder  = '/home/jgra/Coronal Holes/newDataset/' # Update to reflect Dataset Location
traceFolder = os.path.join(ROOT_DIR, 'DatasetTools/DownloadLists/')
CR          = 'CR2100' # Update to reflect chosen Carrington Rotation

# SaveFolder - Update to reflect location where data will be saved 
saveFolder = '/mnt/coronal_holes/Code Paper I Observations/Scaled/'

# Ensure Consistent ACWE Prefix
if resize_param == 8:
    acwePrefix = 'ACWE.'
else:
    acwePrefix = 'ACWEresize_param' + str(resize_param) + '.'

# Run Parameters
overwrite = False # If True run from beginning 

# ACWE Parameters 
acweChoice = '193'
noScale = True            # These values are the default values taken from:
foreground_weight = 1     # L. E. Boucheron, M. Valluri, and R. T. J. McAteer, 
background_weight = 1/50. # "Segmentation of Coronal Holes Using Active 
alpha = 0.3               # Contours Without Edges," Solar Physics, 
narrowband = 2            # vol. 291, pp. 2353-2372, 2016
N = 10

acweVerbose = False           # Show process mask generation 
correctLimbBrightening = True # Correct for Limb Brightening
rollingAlpha = 0.01           # Incrementally Increase Alpha when Failed Threshold
fillInitHoles = True          # Fill holes in mask before running ACWE

# Inform user
verbose = True

# ACWE on 211 Data
# acweChoice = 211; alpha = 0.3; background_weight = 1/100. # Old - refine parameters

# # Time ACWE
# timeFile = os.path.join(ROOT_DIR,'Scaled/') + CR
# if resize_param == 8:
#     timeFile = timeFile + '_timeStandard.csv'
# else:
#     timeFile = timeFile + '_time' + acwePrefix + 'csv'

# In[3]:
# Open file and get list of images

# Open Carrington Rotation Document
CarringtonFile = traceFolder + CR + '.csv'
data = pd.read_csv(CarringtonFile,header=0)
keys = data.keys()

# Find Image group we want 
acweChoice = np.where(keys == acweChoice)[0][0]

# In[4]:
# Prepare for ACWE

# Create or find save location
crSaveFolder = saveFolder + CR + '/'
if not os.path.exists(crSaveFolder):
    os.makedirs(crSaveFolder)
    
# # Prepare time file
# if not os.path.exists(timeFile):
#     with open(timeFile,'w+') as f:
#         f.write('file,time\n')
    
# In[5]:
# Perform ACWE

# Cycle Through Dataset
for file in data[keys[acweChoice]]:#[len(data[keys[acweChoice]])-1:0:-1]:
    
    # Placement Folder
    acweFolder = file.split('/')[0] + '/'
    if not os.path.exists(crSaveFolder + acweFolder):
        os.mkdir(crSaveFolder + acweFolder)
    
    # ACWE Name
    acweFile = acwePrefix + os.path.basename(file)
    acweFile = acweFolder + acweFile + '.npz'
    
    # Check for File
    if overwrite or not os.path.exists(crSaveFolder + acweFile):
        
        # Inform User
        if verbose:
            print('Generating', os.path.basename(acweFile))
            print('    Opening EUV Image')
        
        # Account for possible timeout when updating header
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
        
        # Inform user
        if verbose:
            print('    Running ACWE')
            
        # # Time
        # start = time.time()
        
        # Run ACWE
        seg,alphar,m = acweFunctions_v6.run_acwe(I,H,resize_param,
                                                 foreground_weight,
                                                 background_weight,
                                                 alpha,narrowband,
                                                 N,acweVerbose,
                                                 correctLimbBrightening,
                                                 rollingAlpha,
                                                 fillInitHoles)
        
        # # Time
        # end = time.time()
        # timeTotal = end-start
        # timeTotal = str(timeTotal)
        
        # Inform User
        if verbose:
            print('    Saving Results\n\n')
        
        # Save Result
        init_mask_method = 'alpha*mean(qs)'
        acweSaveSeg_v5.saveSeg(crSaveFolder + acweFile,seg,H,
                               correctLimbBrightening,resize_param,
                               foreground_weight,background_weight,m,
                               init_mask_method,fillInitHoles,alpha,
                               alphar,narrowband,N)
        
        # # Time
        # row = acweFile + ',' + timeTotal + '\n'
        # with open(timeFile,'a+') as f:
        #     f.write(row)
        
# In[6]:
# End Process

print('**Process Complete**')
