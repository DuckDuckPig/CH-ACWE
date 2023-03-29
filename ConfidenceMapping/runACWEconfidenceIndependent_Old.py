#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Description:
    Run ACWE on specified folder, generating confidence maps and placing the
    results in a mirror folder structure. This version generates the standard
    512x512 Confidence Map outputs using the old, no-optimized code, the 
    traditional thresholding method and a rolling alpha.

Created on Fri Feb 11 20:13:07 2022

@author: jgra
"""

# In[1]:
# Import Libraries and Tools
import os 
import sys
import pandas as pd
import numpy as np
import time as sleep
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

# Dataset
# Dataset folders
dataFolder  = '/home/jgra/Coronal Holes/newDataset/'
traceFolder = os.path.join(ROOT_DIR, 'DatasetTools/DownloadLists/')
CR          = 'CR2133'

# SaveFolder
saveFolder = '/mnt/coronal_holes/Code Paper I Observations/ConMapStandardOld/'

# ACWE Prefix
acwePrefix = 'ACWEconMapOld.' # Prefix for ACWE
overwrite = False    # If True run from beginning 

# ACWE Parameters 
acweChoice = '193'
noScale = True
resize_param = 8          # These values are the default values taken from:
foreground_weight = 1     # L. E. Boucheron, M. Valluri, and R. T. J. McAteer, 
background_weight = 1/50. # "Segmentation of Coronal Holes Using Active 
alpha = 0.3               # Contours Without Edges," Solar Physics, 
narrowband = 2            # vol. 291, pp. 2353-2372, 2016
N = 10

# Range of background weights used in Confidence Map
background_weight = 1/np.array(range(10,101))

acweVerbose = False           # Show process mask generation 
correctLimbBrightening = True # Correct for Limb Brightening
rollingAlpha = 0.01           # Incrementally Increase Alpha when Failed Threshold
fillInitHoles = True          # Fill holes in mask before running ACWE

# Inform user
verbose = True

# Pause between failed request
sleepTime = 30

# # Time ACWE
# timeFile = os.path.join(ROOT_DIR,'ConfidenceMapping/') + CR + '_timeOld.csv'

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
        loaded = False
        while not loaded:
            
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
                
                # Exit Loop Last
                loaded = True
                
            except:
                # Inform User
                if verbose:
                    print('    Error Sleeping',sleepTime,'Seconds...')
                
                # Wait Then Try Again
                sleep.sleep(sleepTime)
                
                # Inform User
                if verbose:
                    print('    Opening EUV Image')
    
        # Inform user
        if verbose:
            print('    Running ACWE')
            
        # # Time
        # start = time.time()
        
        # Prepare for Loop
        outputShape = np.asarray(I.shape)/resize_param
        outputShape = np.hstack([len(background_weight),outputShape]).astype(int)
        SegsOld = np.empty(outputShape); SegsOld[:] = np.nan
        
        for i in range(len(background_weight)):
            
            seg,alphar,m = acweFunctions_v6.run_acwe(I,H,resize_param,
                                                     foreground_weight,
                                                     background_weight[i],
                                                     alpha,narrowband,
                                                     N,acweVerbose,
                                                     correctLimbBrightening,
                                                     rollingAlpha,
                                                     fillInitHoles)
            
            SegsOld[i] = seg.astype(int) * 1
        
        # Time
        # end = time.time()
        # timeTotal = end-start
        # timeTotal = str(timeTotal)
        
        # Inform User
        if verbose:
            print('    Saving Results\n\n')
        
        # Save Result
        init_mask_method = 'alpha*mean(qs)'
        
        acweSaveSeg_v5.saveSeg(crSaveFolder + acweFile,SegsOld,H,
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
