#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Description:
    The Following Code Repeats the process from the script <FitsACWE_I.py>
    on level 1.5 dataproducts and across all image groups in the selected
    doucment.
    
Created on Fri Feb 18 14:39:27 2022
updated on Thu Mar  2 15:36:20 2023

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
import pngScale2

# Root directory of the project
ROOT_DIR = os.path.abspath("../")

# Import ACWE Tools
sys.path.append(ROOT_DIR)
from ACWE_python_spring_2023 import acweFunctions_v6, acweSaveSeg_v5

# In[2]:
# Key Varibles

# Dataset
# Dataset folders
dataFolder  = '/home/jgra/Coronal Holes/newDataset/' # Update to reflect Dataset Location
traceFolder = os.path.join(ROOT_DIR, 'DatasetTools/DownloadLists/')
CR          = 'CR2133' # Update to reflect chosen Carrington Rotation

# SaveFolder
saveFolder = '/mnt/coronal_holes/Code Paper I Observations/IntInv/'

# ACWE Prefix
acwePrefix = 'ACWE_IntInv_' # Prefix for ACWE
overwrite = False    # Run from begining 

# ACWE Parameters 
acweChoice = '193'
noScale = True
resize_param = 8          # These values are the defualt values taken from:
foreground_weight = 1     # L. E. Boucheron, M. Valluri, and R. T. J. McAteer, 
background_weight = 1/50. # "Segementation of Coronal Holes Using Active 
alpha = 0.3               # Contours Without Edges," Solar Physics, 
narrowband = 2            # vol. 291, pp. 2353-2372, 2016
N = 10

acweVerbose = False           # Show process mask generation 
correctLimbBrightening = True # Correct for Limb Brightening
rollingAlpha = 0              # Ensure alpha of 0.3 is utalized
fillInitHoles = True          # Fill holes in mask before running ACWE

# Inform user
verbose = True

# ACWE on 211 Data
# acweChoice = 211; alpha = 0.3; background_weight = 1/100. # Old - refine parameters

# Pause between failed request
sleepTime = 30

# Intensity Scaling Parameters
# scale types = [[min,max,scaletype,"restore"forAcWE,scaleName]]
scaleTypes = [["Imin","Imax","linear",False,"LinearCompressFull."],
              ["Imin","Imax","linear",True ,"LinearCompressFullRestored."],
              [0     ,"Imax","Linear",False,"LinearCompress0toMax."],
              [0     ,"Imax","linear",True ,"LinearCompress0toMaxRestored."],
              ["Smin","Smax","linear",False,"LinearCompressSolarLimits."],
              ["Smin","Smax","linear",True ,"LinearCompressSolarLimitsRestored."],
              [20.0  ,2500.0,"linear",False,"LinearCompressDefault."],
              [20.0  ,2500.0,"linear",True ,"LinearCompressDefaultRestored."],
              ["Imin","Imax","log10", False,"Log10CompressFull."],
              ["Imin","Imax","log10", True ,"Log10CompressFullRestored."],
              [0     ,"Imax","log10" ,False,"Log10Compress0toMax."],
              [0     ,"Imax","log10" ,True ,"Log10Compress0toMaxRestored."],
              ["Smin","Smax","log10" ,False,"Log10CompressSolarLimits."],
              ["Smin","Smax","log10" ,True ,"Log10CompressSolarLimitsRestored."],
              [20.0  ,2500.0,"log10" ,False,"Log10CompressDefault."],
              [20.0  ,2500.0,"log10" ,True ,"Log10CompressDefaultRestored."]]

# In[3]:
def make_circle_mask(c,im_dims,r):
    # defomes a binary image with image dimensions im_dims of a circle with 
    # center c and radius r
    cx = c[0]
    cy = c[1]
    ix = im_dims[0]
    iy = im_dims[1]
    x,y = np.meshgrid(np.arange(-(cx),(ix-cx),1),np.arange(-(cy),(iy-cy),1))
    c_mask = (x**2+y**2)<=r**2
    return c_mask

def loadfits(dataFolder,CR,file,loaded):
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
    return I,H,loaded

# In[4]:
# Open file and get list of images

# Open Carrington Rotation Document
CarringtonFile = traceFolder + CR + '.csv'
data = pd.read_csv(CarringtonFile,header=0)
keys = data.keys()

# Find Image group we want 
acweChoice = np.where(keys == acweChoice)[0][0]

# In[5]:
# Prepare for ACWE

# Create or find save location
crSaveFolder = saveFolder + CR + '/'
if not os.path.exists(crSaveFolder):
    os.makedirs(crSaveFolder)

# In[6]:
# Perform ACWE

# Cycle Through Dataset
for file in data[keys[acweChoice]]:#[len(data[keys[acweChoice]])-1:0:-1]:
    
    # Placement Folder
    acweFolder = file.split('/')[0] + '/'
    if not os.path.exists(crSaveFolder + acweFolder):
        os.mkdir(crSaveFolder + acweFolder)
        
    # Account for possible timeout when updating header
    loaded = False
    
    
    # Cycle Through Defined Scales
    for scaleType in scaleTypes:
        
        # ACWE Name
        acweFile = acwePrefix + scaleType[-1] + os.path.basename(file)
        acweFile = acweFolder + acweFile + '.npz'
        
        # Run ACWE if Needed
        if overwrite or not os.path.exists(crSaveFolder + acweFile):
            
            # Load Data
            if not loaded:
                
                if verbose:
                    print('Generating results for', os.path.basename(file))
                    print('    Opening EUV Image')
                    
                I,H,loaded = loadfits(dataFolder,CR,file,loaded)
                        
            # Inform User
            if verbose:
                print('        Performing ACWE On',scaleType[-1][:-1])
            
            # Update Scale As Needed
            Imin = scaleType[0];Imax = scaleType[1]
            if Imin == 'Imin':   # Image Minimum
                Imin = np.min(I)
            elif Imin == 'Smin': # Minimum Value on Solar Disk
                Ishape = np.asarray(I.shape)
                c_mask = make_circle_mask(Ishape/2.,Ishape,H['R_SUN'])
                Imin = np.min(I[c_mask])
            if Imax == 'Imax':   # Image Maximum
                Imax = np.max(I)
            elif Imax == 'Smax': # Maximum Value on Solar Disk
                Ishape = np.asarray(I.shape)
                c_mask = make_circle_mask(Ishape/2.,Ishape,H['R_SUN'])
                Imax = np.max(I[c_mask])
                
            # Scale (and Unscale)
            Itmp = pngScale2.scale(I,scaleType[2],Imin,Imax)
            if scaleType[3]:
                Itmp = pngScale2.unScale(Itmp,scaleType[2],Imin,Imax)
            
            # Save Scaling Inforamtion:
            image_preprocess = {'min':[scaleType[0],Imin],
                                'max':[scaleType[1],Imax],
                                'scaleType':scaleType[2],
                                'ReverseScale':scaleType[3]}
            
            if verbose:
                print('            Min:',image_preprocess['min'])
                print('            Max:',image_preprocess['max'])
                print('            ScaleType:',image_preprocess['scaleType'])
                print('            Reverse:',image_preprocess['ReverseScale'])
                print('            NaNs:',np.sum(np.isnan(Itmp)))
            
            try:
                # Run ACWE
                seg,m = acweFunctions_v6.run_acwe(Itmp,H,resize_param,
                                                  foreground_weight,
                                                  background_weight,
                                                  alpha,narrowband,
                                                  N,acweVerbose,
                                                  correctLimbBrightening,
                                                  rollingAlpha,
                                                  fillInitHoles)
                
                # Inform User
                if verbose:
                    print('            Saving Results')
                
                # Save Result
                init_mask_method = 'alpha*mean(qs)'
                alphar = alpha * 1
                acweSaveSeg_v5.saveSeg(crSaveFolder + acweFile,seg,H,
                                       correctLimbBrightening,resize_param,
                                       foreground_weight,background_weight,m,
                                       init_mask_method,fillInitHoles,alpha,
                                       alphar,narrowband,N,image_preprocess)
            except:
                
                # Inform User
                if verbose:
                    print('            Error Saving Empty')
                    
                # Save Result
                init_mask_method = 'alpha*mean(qs)'
                segShape = (np.asarray(I.shape)/resize_param).astype(int)
                seg = np.empty(segShape); seg[:] = np.nan
                alphar = alpha * 1
                m = np.empty(segShape); m[:] = np.nan
                acweSaveSeg_v5.saveSeg(crSaveFolder + acweFile,seg,H,
                                       correctLimbBrightening,resize_param,
                                       foreground_weight,background_weight,m,
                                       init_mask_method,fillInitHoles,alpha,
                                       alphar,narrowband,N,image_preprocess)
                
    # Inform User When Aplicible
    if verbose and loaded:
        print('\n\n')
        
# In[6]:
# End Process

print('**Process Complete**')
