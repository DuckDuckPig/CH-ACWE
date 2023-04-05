#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Description: Calculate and report the change in area and change in mean 
             intensity for the foreground region.

Created on Tue Nov 29 11:42:43 2022

@author: jgra
"""

# In[1]
# Import Libraires and tools
import os
import sys
import pandas as pd
import numpy as np
import glob
from astropy.io import fits
import sunpy.map
from aiapy.calibrate import register, update_pointing
import warnings
warnings.filterwarnings("ignore")

# ACWE utilities
# Root directory of the project
ROOT_DIR = os.path.abspath("../../")

# Import ACWE Tools
sys.path.append(ROOT_DIR)
from ACWE_python_spring_2023 import acweFunctions_v6 as af6, acweSaveSeg_v5
from ACWE_python_spring_2023.ACWE_python_v3 import correct_limb_brightening

# In[2]
# Key Varibles

# Dataset Selection
traceFolder = os.path.join(ROOT_DIR, 'DatasetTools/DownloadLists/')
CR = 'CR2133'

# Data
dataset = '/home/jgra/Coronal Holes/newDataset/'
conMap  = '/mnt/coronal_holes/Code Paper I Observations/ConMapStandardDefault/'

# ACWE Parameters 
acweChoice = '193'

# Inform User
verbose = True

# In[3]
# Open file and get list of images

# Open Carrington Rotation Document
CarringtonFile = traceFolder + CR + '.csv'
data = pd.read_csv(CarringtonFile,header=0)
keys = data.keys()

# Find Image group we want 
acweChoice = np.where(keys == acweChoice)[0][0]

# In[4]:
# Prepare for analysis

# Define Folder locations
ImageFolder  = dataset + CR + '/'
conMapFolder = conMap  + CR + '/'

# Determin dimensions of Confidence Map
file = data[keys[acweChoice]][0]
conMap = glob.glob(conMapFolder + '*/*' + os.path.basename(file) + '*')[0]
_,_,SEG = acweSaveSeg_v5.openSeg(conMap)

# Create Placeholder
segArea      = np.empty([len(data),len(SEG)]); segArea[:]      = np.nan
segIntnMean  = np.empty([len(data),len(SEG)]); segIntnMean[:]  = np.nan
segIntnStd   = np.empty([len(data),len(SEG)]); segIntnStd[:]   = np.nan
segIntnMin   = np.empty([len(data),len(SEG)]); segIntnMin[:]   = np.nan
segIntnMax   = np.empty([len(data),len(SEG)]); segIntnMax[:]   = np.nan
bkgWeights   = np.empty([len(data),len(SEG)]); bkgWeights[:]   = np.nan
seedArea     = np.empty(len(data));            seedArea[:]     = np.nan
seedIntnMean = np.empty(len(data));            seedIntnMean[:] = np.nan
seedIntnStd  = np.empty(len(data));            seedIntnStd[:]  = np.nan
seedIntnMin  = np.empty(len(data));            seedIntnMin[:]  = np.nan
seedIntnMax  = np.empty(len(data));            seedIntnMax[:]  = np.nan
IOO          = np.empty([len(data),len(SEG)]); IOO[:]          = np.nan

# In[5]:
# Perform analysis

# Cycle Through Dataset
for i in range(len(data[keys[acweChoice]])):
    
    # Get file
    file = data[keys[acweChoice]][i]
    
    # Inform User
    if verbose:
        print()
        print(os.path.basename(file))
        print('    Assesing Confidence Map')
        
    # Find Confidence Map
    conMap = glob.glob(conMapFolder + '*/*' + os.path.basename(file) + '*')[0]
    
    # Open Confidence map
    H,AH,SEG = acweSaveSeg_v5.openSeg(conMap)
    MASK   = AH['INIT_MASK']
    
    # Area of Init Mask
    seedArea[i] = np.sum(MASK.astype(int))
    
    # List of BackgroundWeights
    bkgWeights[i] = AH['BACKGROUND_WEIGHT']
    
    # Area of each segmentation and intersection with original
    for j in range(len(SEG)):
        
        # Segmentation Area
        segArea[i,j] = np.sum(SEG[j].astype(int))
        
        # Percet of orignal mask retained
        num = SEG[j].astype(bool) & MASK.astype(bool)
        num = np.sum(num.astype(int))
        den = np.sum(MASK.astype(int))
        IOO[i,j] = (num.astype(float))/(den.astype(float))
    
    # Inform User
    if verbose:
        print('    Opening original Image for Evaluation')
    
    # Open Original Image
    success = False
    while not success:
        try:
            # Extract Image and Header Data
            hdulist = fits.open(ImageFolder+file)
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
                I = J*1
                H = dict(h)
            success = True
        except:
            pass
    
    if verbose:
        print('    Evaluating Original Image')
    
    # Resize image
    I,im_size,sun_radius,sun_center = af6.resize_EUV(I,H,AH['RESIZE_PARAM'])
    
    # Correct for limb Brightening
    if AH['CORRECT_LIMB_BRIGHTENING']:
        I = correct_limb_brightening.correct_limb_brightening(I,sun_center,
                                                              sun_radius)
    
    # Seed Intensity Stats
    seedIntnMean[i] = np.mean(I[MASK.astype(bool)])
    seedIntnStd[i]  = np.std(I[MASK.astype(bool)])
    seedIntnMin[i]  = np.min(I[MASK.astype(bool)])
    seedIntnMax[i]  = np.max(I[MASK.astype(bool)])
    
    # Intensity each segmentation
    for j in range(len(SEG)):
        segIntnMean[i,j] = np.mean(I[SEG[j].astype(bool)])
        segIntnStd[i,j]  = np.std(I[SEG[j].astype(bool)])
        segIntnMin[i,j]  = np.min(I[SEG[j].astype(bool)])
        segIntnMax[i,j]  = np.max(I[SEG[j].astype(bool)])

# In[6]:
# Save Results for Graphing

if verbose:
    print()
    print('Saving Results')
    print()
resultsFolder = 'ConfidenceMapping/AnalysisGrowthAndIntensity/Results/'
resultsFolder = os.path.join(ROOT_DIR,resultsFolder)
if not os.path.exists(resultsFolder):
    os.mkdir(resultsFolder)
filename = CR + '.GrowthAndIntensity.npz'
filename = resultsFolder + filename
np.savez(filename,seedArea,seedIntnMean,seedIntnStd,seedIntnMin,seedIntnMax,
         segArea,segIntnMean,segIntnStd,segIntnMin,segIntnMax,bkgWeights,IOO)

# In[7]:
# End Program
print('**Process Complete**')