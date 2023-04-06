#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Description:
    Comparison of Segmentations as a function of lambda_i/Lambda_0 value
    for old and new confidence map methods and comparison between old and new
    confidence map methods.
    
Created on Mon Feb 14 11:00:29 2022

@author: jgra
"""
# In[1]:
# Import Libraries and Tools
import os
import sys
import pandas as pd
import numpy as np
import glob
from skimage.metrics import mean_squared_error as mse
from skimage.metrics import normalized_root_mse as nrmse
from skimage.metrics import structural_similarity as ssim


# ACWE utilities
# Root directory of the project
ROOT_DIR = os.path.abspath("../../")

# Import ACWE Tools
sys.path.append(ROOT_DIR)
from ACWE_python_spring_2023 import acweSaveSeg_v5
from Metrics import consistancyErrorMetricsIII as cem, JaccardIndexMetric as jim

# In[2]:
# Key Variables

# Dataset
traceFolder = os.path.join(ROOT_DIR, 'DatasetTools/DownloadLists/')
CR = 'CR2133'

# ACWE Folders
saveFolderOld  = '/mnt/coronal_holes/Code Paper I Observations/ConMapStandardOld/'
saveFolderNew = '/mnt/coronal_holes/Code Paper I Observations/ConMapStandardDefault/'

# # ACWE Compare Prefix
# oldACWEprefix = 'ACWEconMapOld.'
# mewACWEprefix = 'ACWEconMap.'

# ACWE Parameters 
acweChoice = '193'

# Inform user
verbose = True

# In[3]:
# Open file and get list of images

# Open Carrington Rotation Document
CarringtonFile = traceFolder + CR + '.csv'
data = pd.read_csv(CarringtonFile,header=0)
keys = data.keys()

# Find Image group we want 
acweChoice = np.where(keys == acweChoice)[0][0]

# In[4]:
# Prepare for analysis

# Create or find save location
crSaveFolder1 = saveFolderOld  + CR + '/'
crSaveFolder2 = saveFolderNew + CR + '/'

# Open Sample
acweFolder = data[keys[acweChoice]][0].split('/')[0] + '/'
oldConMapFolder = crSaveFolder1 + acweFolder
oldConMap = glob.glob(oldConMapFolder + '*' + os.path.basename(data[keys[acweChoice]][0]) + '*')[0]
oH,oAH,oSEG = acweSaveSeg_v5.openSeg(oldConMap)

#  Create Placeholder for stats
oldSegArea       = np.empty([len(data),len(oSEG)]); oldSegArea[:]       = np.nan
oldForwardUnique = np.empty([len(data),len(oSEG)]); oldForwardUnique[:] = np.nan
oldBackUnique    = np.empty([len(data),len(oSEG)]); oldBackUnique[:]    = np.nan
newSegArea       = np.empty([len(data),len(oSEG)]); newSegArea[:]       = np.nan
newForwardUnique = np.empty([len(data),len(oSEG)]); newForwardUnique[:] = np.nan
newBackUnique    = np.empty([len(data),len(oSEG)]); newBackUnique[:]    = np.nan
MSE   = np.empty(len(data)); MSE[:]   = np.nan
NRMSE = np.empty(len(data)); NRMSE[:] = np.nan
SSIM  = np.empty(len(data)); SSIM[:]  = np.nan
GCE   = np.empty(len(data)); GCE[:]   = np.nan
LCE   = np.empty(len(data)); LCE[:]   = np.nan
IOUw  = np.empty(len(data)); IOUw[:]  = np.nan 

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
        
    # Find Confidence Maps
    acweFolder = file.split('/')[0] + '/'
    oldConMapFolder = crSaveFolder1 + acweFolder
    oldConMap = glob.glob(oldConMapFolder + '*' + os.path.basename(file) + '*')[0]
    newConMapFolder = crSaveFolder2 + acweFolder
    newConMap = glob.glob(newConMapFolder + '*' + os.path.basename(file) + '*')[0]
    
    # Open Files
    oH,oAH,oSEG = acweSaveSeg_v5.openSeg(oldConMap)
    nH,nAH,nSEG = acweSaveSeg_v5.openSeg(newConMap)
    
    # Inform User
    if verbose:
        print('    Finding Overlaps')
        
    # Calculate Overlap
    for j in range(len(oSEG)):
        oldSegArea[i,j] = np.sum(oSEG[j].astype(int))
        if j > 0:
            jof = np.logical_and(oSEG[j],np.logical_not(oSEG[j-1])) * 1
            oldForwardUnique[i,j] = np.sum(jof)
        if j < (len(oSEG) - 1):
            job = np.logical_and(oSEG[j],np.logical_not(oSEG[j+1])) * 1
            oldBackUnique[i,j] = np.sum(job)
    for j in range(len(nSEG)):
        newSegArea[i,j] = np.sum(nSEG[j].astype(int))
        if j > 0:
            jnf = np.logical_and(nSEG[j],np.logical_not(nSEG[j-1])) * 1
            newForwardUnique[i,j] = np.sum(jnf)
        if j < (len(nSEG) - 1):
            jnb = np.logical_and(nSEG[j],np.logical_not(nSEG[j+1])) * 1
            newBackUnique[i,j] = np.sum(jnb)
    
    # Inform User
    if verbose:
        print('    Comparing Maps')
    
    # Compare Maps
    oseg = np.sum(oSEG.astype(int),axis=0)
    nseg = np.sum(nSEG.astype(int),axis=0)
    MSE[i]    = mse(oseg,nseg)
    NRMSE[i]  = nrmse(oseg, nseg)
    SSIM[i],_ = ssim(oseg,nseg,full=True)
    GCE[i],LCE[i] = cem.CE(oseg, nseg)
    IOUw[i] = jim.IOU(oseg,nseg,binary=False)
    
# In[6]:
# Save Results for Graphing
resultsFolder = 'ConfidenceMapping/AnalysisInternal/Results/'
resultsFolder = os.path.join(ROOT_DIR,resultsFolder)
if not os.path.exists(resultsFolder):
    os.mkdir(resultsFolder)
filename = os.path.basename(CarringtonFile) + '.conCompare.npz'
filename = resultsFolder + filename
np.savez_compressed(filename,MSE,NRMSE,SSIM,GCE,LCE,IOUw)
filename = os.path.basename(CarringtonFile) + '.oldConOverlap.npz'
filename = resultsFolder + filename
np.savez_compressed(filename,oldSegArea,oldForwardUnique,oldBackUnique)
filename = os.path.basename(CarringtonFile) + '.newConOverlap.npz'
filename = resultsFolder + filename
np.savez_compressed(filename,newSegArea,newForwardUnique,newBackUnique)

# In[7]:
# End Program
print('**Process Complete**')
