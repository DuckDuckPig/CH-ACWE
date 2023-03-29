#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Description:
    Compare segmentation of the same EUV image taken at different spatial 
    resolutions.

Created on Wed Jul 21 09:06:36 2021

@author: jgra
"""

# In[1]:
# Import Libraries and Tools
import os
import sys
import pandas as pd
import numpy as np
import glob
from skimage.metrics import structural_similarity as ssim

# ACWE utilities
# Root directory of the project
ROOT_DIR = os.path.abspath("../../")

# Import ACWE Tools
sys.path.append(ROOT_DIR)
from ACWE_python_spring_2023 import acweSaveSeg_v5, acweRestoreScale
from Metrics import consistancyErrorMetricsIII as cem, JaccardIndexMetric as jim

# In[2]:
# Key Variables

# Dataset
traceFolder = os.path.join(ROOT_DIR, 'DatasetTools/DownloadLists/')
CR = 'CR2133'  # Update to reflect chosen Carrington Rotation

# ACWE Folders
saveFolderScaled  = '/mnt/coronal_holes/Code Paper I Observations/Scaled/'   # Update to reflect Dataset Location
saveFolderDefault = '/mnt/coronal_holes/Code Paper I Observations/Standard/' # Update to reflect Dataset Location

# ACWE Compare Prefix
coreACWEprefix = 'ACWEresize_param1.' # Prefix for ACWE which we
                                      # will compare all others to
otherPrefix = ['ACWEresize_param2.','ACWEresize_param4.','ACWE.']

# ACWE Parameters 
acweChoice = '193'

# Inform user
verbose = True

# Resize Parameters
split = 0.5 # Breakpoint for upscaling

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
crSaveFolder1 = saveFolderScaled  + CR + '/'
crSaveFolder2 = saveFolderDefault + CR + '/'

# list of upscale methods
upscale = ['Nearest-neighbor','Bi-linear','Bi-quadratic','Bi-cubic','Bi-quartic','Bi-quintic']

# Create placeholder for stats
outputShape = [len(data),len(otherPrefix),len(upscale)]
IOU  = np.empty(outputShape); IOU[:]  = np.nan
SSIM = np.empty(outputShape); SSIM[:] = np.nan
GCE  = np.empty(outputShape); GCE[:]  = np.nan
LCE  = np.empty(outputShape); LCE[:]  = np.nan
alphas = np.empty([len(data),len(otherPrefix)+1,2]);alphas[:] = np.nan

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
        
    # Find all images
    acweFolder = file.split('/')[0] + '/'
    ACWEfiles = [sorted(glob.glob(crSaveFolder1+acweFolder+'*'+os.path.basename(file)+'*')),
                 glob.glob(crSaveFolder2+acweFolder+'*'+os.path.basename(file)+'*')]
    ACWEfiles = np.hstack(ACWEfiles)
    
    # Find coreACWEprefix
    coreACWE = None
    for ACWEfile in ACWEfiles:
        if coreACWEprefix in ACWEfile:
            coreACWE = ACWEfile
    
    # Compare Images
    if coreACWE != None:
        
        # Open File
        _,AHC,SegC = acweSaveSeg_v5.openSeg(coreACWE)
        cResize = AHC['RESIZE_PARAM']
        
        # Check for Errors
        alphas[i,0,0] = AHC['INIT_ALPHA']
        alphas[i,0,1] = AHC['ALPHA']
        
        # Cycle through ACWE Files
        for ACWEfile in ACWEfiles:
            
            # Organize
            for j in range(len(otherPrefix)):
                
                # Compare
                if otherPrefix[j] in ACWEfile:
                    
                    # Inform user
                    if verbose:
                        name = os.path.basename(ACWEfile).split('.')[0]
                        print('    Running on',name)
                    
                    # Open File
                    _,AH,Seg = acweSaveSeg_v5.openSeg(ACWEfile)
                    
                    # Check for Errors
                    alphas[i,j+1,0] = AH['INIT_ALPHA']
                    alphas[i,j+1,1] = AH['ALPHA']
                    
                    # Determin Resize Parameter
                    resizeParam = AH['RESIZE_PARAM']/cResize
                    
                    # Cycle Through Resize Methods
                    for m in range(len(upscale)):
                        
                        # Resize Image Using Specified Method
                        SegResized = acweRestoreScale.upscale(Seg,
                                                              {'RESIZE_PARAM':\
                                                               resizeParam},
                                                              upscale[m],
                                                              split,False)
                        SegResized = SegResized.astype(bool)
                        
                        # Calculate IOU
                        IOU[i,j,m] = jim.IOU(SegC,SegResized,True)
                        
                        # SSIM
                        SSIM[i,j,m],_ = ssim(SegResized,SegC,full=True)
                        
                        # GCE & LCE
                        GCE[i,j,m],LCE[i,j,m] = cem.CE(SegResized,SegC)
                        
# In[6]:
# Save Results for Graphing
resultsFolder = 'Scaled/Analysis/Results/'
resultsFolder = os.path.join(ROOT_DIR,resultsFolder)
if not os.path.exists(resultsFolder):
    os.mkdir(resultsFolder)
filename = os.path.basename(CarringtonFile) + '.ScaleStats.npz'
filename = resultsFolder + filename
np.savez(filename,alphas,IOU,SSIM,GCE,LCE)

# In[7]:
# End Program
print('**Process Complete**')
