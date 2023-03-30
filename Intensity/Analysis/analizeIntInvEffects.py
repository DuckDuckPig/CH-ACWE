#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Description: Return similiarty between original segmentation and segmetnations
             generated from reduced intensity dynmatic range.

Created on Mon Feb 21 15:37:57 2022
Edited on  Thu Mar  2 18:22:39 2023

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
from ACWE_python_spring_2023 import acweSaveSeg_v5
from Metrics import consistancyErrorMetricsIII as cem, JaccardIndexMetric as jim


# In[2]:
# Key Varibles

# Dataset
traceFolder = os.path.join(ROOT_DIR, 'DatasetTools/DownloadLists/')
CR = 'CR2133'

# ACWE Folders
saveFolderIntInv  = '/mnt/coronal_holes/Code Paper I Observations/IntInv/'
saveFolderDefault = '/mnt/coronal_holes/Code Paper I Observations/Standard/'

# ACWE Compare Prefix
coreACWEprefix = 'ACWE.' # Prefix for ACWE which we
                         # will compare all others to
otherPrefix = ["LinearCompressFull.",
               "LinearCompressFullRestored.",
               "LinearCompress0toMax.",
               "LinearCompress0toMaxRestored.",
               "LinearCompressSolarLimits.",
               "LinearCompressSolarLimitsRestored.",
               "LinearCompressDefault.",
               "LinearCompressDefaultRestored.",
               "Log10CompressFull.",
               "Log10CompressFullRestored.",
               "Log10Compress0toMax.",
               "Log10Compress0toMaxRestored.",
               "Log10CompressSolarLimits.",
               "Log10CompressSolarLimitsRestored.",
               "Log10CompressDefault.",
               "Log10CompressDefaultRestored."]

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
crSaveFolder1 = saveFolderIntInv  + CR + '/'
crSaveFolder2 = saveFolderDefault + CR + '/'

outputShape = [len(data),len(otherPrefix)]
IOU  = np.empty(outputShape); IOU[:]  = np.nan
SSIM = np.empty(outputShape); SSIM[:] = np.nan
GCE  = np.empty(outputShape); GCE[:]  = np.nan
LCE  = np.empty(outputShape); LCE[:]  = np.nan
sucess = np.empty(outputShape); sucess[:] = np.nan
alphas = np.empty(outputShape); alphas[:] = np.nan

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
                 sorted(glob.glob(crSaveFolder2+acweFolder+'*'+os.path.basename(file)+'*'))]
    ACWEfiles = np.hstack(ACWEfiles)
    
    # Find coreACWEprefix
    coreACWE = None
    for ACWEfile in ACWEfiles:
        if coreACWEprefix in ACWEfile:
            coreACWE = ACWEfile
    
    # Compare Images
    if coreACWE != None:
        
        # Open File
        _,_,SegC = acweSaveSeg_v5.openSeg(coreACWE)
        
        # Cycle through ACWE Files
        for ACWEfile in ACWEfiles:
            
            # Organize
            for j in range(len(otherPrefix)):
                
                # Compare
                if otherPrefix[j] in ACWEfile:
                    
                    # Inform user
                    if verbose:
                        name = os.path.basename(ACWEfile).split('.')[0]
                        print('    Running on',j,name)
                    
                    # Open File
                    _,AH,Seg = acweSaveSeg_v5.openSeg(ACWEfile)
                    
                    # Check for Errors
                    if np.sum(np.isnan(Seg).astype(int)) == 0:
                        
                        # Calculate IOU
                        IOU[i,j] = jim.IOU(SegC,Seg,True)
                        
                        # SSIM
                        SSIM[i,j],_ = ssim(Seg,SegC,full=True)
                        
                        # GCE & LCE
                        GCE[i,j],LCE[i,j] = cem.CE(Seg,SegC)
                        
                        # Save Sucess and Alpha Parameter
                        if np.sum(Seg) != 0:
                            # Earmark that this segmentation is not broken
                            sucess[i,j] = 1
                        else:
                            # Earmark that segmentation could not be generated
                            sucess[i,j] = 0.5
                        alphas[i,j] = AH['ALPHA']
                    
                    else: # Error Found
                        
                        #  Earmark that this segmentaiton is broken
                        sucess[i,j] = 0
                        alphas[i,j] = AH['ALPHA']
                    
# In[6]:
# Save Results for Graphing
if not os.path.exists('Results/'):
    os.mkdir('Results/')
filename = 'Results/' + os.path.basename(CarringtonFile) + '.intInvStats.npz'
np.savez(filename,sucess,alphas,IOU,SSIM,GCE,LCE)

# In[7]:
# End Program
print('**Process Complete**')
