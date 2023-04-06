#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Description:
    Calculate skewness for each CH group as a function of confidence for each
    confidence map in the specified Carrington rotation. Results are saved in 
    a mirrored folder structure in the user-specified location as .npz files.
    the results can be viewed using the 'skew check*.ipynb' notebooks.
    
Created on Thu Dec 16 13:50:02 2021

@author: jgra
"""

# In[1]:
# Import Libraries and Tools
import os
import sys
import pandas as pd
import numpy as np
import skimage.morphology
import skimage.measure
from astropy.io import fits
import sunpy.map
from reproject import reproject_interp
import scipy.stats

import warnings
warnings.filterwarnings("ignore")

# ACWE utilities
# Root directory of the project
ROOT_DIR = os.path.abspath("../../")

# Import ACWE Tools
sys.path.append(ROOT_DIR)
from ACWE_python_spring_2023 import acweSaveSeg_v5, acweRestoreScale

# In[2]:
# Key Variables

# Dataset
# Dataset Folders
dataFolder  = '/home/jgra/Coronal Holes/newDataset/'
traceFolder = os.path.join(ROOT_DIR, 'DatasetTools/DownloadLists/')
CR          = 'CR2101'

# SaveFolder
saveFolder = '/mnt/data/jgraCoronalHoles/CodeVI Observations/ConMapStandardDefault_5_1/'
skewFolder = '/mnt/data/jgraCoronalHoles/CodeVI Observations/ConMapStandardDefaultSkew/'

# ACWE Choice
acwePrefix  = 'ACWEconMap.' # Prefix for ACWE
acweChoice  = '193'
magnetogram = 'Magnetogram'

# Save Result to
skewPrefix = 'Skewness.'
overwrite = False    # Run from begining 

# Inform User
verbose = True

# In[3]:
# Open File and Get List of Images

# Open Carrington Rotation Document
CarringtonFile = traceFolder + CR + '.csv'
data = pd.read_csv(CarringtonFile,header=0)
keys = data.keys()

# Find Image Group we Want 
acweChoice  = np.where(keys == acweChoice)[0][0]
magnetogram = np.where(keys == magnetogram)[0][0]

# In[4]:
# Prepare to Calculate Skewness

# Find Folders
crSaveFolder   = saveFolder + CR + '/'
crSourceFolder = dataFolder + CR + '/'
crSkewFolder   = skewFolder + CR + '/'

# Create any Missing Folders
if not os.path.exists(crSkewFolder):
    os.mkdir(crSkewFolder)

# In[5]:
# Calculate Skewness

# Cycle Through Dataset
for i in range(len(data)):#-1,-1,-1):
    
    # Inform User
    if verbose:
        print(data[keys[0]][i])
    
    # Get AIA File Name
    file = data[keys[acweChoice]][i]
    
    # Placement Folders
    acweFolder = file.split('/')[0] + '/'
    if not os.path.exists(crSkewFolder + acweFolder):
        os.mkdir(crSkewFolder + acweFolder)
    
    # ACWE Name
    acweFile = acwePrefix + os.path.basename(file)
    acweFile = acweFolder + acweFile + '.npz'
    
    # Result Name
    skewFile = skewPrefix + acwePrefix + os.path.basename(file)
    skewFile = acweFolder + skewFile + '.npz'
    
    # Check for File
    if overwrite or not os.path.exists(crSkewFolder + skewFile):
        
        # Inform User
        if verbose:
            print('    Opening and resizing',os.path.basename(acweFile))
        
        # Open Confidence Map
        H,AH,SEG = acweSaveSeg_v5.openSeg(crSaveFolder+acweFile)
        
        # Upscale
        SEG = acweRestoreScale.upscaleConMap(SEG, AH)
        
        # Combine
        SEG[np.where(np.isnan(SEG))]=0
        SEG = np.sum(SEG,axis=0)
        
        # Inform User
        if verbose:
            print('    Identifying Clusters')
        
        # Find and Label CH Regions
        ACWE_clusters = SEG * 1
        ACWE_clusters[np.where(ACWE_clusters>0)] = 1
        ACWE_clusters = skimage.morphology.dilation(ACWE_clusters.astype(bool),
                                                    selem=np.ones([40,40]))
        ACWE_clusters = skimage.measure.label(ACWE_clusters,connectivity=2)
        
        # Split CH Regions
        SEG2 = np.zeros([np.max(ACWE_clusters),len(SEG),len(SEG[0])])
        for j in range(1,np.max(ACWE_clusters)+1):
            cluster = ACWE_clusters * 1
            cluster[np.where(cluster!=j)] = 0
            cluster = cluster/np.max(cluster)
            SEG2[j-1] = SEG*cluster
        
        # Inform User
        if verbose:
            print('    Preparing Magnetogram')
        
        # Open Magnetogram
        hdulist = fits.open(crSourceFolder + data[keys[magnetogram]][i])
        hdulist.verify('silentfix') #necessary for successful data read
        h_mag = hdulist[1].header
        J_mag = hdulist[1].data
        hdulist.close()
        
        # Create Maps
        mapHMI = sunpy.map.Map((J_mag,h_mag))
        mapHMI.plot_settings['cmap']='hmimag'
        mapACWE = sunpy.map.Map(SEG,H)
        
        # Reproject Magnetogram
        hmiReproject,footprint = reproject_interp(mapHMI,mapACWE.wcs,
                                                  mapACWE.data.shape)
        
        # Generate Weights to Address Projection Effects
        # Based on: https://docs.sunpy.org/en/stable/generated/gallery/map_transformations/reprojection_aia_euvi_mosaic.html#improving-the-output
        Weights = sunpy.map.all_coordinates_from_map(mapACWE)
        Weights = Weights.transform_to('heliocentric').z.value
        Weights = (Weights/np.nanmax(Weights)) ** 3
    
        # Make Weighted Magentogram Map
        hmiReprojectWeighted = hmiReproject * Weights
        
        # Inform User
        if verbose:
            print('    Calculating Skewness')
        
        # Prepare to Calculate Skewness
        skew = np.empty([len(SEG2),2,len(AH['BACKGROUND_WEIGHT'])])
        skew[:] = np.nan
        
        # Cycle Through Regions
        for j in range(len(SEG2)):
            
            # Skewness - Both Unweighted (M) and Weighted (W)
            skwM = []
            skwW = []
            
            # Cycle Through Confidence Levels
            for k in range(len(AH['BACKGROUND_WEIGHT'])):
                
                # Threshold by Confidence
                CH_kept = SEG2[j]*1
                CH_kept[np.where(CH_kept < k)] = 0 
                CH_kept = CH_kept.astype(bool)
                
                # Flatten Selected Region
                flattenNormal = hmiReproject[CH_kept].flatten()
                flattenWeight = hmiReprojectWeighted[CH_kept].flatten()
                
                # Remove NAN Values, if any Persist
                flattenNormal = flattenNormal[np.logical_not(np.isnan(flattenNormal))]
                flattenWeight = flattenWeight[np.logical_not(np.isnan(flattenWeight))]
                
                # Calculate Skewness of CH Regions
                skwM.append(scipy.stats.skew(flattenNormal))
                skwW.append(scipy.stats.skew(flattenWeight))
            
            # Store Results
            skew[j,0] = skwM # Unweighted Skewness
            skew[j,1] = skwW # Weighted Skewness
        
        # Inform User
        if verbose:
            print('    Saving Results')
        
        # Save Results
        backgroundWeights = np.ones(len(SEG2)) * len(AH['BACKGROUND_WEIGHT'])
        np.savez_compressed(crSkewFolder+skewFile,SEG2,skew,backgroundWeights)
    
# In[6]:
# End Process

print('**Process Complete**')
