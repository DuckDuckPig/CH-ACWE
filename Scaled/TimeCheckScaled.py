#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  7 08:16:29 2023

@author: jgra
"""
# In[1]
# Import Libraries and Tools
import pandas as pd
import numpy as np
import glob
import os

# In[2]
# Key Varibles

# Dataset
CR = 'CR2099'
StandardFolder = '../Standard/'

# Time File titles
fullScaleFile = CR + '_timeACWEresize_param1.csv'
timePrefix = '_time'

# In[3]
# Prepare to Cycle Through Files

# Find Remaining Files
scaleFiles = sorted(glob.glob(CR+timePrefix+'*'))[1:] 
scaleFiles.append( glob.glob(StandardFolder + CR + timePrefix + '*')[0])

# Open Inital File
fullScaleData = pd.read_csv(fullScaleFile,header=0)
fullScaleKeys = fullScaleData.keys()
fullScaleTimeList = fullScaleData[fullScaleKeys[1]]
print('Time, seconds')
print('   ',os.path.basename(fullScaleFile))
print('        Min: ',np.min(fullScaleTimeList))
print('        Mean:',np.mean(fullScaleTimeList),'±',np.std(fullScaleTimeList))
print('        Max: ',np.max(fullScaleTimeList))
print()

# In[4]
# Cycle Through Files
for scaleFile in scaleFiles:
    
    # Open File
    scaleData = pd.read_csv(scaleFile,header=0)
    scaleKeys = scaleData.keys()
    scaleTimeList = scaleData[scaleKeys[1]]
    
    print('   ',os.path.basename(scaleFile))
    print('        Min: ',np.min(scaleTimeList))
    print('        Mean:',np.mean(scaleTimeList),'±',np.std(scaleTimeList))
    print('        Max: ',np.max(scaleTimeList))
    print()
    print('        time mean(old/new) ± std(old/new):',np.mean(fullScaleTimeList/scaleTimeList),
          '±',np.std(fullScaleTimeList/scaleTimeList))
    print()
