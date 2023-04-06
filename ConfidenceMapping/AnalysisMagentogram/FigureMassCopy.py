#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  9 14:21:16 2023

@author: jgra
"""

# In[1]
# Import Libraries and Tools
import glob
import os
import shutil

# In[2]
# Key Varibles

figureFolder = '/mnt/coronal_holes/Code Paper I Observations/ConMapStandardDefaultSkew/FigureVersions/'
copyFolder   = '/mnt/coronal_holes/Code Paper I Observations/ConMapStandardDefaultSkew/FigureVersionsCenter/'

CR = 'CR2133'

fileSuffex = '.png'

verbose = True

# In[3]
# prepare for coping

# Get file names
files = sorted(glob.glob(figureFolder + CR + '/*/*' + fileSuffex))

# Make sure folder exists
if not os.path.exists(copyFolder + CR + '/'):
    os.mkdir(copyFolder + CR + '/')

# In[4]
# Move files
for file in files:
    
    if verbose:
        print('Copying',os.path.basename(file))
    
    newLocation = copyFolder + CR + '/' + os.path.basename(file)
    
    shutil.copy2(file,newLocation)
    
# In[5]
# Inform User
print('**Process Complete**')