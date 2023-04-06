#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 20 13:46:48 2023

@author: jgra
"""

# In[1]
# Import key libaries and tools
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import glob
from astropy.io import fits
import sunpy.map
from aiapy.calibrate import register, update_pointing

plt.rcParams.update({'font.size': 18})

import warnings
warnings.filterwarnings("ignore")

# ACWE utilities
# Root directory of the project
ROOT_DIR = os.path.abspath("../../")

# Import ACWE Tools
sys.path.append(ROOT_DIR)
from ACWE_python_spring_2023 import acweSaveSeg_v5

# In[2]
# Key Varibles
sFolder = '/mnt/coronal_holes/Code Paper I Observations/ConMapStandardDefaultSkew/'
cFolder = '/mnt/coronal_holes/Code Paper I Observations/ConMapStandardDefault/'
aFolder = '/home/jgra/Coronal Holes/newDataset/'
CR      = 'CR2133'
files   = sorted(glob.glob(sFolder + CR + '/*/*.npz'))
conMaps = sorted(glob.glob(cFolder + CR + '/*/*.npz'))
AIAs    = sorted(glob.glob(aFolder + CR + '/*/*.193.image_lev1.fits'))
saveTo  = '/mnt/coronal_holes/Code Paper I Observations/ConMapStandardDefaultSkew/FigureVersions/'

# In[3]
# Sanity Check

print()
print(CR)
print('Skew Data',len(files))
print('Confidence Maps',len(conMaps))
print('EUV Files',len(AIAs))
print()
print('Converting Skew Data to Figures:')

# In[4]
# Functions

# Circle Mask Function
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

# Open File and Prepare for ACWE
def openAIA(filename):
    
    # Extract Image and Header Data
    hdulist = fits.open(filename)
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
        
        # Prepare Display Version
    Idsp = np.clip(I,20,2500)
    Idsp = np.log10(Idsp)
    Idsp = Idsp - np.min(Idsp)
    Idsp = Idsp/np.max(Idsp)
    
    return I,Idsp,H

# In[5]:
for i in range(len(files)):
    
    # Find Files
    file   = files[i]
    conMap = conMaps[i]
    AIA    = AIAs[i]
    
    # File Time
    fileTime = file.split('/')[-2]#.replace('-','')
    
    #  Retreve Data - Skew
    data = np.load(file, allow_pickle=True)
    lst = data.files

    SEG2 = data[lst[0]]
    Skew = data[lst[1]]
    BW   = data[lst[2]]
    
    HaveSEG = False

    # Plot Each CH Group
    for j in range(len(SEG2)):
        
        fileFolder = saveTo + CR + '/' + fileTime + '/'
        title = fileFolder + os.path.basename(file) + '.CH_Group' + str(j+1) + '.png'
        if not os.path.exists(title):
            
            if not HaveSEG:
                print('   ',fileTime)
                
                # Open AIA Image
                I,Idsp,H = openAIA(AIA)
                
                # Retreve Data - Origial Con Map
                H,AH,SEG = acweSaveSeg_v5.openSeg(conMap)
                
                # Determine characteristics of image
                im_size = np.asarray(SEG2[0].shape) # size of the image
                try: # SDO/AIA
                    sun_radius = H['R_SUN'] # solar radius from metadata
                except: # STERO A or B
                    sun_radius = (H['RSUN']/H['CDELT1']) # solar radius from metadata
                sun_center = np.asarray([int(round(H['CRPIX1']))-1,int(round(H['CRPIX2']))-1]) # solar center from metadata
    
                # Make Circle Mask
                sd_mask = make_circle_mask(sun_center,im_size,sun_radius)
            
                HaveSEG = True
        
            # Prepare Plot
            plt.figure(figsize=[80,20])
            plt.rcParams.update({'font.size': 60})
            plt.suptitle(os.path.basename(AIA))

            # Plot Image
            plt.subplot(1,3,1)
            plt.imshow(np.flip(Idsp,axis=0),cmap ='gray')
            plt.title('Observation')
            plt.axis(False)

            # Plot Confidence Map
            plt.subplot(1,3,2)
            plt.imshow(np.flip(SEG2[j]/BW[j],axis=0),vmin=0,vmax=1,interpolation='None')
            plt.colorbar(shrink=0.95)
            plt.contour(np.flip(sd_mask,axis=0),0,colors='r')
            plt.axis(False)
            s = 'CH Group ' + str(j+1) + ' Confidence Map'
            plt.rcParams.update({'font.size': 50})
            plt.title(s)

            # Plot Skewness, Weighted & Unweighted
            plt.subplot(1,3,3)
            xaxis = np.asarray(range(1,BW[j].astype(int)+1))/BW[j]
            plt.plot(xaxis,Skew[j,0],label='unweighted',linewidth=4)
            plt.plot(xaxis,Skew[j,1],label='weighted',linewidth=4)
            plt.xlabel('Confidence')
            plt.ylabel('Skew')
            s = 'CH Group ' + str(j+1) + ' Skew'
            plt.title(s)
            plt.legend()
            plt.grid()

            if not os.path.exists(saveTo+CR+'/'):
                os.makedirs(saveTo+CR+'/')
            if not os.path.exists(fileFolder):
                os.mkdir(fileFolder)
            title = fileFolder + os.path.basename(file) + '.CH_Group' + str(j+1) + '.png'
            plt.savefig(title)
            plt.close()

# In[6]
# End Process

print()
print('**Process Complete**')
