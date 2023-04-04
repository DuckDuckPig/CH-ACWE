#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Description:
    Explores the similarity of segmentations over short time periods by 
    comparing segmentations within +- 12 hours of the target time.
Created on Thu Oct 21 10:58:24 2021

@author: jgra
"""
# In[1]:
# Import Libraries and Tools
import os
import sys
import datetime
import numpy as np
import pandas as pd
from skimage.metrics import structural_similarity as ssim
import warnings
warnings.filterwarnings('ignore')

# taken from repojected_map
from reproject import reproject_interp

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS

import sunpy.map
from sunpy.coordinates import Helioprojective, RotatedSunFrame, transform_with_sun_center

# ACWE utilities
# Root directory of the project
ROOT_DIR = os.path.abspath("../../")

# Import ACWE Tools
sys.path.append(ROOT_DIR)
from ACWE_python_spring_2023 import acweSaveSeg_v5, acweRestoreScale
from Metrics import consistancyErrorMetricsIII as cem, JaccardIndexMetric as jim
from DatasetTools import DataManagmentTools as dmt

# In[2]:
# Key Variables

# Dataset
traceFolder = os.path.join(ROOT_DIR, 'DatasetTools/DownloadLists/')
CR = 'CR2133'

# ACWE Folder
saveFolder = '/mnt/coronal_holes/Code Paper I Observations/Standard/'

# ACWE Parameters 
prefix = 'ACWE.'
acweChoice = '193'

# Inform user
verbose = True

# Temporal summary
past   = 12 # hours
future = 12 # hours

# Resize Parameters
interpolation = 'Bi-linear' # Upscale using this interpolation
split = 0.5                 # Breakpoint for upscaling

# Backup rate
cycleLength = 5

# In[3]:
# Key Functions

def reprject(aiamap1,aiamap2):
    
    # Get Times
    t1 = aiamap1.date
    t2 = aiamap2.date
    
    # Define output frames
    out_frame2 = Helioprojective(observer='earth', obstime=t1)
    
    # Use output frame to define the "Target Frame"
    rot_frame2 = RotatedSunFrame(base=out_frame2, rotated_time=t2)
    
    # Construct WCS objects for the output map
    out_shape2 = aiamap1.data.shape # match shape of aiamap1
    out_center2 = SkyCoord(0*u.arcsec, 0*u.arcsec, frame=out_frame2)
    header2 = sunpy.map.make_fitswcs_header(out_shape2, out_center2,
                                            scale=u.Quantity(aiamap2.scale))
    out_wcs2 = WCS(header2)
    out_wcs2.coordinate_frame = rot_frame2
    
    # Perform Reprjection
    with transform_with_sun_center():
        arr2, _ = reproject_interp(aiamap2, out_wcs2, out_shape2)
        
    # Create the output map and preserve the original map's plot settings
    out_warp2 = sunpy.map.Map(arr2, out_wcs2)
    out_warp2.plot_settings = aiamap2.plot_settings
    
    # Extract rotated image anew
    image2Rot = out_warp2.data
    
    # Remove nans
    image2Rot[np.isnan(image2Rot)]=0
    
    # Convert back to map
    out2 = sunpy.map.Map(image2Rot, out_wcs2)
    out2.plot_settings = aiamap2.plot_settings
    
    # Return result
    return out2

# In[4]:
# Open file and get list of images

# Open Carrington Rotation Document
CarringtonFile = traceFolder + CR + '.csv'
data = pd.read_csv(CarringtonFile,header=0)
keys = data.keys()

# Find Image group we want 
acweChoice = np.where(keys == acweChoice)[0][0]

# In[5]:
# Prepare for analysis

# Create or find save location
crSaveFolder = saveFolder  + CR + '/'

# Final Save Folder and Temp File
resultsFolder = 'TemporalEffects/Analysis/Results/'
resultsFolder = os.path.join(ROOT_DIR,resultsFolder)
filename = os.path.basename(CarringtonFile) + '.TempFx.npz'
filename = resultsFolder + filename
progress = os.path.basename(CarringtonFile) + '.progress.txt'
progress = 'TemporalEffects/Analysis/' + progress
progress = os.path.join(ROOT_DIR,progress)

# Ensure Results Folder Exits
if not os.path.exists(resultsFolder):
    os.mkdir(resultsFolder)

# determin if starting or Restarting
# Case 1: Start Senario
if not (os.path.exists(filename) and os.path.exists(progress)):

    # Stats for final Graphs
    IOUpast    = np.empty([len(data),past]);   IOUpast[:,:]    = np.nan
    IOUfuture  = np.empty([len(data),future]); IOUfuture[:,:]  = np.nan
    SSIMpast   = np.empty([len(data),past]);   SSIMpast[:,:]   = np.nan
    SSIMfuture = np.empty([len(data),future]); SSIMfuture[:,:] = np.nan
    GCEpast    = np.empty([len(data),past]);   GCEpast[:,:]    = np.nan
    GCEfuture  = np.empty([len(data),future]); GCEfuture[:,:]  = np.nan
    LCEpast    = np.empty([len(data),past]);   LCEpast[:,:]    = np.nan
    LCEfuture  = np.empty([len(data),future]); LCEfuture[:,:]  = np.nan
    
    # Start = beginning of data
    start = 0

# Case 2: Restart Senario
else: 
    
    # Open data summary to reload results
    backupdata = np.load(filename, allow_pickle=True)
    lst = backupdata.files

    # Stats for final Graphs
    IOUpast    = backupdata[lst[0]]
    IOUfuture  = backupdata[lst[1]]
    SSIMpast   = backupdata[lst[2]]
    SSIMfuture = backupdata[lst[3]]
    GCEpast    = backupdata[lst[4]]
    GCEfuture  = backupdata[lst[5]]
    LCEpast    = backupdata[lst[6]]
    LCEfuture  = backupdata[lst[7]]
    
    # start = Last Backup + 1
    with open(progress,'r') as f:
        start = int(float(f.read()))
        print('*Restarting from Last Backup*')

# In[6]:
# Perform analysis

# If not complete
if start<len(data[keys[acweChoice]]):
        
    # Cycle Through Dataset
    for i in range(start,len(data[keys[acweChoice]])):
        
        # Inform user
        if verbose:
            print('Record:',data[keys[0]][i])
            
        # convert date to dateime object
        currentTime = dmt.timeUnformat(data[keys[0]][i])
        
        # Get file
        file = data[keys[acweChoice]][i]
        acweFolder = file.split('/')[0] + '/'
        file = acweFolder + prefix + os.path.basename(file) + '.npz'
        
        # Open Main Image
        H0,AH0,SEG0 = acweSaveSeg_v5.openSeg(crSaveFolder+file)
        
        # Proceed if no error
        if AH0['INIT_ALPHA'] == AH0['ALPHA']:
            
            # Upscale
            Seg0 = acweRestoreScale.upscale(SEG0,AH0,interpolation,split,False)
            Seg0 = Seg0.astype(bool)
            
            # Convert to Sunpy Map
            ACWEmap0 = sunpy.map.Map(Seg0.astype(int),H0)
            
            # Past Images
            for j in range(1,past+1):
                
                #Determin Valid Range (i.e. approx 1 hour ago, 2 hours ago, etc.)
                timeBound = currentTime - datetime.timedelta(hours=j)
                timeBound = [timeBound - datetime.timedelta(minutes=6),
                             timeBound + datetime.timedelta(minutes=6)]
                
                # Inform user
                if verbose:
                    if j > 1:
                        print('    Past',j,'hours')
                    else:
                        print('    Past 1 hour')
                
                # Cycle through possible images
                for k in range(1,past+2): # Assuming a 1 Hour Cadence
                    
                    # Valid Images
                    if i - k >=0:
                        pastTime = dmt.timeUnformat(data[keys[0]][i-k])
                        
                        # If Image is Found
                        if pastTime > timeBound[0] and pastTime < timeBound[1]:
                            
                            # Inform User
                            if verbose:
                                print('       ',data[keys[0]][i-k])
                            
                            # Get file
                            file = data[keys[acweChoice]][i-k]
                            acweFolder = file.split('/')[0] + '/'
                            file = acweFolder + prefix + os.path.basename(file) + '.npz'
                            
                            # Open Past Image
                            H1,AH1,SEG1 = acweSaveSeg_v5.openSeg(crSaveFolder+file)
                            
                            # Procede if no error
                            if AH1['INIT_ALPHA'] == AH1['ALPHA']:
                                
                                # Upscale
                                Seg1 = acweRestoreScale.upscale(SEG1,AH1,
                                                                interpolation,
                                                                split,False)
                                
                                # Convert to Sunpy Map
                                ACWEmap1 = sunpy.map.Map(Seg1.astype(int),H1)
                                
                                # Align Image 1
                                ACWEmap2 = reprject(ACWEmap0,ACWEmap1)
                                Seg2 = ACWEmap2.data
                                
                                # Convert back to integer
                                Seg2[np.where(Seg2>split)]=1
                                Seg2[np.where(Seg2!=1)]=0
                                Seg2 = Seg2.astype(bool)
                                
                                # Calulate IOU
                                IOU = jim.IOU(Seg2,Seg0,binary=True)
                                IOUpast[i,j-1] = IOU * 1
                                
                                # Calculate SSIM
                                SSIM,_ = ssim(Seg0.astype(bool),Seg2.astype(bool),full=True)
                                SSIMpast[i,j-1] = SSIM * 1
            
                                # Calcualte GCE and LCE
                                GCE,LCE = cem.CE(Seg0,Seg2)
                                GCEpast[i,j-1] = GCE * 1
                                LCEpast[i,j-1] = LCE * 1
                            
                            # move on to next value
                            break
            # Future Images
            for j in range(1,future+1):
                
                #Determine Valid Range (i.e. approx. 1 hour from now, 2 hours from now, etc.)
                timeBound = currentTime + datetime.timedelta(hours=j)
                timeBound = [timeBound - datetime.timedelta(minutes=6),
                             timeBound + datetime.timedelta(minutes=6)]
                
                # Inform user
                if verbose:
                    if j > 1:
                        print('    Future',j,'hours')
                    else:
                        print('    Future 1 hour')
                
                # Cycle through possible images
                for k in range(1,future+2): # Assuming a 1 Hour Cadence
                    
                    # Valid Images
                    if i + k < len(data):
                        futureTime = dmt.timeUnformat(data[keys[0]][i+k])
                        
                        # If Image is Found
                        if futureTime > timeBound[0] and futureTime < timeBound[1]:
                            
                            # Inform User
                            if verbose:
                                print('       ',data[keys[0]][i+k])
                            
                            # Get file
                            file = data[keys[acweChoice]][i+k]
                            acweFolder = file.split('/')[0] + '/'
                            file = acweFolder + prefix + os.path.basename(file) + '.npz'
                            
                            # Open Future Image
                            H1,AH1,SEG1 = acweSaveSeg_v5.openSeg(crSaveFolder+file)
                            
                            # Proceed if no error
                            if AH1['INIT_ALPHA'] == AH1['ALPHA']:
                                
                                # Upscale
                                Seg1 = acweRestoreScale.upscale(SEG1,AH1,
                                                                interpolation,
                                                                split,False)
                                
                                # Convert to Sunpy Map
                                ACWEmap1 = sunpy.map.Map(Seg1.astype(int),H1)
                                
                                # Align Image 1
                                ACWEmap2 = reprject(ACWEmap0,ACWEmap1)
                                Seg2 = ACWEmap2.data
                                
                                # Convert back to integer
                                Seg2[np.where(Seg2>split)]=1
                                Seg2[np.where(Seg2!=1)]=0
                                Seg2 = Seg2.astype(bool)
                                
                                # Calulate IOU
                                IOU = jim.IOU(Seg2,Seg0,binary=True)
                                IOUfuture[i,j-1] = IOU * 1
                                
                                # Calculate SSIM
                                SSIM,_ = ssim(Seg0.astype(bool),Seg2.astype(bool),full=True)
                                SSIMfuture[i,j-1] = SSIM * 1
            
                                # Calcualte GCE and LCE
                                GCE,LCE = cem.CE(Seg0,Seg2)
                                GCEfuture[i,j-1] = GCE * 1
                                LCEfuture[i,j-1] = LCE * 1
                            
                            # move on to next value
                            break
        # Backup
        if i % cycleLength == 0:
            
            # Inform User
            if verbose:
                print()
                print('Backing up results')
                
            # Save Backup
            np.savez(filename,IOUpast,IOUfuture,SSIMpast,SSIMfuture,GCEpast,
                     GCEfuture,LCEpast,LCEfuture)
            
            # Earmark
            with open(progress,'w') as f:
                start = f.write(str(i + 1))
        
        # Inform User
        if verbose:
            print()

# In[6]:
# Save Results for Graphing
np.savez(filename,IOUpast,IOUfuture,SSIMpast,
         SSIMfuture,GCEpast,GCEfuture,LCEpast,LCEfuture)

# Earmark
with open(progress,'w') as f:
    start = f.write(str(i))

# In[7]:
# End Program
print('**Process Complete**')
