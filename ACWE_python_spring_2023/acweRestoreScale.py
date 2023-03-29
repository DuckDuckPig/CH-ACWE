#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 22 10:45:44 2021

@author: jgra
"""
# In[1]:
# Import Libraries and Tools
import numpy as np
import skimage.transform

# In[2]:
# Define upscale functions

# Define upscale function
def upscaleConMap(SEG,ACWEHEADER,interpolation='Bi-linear',split=0.5,
                  returnInitMask=False):
    '''
    Function takes in the unaltered confidence map and the 
    ACWE header and uses this information to upscale all
    of the individual layers to the desired scale.
    
     Parameters
    ----------
    SEG : [float]
        ACWE segmentation or segmentations
    ACWEHEADER : dict
        ACWE Header, as developed by the saveSeg function
    interpolation : str
        Interpolation method for the resizing process. Valid options are '
        Nearest-neighbor','Bi-linear','Bi-quadratic','Bi-cubic',
        'Bi-quartic', and 'Bi-quintic'.
        
        Default Value: 'Bi-linear'
    split : float
        Value above which all pixels in the upscaled image are assumed to be 
        part of the segmentation.
        
        Default Value: 0.5
    returnInitMask : bool
        Return an upscaled copy of the initial mask. This can be used to 
        identify cases where ACWE has changed target from CH to QS
        
        Default Value: False
    Returns
    -------
    seg : [float]
        Copy of the ACWE segmentation or segmentations upscaled to the size  
        of the original .fits file.
    initMaks : [bool] (optional)
        Upscaled version of initial mask
    '''
    
    # Determine size of upscaled images
    shape = np.asarray(SEG.shape)
    shape[1:] = shape[1:] * ACWEHEADER['RESIZE_PARAM']
    
    # Determin Upscale Method
    upscale = ['Nearest-neighbor','Bi-linear','Bi-quadratic','Bi-cubic',
               'Bi-quartic','Bi-quintic']
    for m in range(len(upscale)):
        if upscale[m] == interpolation:
            break
        
    # Generate placeholder Segmentation
    seg = np.empty(shape); seg[:,:,:] = np.nan
    
    # Populate
    for i in range(shape[0]):
        
        # valid Segmentation
        if np.sum(np.isnan(SEG[i]).astype(int)) == 0:
            
            # upscale
            s = skimage.transform.resize(SEG[i].astype(int),shape[1:],order=m,
                                         preserve_range=True,
                                         anti_aliasing=True)
            
            # Convert back to integer
            s[np.where(s>split)]=1
            s[np.where(s!=1)]=0
            
            # Copy Upscaled Segmentation into correct location
            seg[i] = s*1
            
    # Initial Mask
    if returnInitMask:
        initMask = skimage.transform.resize(ACWEHEADER['INIT_MASK'].astype(int),
                                            shape,order=m,
                                            preserve_range=True,
                                            anti_aliasing=True)
        
        # Convert back to integer
        initMask[np.where(initMask>split)]=1
        initMask[np.where(initMask!=1)]=0
        
        # Return upscaled Segmentation
        return seg,initMask.astype(bool)
    
    else:
            
        # Return upscaled Segmentation
        return seg

def upscale(SEG,ACWEHEADER,interpolation='Bi-linear',split=0.5,
            returnInitMask=False):
    '''
    Function takes in the unaltered segmentation and the 
    ACWE header and uses this information to upscale the
    segmentation back to the size of the original .fits
    file.
    
    Parameters
    ----------
    SEG : [float]
        ACWE segmentation or segmentations
    ACWEHEADER : dict
        ACWE Header, as developed by the saveSeg function
    interpolation : str
        Interpolation method for the resizing process. Valid options are '
        Nearest-neighbor','Bi-linear','Bi-quadratic','Bi-cubic',
        'Bi-quartic', and 'Bi-quintic'.
        
        Default Value: 'Bi-linear'
    split : float
        Value above which all pixels in the upscaled image are assumed to be 
        part of the segmentation.
        
        Default Value: 0.5
    returnInitMask : bool
        Return an upscaled copy of the initial mask. This can be used to 
        identify cases where ACWE has changed target from CH to QS
        
        Default Value: False
    Returns
    -------
    seg : [float]
        Copy of the ACWE segmentation or segmentations upscaled to the size  
        of the original .fits file.
    initMaks : [bool] (optional)
        Upscaled version of initial mask
    '''
    
    # Determine size of upscaled images
    shape = np.asarray(SEG.shape)
    shape = shape * ACWEHEADER['RESIZE_PARAM']
    
    # Determine Upscale Method
    upscale = ['Nearest-neighbor','Bi-linear','Bi-quadratic','Bi-cubic',
               'Bi-quartic','Bi-quintic']
    for m in range(len(upscale)):
        if upscale[m] == interpolation:
            break
        
    # Upscale
    seg = skimage.transform.resize(SEG.astype(int),shape,order=m,
                                   preserve_range=True,anti_aliasing=True)
            
    # Convert back to integer
    seg[np.where(seg>split)]=1
    seg[np.where(seg!=1)]=0
            
    # Initial Mask
    if returnInitMask:
        initMask = skimage.transform.resize(ACWEHEADER['INIT_MASK'].astype(int),
                                            shape,order=m,
                                            preserve_range=True,
                                            anti_aliasing=True)
        
        # Convert back to integer
        initMask[np.where(initMask>split)]=1
        initMask[np.where(initMask!=1)]=0
        
        # Return upscaled Segmentation
        return seg,initMask.astype(bool)
    
    else:
            
        # Return upscaled Segmentation
        return seg
