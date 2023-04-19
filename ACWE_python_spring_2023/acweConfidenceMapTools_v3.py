#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 21 09:36:12 2023

@author: jgra
"""
# In[1]
# Import Libraries and Tools
import numpy as np
from . import acweRestoreScale

# In[2]
# Smart Combine Function (Recommended default)
def smartConMap(SEG,ACWEHEADER,buffer=0.05,normalize=True,restoreScale=True,
                interpolation='Bi-linear',split=0.5,returnInitMask=False,
                returnBackgroundWeights=False):
    """
    Generate Confidence map from input segmentation group. This function will
    attempt to recognize and remove Change of Target Cases

    Parameters
    ----------
    SEG : [float]
        ACWE Segmentations
    ACWEHEADER : dict
        ACWE Header, as developed by the saveSeg function.
    buffer : float, optional
        This buffer ensures that minute changes in the segmentation caused by
        the stopping criteria does not result in a non change of target case 
        being misidentified. This represents an acceptable drop in total
        area (normalized) of the original seed that can be excluded
        when dropping from one background/foreground ratio to the next lowest
        ratio. The default is 0.05.
    normalize : bool, optional
        Return the map with all values normalized to a range of 0 to 1. 
        The default is True.
    restoreScale : TYPE, optional
        Return confidence map that has been resized to the dimensions of the
        original EUV image. The default is True.
    interpolation : str, optional
        Interpolation method for the resizing process. Valid options are '
        Nearest-neighbor','Bi-linear','Bi-quadratic','Bi-cubic',
        'Bi-quartic', and 'Bi-quintic'. The default is 'Bi-linear'.
    split : float, optional
        Value above which all pixels in the upscaled image are assumed to be 
        part of the segmentation. The default is 0.5.
    returnInitMask : bool, optional
        Return an copy of the initial mask at the same scale. The default is 
        False.
    returnBackgroundWeights : bool, optional
        Return an copy of background weights that are part of the final 
        segmentation. The default is False.

    Returns
    -------
    ConMap : [float]
        ACWE confidence map, upscaled if requested.
    initMaks : [bool], optional
        Upscaled version of initial mask
    background_weights : [float], optional
        List of the background weights that were actually used
    """
    
    # Extract Inital Mask
    init_mask = ACWEHEADER['INIT_MASK']
    
    # Extract Background Weights
    background_weights = ACWEHEADER['BACKGROUND_WEIGHT']
    
    # Generate Ordered List of Background Weights
    background_weight_ordered = np.unique(np.sort(background_weights))
    
    # Prepare to generate Confidence Map
    SegNumber   = 0
    IOOlast     = 0
    indexList   = []
    outterBreak = False
    
    # Cycle through Ordered Background Weights
    for background_weight in background_weight_ordered:
        
        # Find and Check
        index = np.where(background_weight==background_weights)[0]
        for i in index:
            
            # Intersection Over Original Mask
            num = init_mask.astype(bool) & SEG[i].astype(bool)
            num = np.sum(num.astype(int))
            den = np.sum(init_mask.astype(int))
            IOOcurrent = float(num)/den
            
            # Check Change of Target (CoT) has Occurred
            # NOTE: A buffer is provided to account for minute changes
            #       due to variances caused by the stopping criteria
            if IOOcurrent + buffer <= IOOlast:
                
                # Ensure full break since Change of Target was found
                outterBreak = True
                break
            
            # Otherwise Add Segmentation to Confidence Mask
            SegNumber = SegNumber + 1
            indexList.append(i)
            IOOlast   = IOOcurrent * 1
            
        # Ensure full break at CoT
        if outterBreak:
            break
    
    # Generate Map
    ConMap             = SEG[indexList]
    background_weights = background_weights[indexList]
    
    # Upscale if requested
    if restoreScale:
        
        # Generate Micro Header with necessary Info
        resize_param = ACWEHEADER['RESIZE_PARAM']
        newHeader = {
                    'RESIZE_PARAM'             : resize_param,
                    'BACKGROUND_WEIGHT'        : background_weight,
                    'INIT_MASK'                : init_mask
                    }
        
        # Upscale
        ConMap,init_mask = acweRestoreScale.upscaleConMap(ConMap,newHeader,
                                                          interpolation,
                                                          split,True)
    # Combine Segmentations
    ConMap = np.sum(ConMap.astype(int),axis=0)
    
    # Normalize Map if User Requests
    if normalize:
        ConMap = ConMap / float(SegNumber)
    
    # Return Confidence Map
    if returnInitMask and returnBackgroundWeights:
        return ConMap,init_mask,background_weights
    elif returnBackgroundWeights:
        return ConMap,background_weights
    elif returnInitMask:
        return ConMap,init_mask
    else:
        return ConMap

# In[3]
# Simple Combine Function
def conMapCombine(SEG,ACWEHEADER,normalize=True,restoreScale=True,
                  interpolation='Bi-linear',split=0.5,returnInitMask=False):
    """
    Generate Confidence map without accounting for Change of Target.

    Parameters
    ----------
    SEG : [float]
        ACWE Segmentations
    ACWEHEADER : dict
        ACWE Header, as developed by the saveSeg function.
    normalize : bool, optional
        Return the map with all values normalized to a range of 0 to 1. 
        The default is True.
    restoreScale : TYPE, optional
        Return confidence map that has been resized to the dimensions of the
        original EUV image. The default is True.
    interpolation : str, optional
        Interpolation method for the resizing process. Valid options are '
        Nearest-neighbor','Bi-linear','Bi-quadratic','Bi-cubic',
        'Bi-quartic', and 'Bi-quintic'. The default is 'Bi-linear'.
    split : float, optional
        Value above which all pixels in the upscaled image are assumed to be 
        part of the segmentation. The default is 0.5.
    returnInitMask : bool, optional
        Return an copy of the initial mask at the same scale. The default is 
        False.

    Returns
    -------
    ConMap : [float]
        ACWE confidence map, upscaled if requested.
    initMaks : [bool], optional
        Upscaled version of initial mask

    """
    
    # Upscale if requested
    if restoreScale:
        
        # Upscale
        ConMap,init_mask = acweRestoreScale.upscaleConMap(SEG,ACWEHEADER,
                                                          interpolation,
                                                          split,True)
        
        # Combine
        ConMap = np.sum(ConMap,axis=0)
    
    else:
        init_mask = ACWEHEADER['INIT_MASK']
        ConMap    = np.sum(SEG,axis=0)
    
    # Normalize Map if User Requests
    if normalize:
        ConMap = ConMap / float(len(ACWEHEADER['BACKGROUND_WEIGHT']))
        
    # Return Results
    if returnInitMask:
        return ConMap,init_mask
    else:
        return ConMap
