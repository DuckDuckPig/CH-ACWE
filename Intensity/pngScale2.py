
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Description:
    Rescale Image to fit in a uint8 array
    
Created on Wed Jun  3 11:59:12 2020
Edited on  Wed Mar  8 14:34:47 2023

@author: jgra
"""

# In[1]:
# Inport Libraries and Tools
import numpy as np
import copy

# In[2]:
def scale(Image, scaleType = 'linear', scaleMin = 20.0, scaleMax = 2500.0):
    '''
    Clips a 2D array donw to the "scaleMin" and "scaleMax" values then performs
    either linear or log10 scalining before projecting the image into uint8.
    
    Inputs:
        Image : [float]
            Array to be scaled
        scaleType : str
            "linear" - each value maps directly - ie a range of 20 to 2500 
            (default inputs) will now be scaled and rounded to fit in a range
            of 0 to 255.
            "log10" - image will be cliped to linear range defined by scaleMin
            and ScaleMax. The result is then log10 scaled, crushing the high 
            values and stretching the low values. This result is then mapped
            to a range of 0 to 255. This image type if often used when 
            visulizing solar data, due to how bright the highest values can be.
        scaleMin : float
            Lower cuttoff for cliping image prior to scalining. Default value 
            20.0. Note that for log10 this function will yield invalid inputs
            if scaleMin <= 0.
        scaleMax : float
            Upper cuttoff for clipping image prior to scalining. Default
            value 2500.0
    Output:
        [uint8]
            Image rescaled (loosy).
    '''
    bitRange = 255
    iType    = 'uint8'
    
    # Copy Image
    I = copy.deepcopy(Image)
    
    # Clip to specified Range
    I = np.clip(I,scaleMin,scaleMax)
    
    # Log10 scale I, min, and max - if requested
    if scaleType == 'log10':
        if scaleMin <= 0: # Prevent log10(0) error
            offset = 1 - scaleMin
        else:
            offset = 0
        I = I + offset
        I   = np.log10(I)
        Min = np.log10(scaleMin + offset)
        Max = np.log10(scaleMax + offset)
    else:
        Min = scaleMin * 1
        Max = scaleMax * 1
        
    # Scale Image
    I = I - Min                   # Shift so min is 0
    I = I/(float(Max - Min))      # Scale 0-1
    I = np.round(I * bitRange)    # Rescale to full range
    I = np.clip(I,0,bitRange)     # Insure no overflow
    
    # Convert
    I = I.astype(iType)
    
    # Return Result
    return I

def unScale(Image, scaleType = 'Linear', scaleMin = 20.0, scaleMax = 2500.0):
    '''
    Takes a 2D array rescaled to fit in uint 8 and attempts to restore the 
    image as best as possible given the lossy nature of the compression.
    
    Inputs:
        Image : [float]
            Array that was scaled
        scaleType : str
            "linear" - each value was maped directly - ie a range of 20 to 2500 
            (default inputs) was scaled and rounded to fit in a range of 0 to 
            255.
            "log10" - image was cliped to linear range defined by scaleMin and 
            ScaleMax. The result was then log10 scaled, crushing the high 
            values and stretching the low values. This result was then mapped
            to a range of 0 to 255. This image type if often used when 
            visulizing solar data, due to how bright the highest values can be.
        scaleMin : float
            Lower cuttof used for clipping image prior to scaling, any value
            below this cuttoff can neither be resored nor approimated. Note 
            that for log10 this function will yield invalid inputs if 
            scaleMin <= 0.
        scaleMax : float
            Upper cuttoff used for clipping image prior to scalining. Default
            value 2500.0.
    Output:
        [uint8]
            Approimation of Image in original scale.
    '''
    bitRange = 255
    
    # Copy Image
    I = Image * 1.0
    
    # Find pre-scale min and max
    if scaleType == 'log10':
        if scaleMin <= 0: # account for adjustement made to prevent log10(0) error
            offset = 1 - scaleMin
        else:
            offset = 0
        Min = np.log10(scaleMin + offset)
        Max = np.log10(scaleMax + offset)
    else:
        Min = scaleMin * 1
        Max = scaleMax * 1
        
        
    # unscale
    I = I / float(bitRange)  # Sclae 0-1
    I = I * (Max-Min)        # Restore to full range
    I = I + Min              # Unshift
    if scaleType == 'log10': # Undo Log10 Scaling - if requessted
        I = np.power(10,I)
        I = I - offset
    
    # Return Result
    return I
