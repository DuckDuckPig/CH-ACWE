#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 13:52:26 2023

@author: jgra
"""

# In[1]
# Import Libraries and Tools
import numpy as np

# In[2]
# Intersection Over Union
def IOU(A,B,binary=True):
    """
    Return the jaccard Index or IOU of the two segmentations.

    Parameters
    ----------
    A : [bool] OR [float]
        Segmentation 1.
    B : [bool] OR [float]
        Segmentation 2.
    binary : TYPE, optional
        True indicates that the two maps are to be treated as binary 
        segmentations. When false the weighted IOU is calculated and returned
        instead. The default is True.

    Returns
    -------
    IOU : TYPE
        Intersection over union of the two images or segmentations.

    """
    
    # Binary IOU
    if binary:
        
        # Cast to Binary
        a = A.astype(bool)
        b = B.astype(bool)

        # Calculate IOU
        num = a & b
        num = np.sum(num.astype(int))
        den = a | b
        den = np.sum(den.astype(int))
        IOU = float(num)/den
    
    # Weighted IOU
    else:
        
        # Combine to optimize process
        stack = np.dstack([A,B])
        
        # Calculate weighted IOU
        num = np.sum(np.min(stack,axis=2))
        den = np.sum(np.max(stack,axis=2))
        IOU = float(num)/den
    
    # Return Results
    return IOU
