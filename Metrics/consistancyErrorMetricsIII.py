#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Description:
    Global Consistency Error and Local Consistency Error from [1] defined
    using an optimized implementation.

Reference:
[1] D. Martin, C. Fowlkes, D. Tal, and J. Malik, "A database of human 
    segmented natural images and its application to evaluating segmentation 
    algorithms and measuring ecological statistics," in Proceedings Eighth 
    IEEE International conference on Computer Vision. ICCV 2001, vol. 2, 2001, 
    pp. 416-423 vol.2.

Created on Wed Aug 12 11:12:56 2020

@author: jgra
"""

import numpy as np

def CE(S1,S2):
    '''
    Calculate and return both global consistency error (GCE) and local 
    consistency error (LCE).
    
    Parameters
    ----------
    S1 : [float]
        Segmentation 1
    S2 : [float]
        Segmentation 2
    Returns
    -------
    GCE : float
        Global consistency error
    LCE : float
        Local consistency error
    '''
    
    # Key Variables
    E1 = [] # holds all local errors in direction of S1 to S2
    E2 = [] # holds all local errors in direction of S2 to S1
    N = len(S1) * len(S1[0])
    
    # Determine number of regions in each segmentation
    s1_span = range(np.min(S1),np.max(S1)+1)
    s2_span = range(np.min(S2),np.max(S2)+1)
    
    # Compare each segmentation head on
    for i in s1_span:
        for j in s2_span:
            # Calculate cardinality of intersect
            intersect = np.sum(np.logical_and(S1==i,S2==j))
            # determine if local error is a valid measure
            if intersect > 0:
                # Calculate cardinality of set difference
                E1.append(np.sum(np.logical_and(S1==i,S2!=j)))
                E2.append(np.sum(np.logical_and(S2==j,S1!=i)))
                # Normalize by cardinaity of region
                E1[-1] = E1[-1]/float(np.sum(S1==i))
                E2[-1] = E2[-1]/float(np.sum(S2==j))
                # Expand by cardianity of intersect
                E1[-1] = E1[-1] * intersect
                E2[-1] = E2[-1] * intersect
    
    # Vstack E1 & E2
    E = np.vstack([E1,E2])
    
    # Determin GCE
    GCE = np.sum(E,axis=1) # Sum over all local error
    GCE = np.min(GCE)      # Take the min of the two values
    GCE = GCE/float(N)     # Normalize by size of image
    
    # Determin LCE
    LCE = np.min(E,axis=0) # Take min of each local error
    LCE = np.sum(LCE)      # Sum result
    LCE = LCE/float(N)     # Normalize by size of image
    
    # Return result
    return GCE, LCE
        
def GCE(S1,S2):
    '''
    Return global consistency error (GCE) for two segmentations
    
    Parameters
    ----------
    S1 : [float]
        Segmentation 1
    S2 : [float]
        Segmentation 2
    Returns
    -------
    GCE : float
        Global consistency error
    '''
    
    GCE,_=CE(S1,S2)
    
    return GCE

def LCE(S1,S2):
    '''
    Return local consistency error (LCE) for two segmentations
    
    Parameters
    ----------
    S1 : [float]
        Segmentation 1
    S2 : [float]
        Segmentation 2
    Returns
    -------
    LCE : float
        Local consistency error
    '''
    
    _,LCE=CE(S1,S2)
    
    return LCE
