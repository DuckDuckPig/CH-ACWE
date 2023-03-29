#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Description:
    Global Consistency Error and Local Consistency Error as directly defined 
    in [1].

Reference:
[1] D. Martin, C. Fowlkes, D. Tal, and J. Malik, "A database of human 
    segmented natural images and its application to evaluating segmentation 
    algorithms and measuring ecological statistics," in Proceedings Eighth 
    IEEE International conference on Computer Vision. ICCV 2001, vol. 2, 2001, 
    pp. 416-423 vol.2.
Created on Fri Aug  7 11:05:12 2020

@author: jgra
"""

# In[1]:
# Import Libraries and tools
import numpy as np

# In[2]:
# Define Local Error
def E(S1,S2,p):
    
    # Define version of image that only retains the set of interest
    R1  = np.zeros(S1.shape).astype(bool)
    R2c = np.zeros(S2.shape).astype(bool)
    R1[np.where(S1 == S1[p[0],p[1]])]  = True # Only keep set that contains the pixel at p
    R2c[np.where(S2 != S2[p[0],p[1]])] = True # Only keep complement of set that contains pixel at p
    
    # Calculate numerator and denominator of Local Error
    num = np.count_nonzero(np.logical_and(R1,R2c))
    den = np.count_nonzero(R1)
    
    # Return Local Error
    return float(num)/den

# In[3]:
# Global Consistency Error
def GCE(S1,S2):
    
    # Prepare for Loop
    E1 = 0
    E2 = 0
    n  = len(S1) * len(S1[0])
    
    # Calculate and dum local errors
    for i in range(len(S1)):
        for j in range(len(S1[0])):
            E1 += E(S1,S2,[i,j])
            E2 += E(S2,S1,[i,j])
            
    # Return GCE
    return 1.0/n * np.min([E1,E2])

# In[4]:
# Local Consistency Error
def LCE(S1,S2):
    
    # Prepare for Loop
    Eout = 0
    n    = len(S1) * len(S1[0])
    
    # Calculate local error and sum minimum of the two
    for i in range(len(S1)):
        for j in range(len(S1[0])):
            E1 = E(S1,S2,[i,j])
            E2 = E(S2,S1,[i,j])
            Eout += np.min([E1,E2])
            
    # Return LCE
    return 1.0/n * Eout