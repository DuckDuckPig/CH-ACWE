#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 20 10:57:57 2023

@author: jgra
"""

# In[1]
# Import libraries and tools
import pandas as pd
import DataManagmentTools as dmt

# In[2]
# Dataset

traceFolder = 'DownloadLists/'
CR          = 'CR2133'

# Record Time
time = 'T_REC'

# In[3]

# Open Carrington Rotation Document
CarringtonFile = traceFolder + CR + '.csv'
data = pd.read_csv(CarringtonFile,header=0)
keys = data.keys()

# Get times
times = data[time]

# In[4]
# Determine Largest gap within dataset

# Cycle through times
for i in range(len(times)-1):
    
    # Time Data
    time1 = times[i]
    time2 = times[i+1]
    
    # Convert to datetime object
    timeobj1 = dmt.timeUnformat(time1)
    timeobj2 = dmt.timeUnformat(time2)
 
    # Calculate Time Difference/Gap
    gap = timeobj2 - timeobj1
    
    # Update Biggest Gap
    if i == 0:
        biggestGap = gap * 1
    elif gap > biggestGap:
        biggestGap = gap * 1
        

# In[5]
# Inform User
print('Largest time gap in',CR,':',biggestGap)
