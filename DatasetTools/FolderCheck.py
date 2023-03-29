#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Description: Inform user if data are missing

Created on Thu Jul  8 09:11:31 2021

@author: jgra
"""
# In[1]:
# Import Libraries and Tools
import pandas as pd
import os

# In[2]:
# Key Variables

# Dataset folders
dataFolder  = '/home/jgra/Coronal Holes/newDataset/' # Specify folder path here.

# Dataset File
traceFile = 'DownloadLists/CR2099.csv' # Change to match 'CR*.csv; file name and location.

# In[3]:
# Read Dataset File

data = pd.read_csv(traceFile,header=0)
keys = data.keys()

# In[4]:
# Find Full File Path

# Format to match folder setup
crFolder = os.path.basename(traceFile).split('.')[0] # CR21001.csv contains
crFolder = dataFolder + crFolder + '/'               # data on folder CR2101

# In[5]:
# Check For Each Folder

for date in data[keys[1]]:
    
    # Format Date to match folder setup
    folder = date.replace(':','') # All observations for T_REC 2010-09-16T07:00:02Z
    folder = crFolder + folder    # are stored in folder 2010-09-16T070002Z
    
    # Inform user if file is missing
    if not os.path.exists(folder):
        print(os.path.basename(folder), 'Missing')

# In[6]:
# End Process

print('**Process Complete**')
