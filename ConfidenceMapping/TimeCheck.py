#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  2 10:55:52 2023

@author: jgra
"""
# In[1]
# Import Libraries and Tools
import pandas as pd
import numpy as np

# In[2]
# Key Varibles

# Dataset
CR = 'CR2099'

# Confidence Map titles
oldDataFile = CR + '_timeOld.csv'
newDataFile = CR + '_timeDefault.csv'

# In[3]
# Open files
oldData = pd.read_csv(oldDataFile,header=0)
oldKeys = oldData.keys()
oldDataTimeList = oldData[oldKeys[1]]

newData = pd.read_csv(newDataFile,header=0)
newKeys = newData.keys()
newDataTimeList = newData[newKeys[1]]

# In[4]
# Inform User
print('Time (mean ± std)')
print('   Old Method:',np.mean(oldDataTimeList),'±',np.std(oldDataTimeList))
print('   New Method:',np.mean(newDataTimeList),'±',np.std(newDataTimeList))
print('   mean(old)/mean(new):',np.mean(oldDataTimeList)/np.mean(newDataTimeList))
print()
print('time mean(old/new) ± std(old/new):',np.mean(oldDataTimeList/newDataTimeList),
      '±',np.std(oldDataTimeList/newDataTimeList))