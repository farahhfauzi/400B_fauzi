#!/usr/bin/env python
# coding: utf-8

# # Homework 2 (Re-do)

# ReadFile.py

# In[1]:


'''
Re-do of Homework2
Farah Fauzi
'''

import numpy as np
import astropy.units as u

def Read(filename):
    file = open(filename,'r')
 
    line1 = file.readline()                 #read the first line in MW_000.txt file
    label, value = line1.split()
    time = float(value)*u.Myr               #store the line read as time

    line2 = file.readline()                 #read the 2nd line in MW_000,txt file
    label, value = line2.split()
    total = float(value)                    #store the line read as total

    file.close()                            #close file

    data = np.genfromtxt(filename,dtype=None,names=True,skip_header=3) #store remainding data, skip first 3 rows

    #print (data['type'][1])                 #test print out particle type of 2nd particle

    return time, total, data

Read("MW_000.txt")


# In[ ]:




