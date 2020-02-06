#!/usr/bin/env python
# coding: utf-8

# # Homework 2 (Re-do)

# ParticleProperties.py

# In[17]:


'''
Redo of Homework 2, 
ParticleProperties
Farah Fauzi
'''

import numpy as np
import astropy.units as u
from ReadFile import Read                    #to read the file again

def ParticleInfo(readFile,p_type,p_num):
    time, total, data = Read("MW_000.txt")   #read in the file

    index = np.where(data['type']==p_type)   #array to store indexes of particles of desired p_type
    
    mag_dist = np.sqrt((data['x'][index][p_num-1])**2 + (data['y'][index][p_num-1])**2 + (data['z'][index][p_num-1])**2)*u.kpc
    dist = np.around(mag_dist,3)             #magnitude of distance is rounded to 3 decimal points
    print "Distance magnitude: ", dist      #print mag of distance in 3 decimal points

    mag_vel = np.sqrt((data['vx'][index][p_num-1])**2 + (data['vy'][index][p_num-1])**2 + (data['vz'][index][p_num-1])**2)*u.km/u.s
    vel = np.around(mag_vel,3)               #magnitude of velocity is rounded to 3 decimal points
    print "Velocity magnitude: ", vel        #print mag of velocity in 3 decimal points

    mass = data['m'][index][p_num-1]*1e10*u.M_sun
    print "Mass of the particles: ",mass       #print mass in the unit of Msun/1e10

    dist_new = mag_dist.to(u.lyr)
    dist_new = np.around(dist_new,3)
    print "Distance magnitude in lightyear: ",dist_new

ParticleInfo("MW_000.txt",2.00000,100)


# In[ ]:




