#!/usr/bin/env python
# coding: utf-8

# In[1]:


'''
Homework 3
Farah Fauzi
'''
import numpy as np
import astropy.units as u
from ReadFile import Read


# In[10]:


# function ComponentMass return total mass of galaxy component that takes input filename and particle type
# return mass in units of 1e12*Msun
def ComponentMass(ptype):
    MW = []                         #empty array for each galaxy
    M31 = []
    M33 = []
    
    #summing up mass of each component for different galaxy
    #read data from filename
    #append data into array for each galaxy
    #print total mass of component for each galaxy
    filename = "MW_000.txt"
    time, total, data = Read(filename)
    for p_type in range(1,4):
        index = np.where(data['type']==p_type)
        #get total mass for all p_type
        totMass1 = sum(data['m'][index]*1e-2*u.M_sun)
        newMass1 = np.round(totMass1,3)
        #get total mass for desired p_type for MW
        if (p_type==ptype):
            totMassH = sum(data['m'][index]*1e-2*u.M_sun)
            newMassH = np.round(totMassH,3)
        #insert total mass of component into MW array
        MW.append(newMass1)
        p_type+=1
        
    filename = "M31_000.txt"
    time, total, data = Read(filename)
    for p_type in range(1,4):
        index = np.where(data['type']==p_type)
        #get total mass for all p_type
        totMass2 = sum(data['m'][index]*1e-2*u.M_sun)
        newMass2 = np.round(totMass2,3)
        #get total mass for desired p_type for M31
        if (p_type==ptype):
            totMassD = sum(data['m'][index]*1e-2*u.M_sun)
            newMassD = np.round(totMassD,3)
        #insert total mass of component into M31 array
        M31.append(newMass2)
        p_type+=1
        
    filename = "M33_000.txt"
    time, total, data = Read(filename)
    for p_type in range(1,4):
        index = np.where(data['type']==p_type)
        #get total mass for all p_type
        totMass3 = sum(data['m'][index]*1e-2*u.M_sun)
        newMass3 = np.round(totMass3,3)
        #get total mass for desired p_type for M33
        if (p_type==ptype):
            totMassB = sum(data['m'][index]*1e-2*u.M_sun)
            newMassB = np.round(totMassB,3)
        #insert total mass of component into M33 array
        M33.append(newMass3)
        p_type+=1
    
    #total mass of components for each galaxy
    GalMassMW=sum(MW)
    MW.append(GalMassMW)
    print (MW)
    GalMassM31=sum(M31)
    M31.append(GalMassM31)
    print (M31)
    GalMassM33=sum(M33)
    M33.append(GalMassM33)
    print (M33)
    
    #fbar ratio calculation
    fbarMW = np.around((MW[1]+MW[2])/(MW[3]),3)
    print (fbarMW)
    fbarM31 = np.around((M31[1]+M31[2])/(M31[3]),3)
    print (fbarM31)
    fbarM33 = np.around((M33[1]+M33[2])/(M33[3]),3)
    print (fbarM33)
    
    print newMassH, newMassD, newMassB #print total mass of galaxy component
    
    return newMassH, newMassD, newMassD #return total mass of galaxy component

#if desired ptype=2
ptype = 1
ComponentMass(ptype)


# In[ ]:




