#!/usr/bin/env python


import numpy as np
import astropy.units as u

def Read(filename):
    file = open(filename,'r')
    	 
    line1 = file.readline()                 #read the first line in MW_000.txt file
    label, value = line1.split()
    time = float(value)*u.Myr

    line2 = file.readline()                 #read the 2nd line in MW_000,txt file
    label, value = line2.split()
    total = float(value)

    file.close()

    data = np.genfromtxt(filename,dtype=None,names=True,skip_header=3) #store remainding data
    #print (data['type'][1])                 #print out particle type of 2nd particle
    
    return time, total, data

Read("MW_000.txt")

