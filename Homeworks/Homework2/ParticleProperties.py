import numpy as np
import astropy.units as u
from ReadFile import Read     #to read the file again

def ParticleInfo(readFile,p_type,p_total):

    time, total, data = Read("MW_000.txt")

    index = np.where(data['type']==p_type)

    mag_dist = np.sqrt((data['x'][index][p_total])**2 + (data['y'][index][p_total])**2 + (data['z'][index][p_total])**2)*u.kpc
    dist = np.around(mag_dist,3)
    print(dist)

    
    mag_vel = np.sqrt((data['vx'][index][p_total])**2 + (data['vy'][index][p_total])**2 + (data['vz'][index][p_total])**2)*u.km/u.s
    vel = np.around(mag_vel,3)
    print(vel)

    
    mass = data['m'][index][p_total]*u.M_sun
    print(mass)

    dist_new = mag_dist.to(u.lyr)
    dist_new = np.around(dist_new,3)
    print(dist_new)

ParticleInfo("MW_000.txt",2.00000,99)
