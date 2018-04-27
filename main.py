# -*- coding: utf-8 -*-
"""
Created on Thu Dec  7 18:09:29 2017

@author: majid
"""
import numpy as np

from timestep import timestep
from params import *


#=================================
pi = 3.141592653
G = 6.673e-8
sigma_SB = 5.6704e-5
D = 5.394e2
freeze = 273.0	
Msol = 1.99e33
Lsol = 3.826e33
AU = 1.496e13
earth_yr = 3.1556926e7 
q0 = Lsol/(4.0*pi*AU*AU)	
p0 = 1 
C_atm0 = 10.1e9
    
num_point = 144
star_mass	= 1.0
Lstar = 1.0         
semi_major_axis = 1.0		
e = 0.0	
phi_peri = 0.0         
phi = 0.0		
rot_period = 1		
spin_obliq = 0.0		
azim_obliq = 0.0       
f_ocean = 0.7					
maxtime = 1		
dumpfreq = 0.01		
p = 1
time = 0.0
deltat = 0.0 
To = 288
T = np.ones(num_point+1)*To
#=========================================
latt, latdeg, delta_latt = mesh_1D(num_point)


# some variables are not the same
phi_peri,spin_obliq,azim_obliq,phi,phidot,r = orbit_single_star(semi_major_axis, star_mass,earth_yr,e ,phi, phi_peri,spin_obliq,azim_obliq)

#correct
insol = INSOLATION(num_point,q0,Lstar,r,latt, spin_obliq, phi, phi_peri, azim_obliq)

# correct!
albedo = ALBEDO(num_point,spin_obliq,phi,phi_peri,azim_obliq,latt,f_ocean,T)

# correct
diff = DIFFUSION(D, rot_period, p, p0)

#correct
boil = water_boil(p)

# correct
C_to = C_tot(num_point,p, p0, C_atm0,f_ocean,freeze,T)

#correct
tau = tau_ir(p, p0, T, freeze)

#tau! BTW correct
infrared = OLR(sigma_SB,T,tau)

# correct
qq = Q(insol, albedo, infrared)

#
ha = hab(num_point, T,freeze, boil)


deltat = timestep(diff, latt, delta_latt, num_point, C_to )
"""

phi += phidot*deltat
diff, deltax, x,  insol[i], f_ice[i],  C_tot[i],  albedo[i], tau_ir[i],infrared[i],Q[i],hab[i] = value_update(T,phi,deltat)

deltat = timestep(diff, x, deltax, nx, C_tot)
phi += phidot*deltat
diff, deltax, x,  insol[i], f_ice[i],  C_tot[i],  albedo[i], tau_ir[i],infrared[i],Q[i],hab[i] = value_update(T,phi,deltat)

T_max = freeze + 0.9*(boil - freeze)

timeyr = time/yr

T_mean_old = np.mean(T)

#mity = time.clock()
while timeyr < maxtime :
    
    diff, deltax, x,  insol[i], f_ice[i],  C_tot[i],  albedo[i], tau_ir[i],infrared[i],Q[i],hab[i] = value_update(T,phi,deltat)
    time =time+deltat
    timeyr = time/yr
    







timeyr = time/yr
"""
