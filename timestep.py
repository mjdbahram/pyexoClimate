# -*- coding: utf-8 -*-
"""
Created on Thu Dec  7 18:15:36 2017

@author: majid
"""
import numpy as np



def timestep(diff, latt, delta_latt, num_point, C_tot ):
    x      = np.sin(latt)
    delta_x = np.cos(latt)*delta_latt
    t_try = np.zeros(num_point+1)
    for i in range(num_point+1):
        if i == num_point :
            Dplus = diff*(1.0 - x[i]*x[i])
        if i!=num_point :
            Dplus = 0.5*diff*((1.0-x[i+1]*x[i+1]) + (1.0-x[i]*x[i]))
            
        if i!= 0 and i != num_point :
            t_try[i] = (delta_x[i]*delta_x[i]*C_tot[i])/(2.0*Dplus)
        else :
            t_try[i] = 1.0e30
            
    deltat = 0.9*np.min(t_try)
    return deltat
    
  



