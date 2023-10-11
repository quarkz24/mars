# -*- coding: utf-8 -*-
"""
Created on Tue Oct 10 15:35:41 2023

@author: josep
"""
##FLUX FUNCTION
##intended as a program that I can pull into finite_difference.py
##so no constants/plotting found here
##returns a list of lists, indexable

##import packages

import numpy as np
import matplotlib.pyplot as plt
import time

##constants (can get rid of these when its a functioning function)

a = 149597870.7*1000
q = 1360

#latlist = [0, 0.17, 0.35, 0.52, 0.70, 0.87, 1.05, 1.22, 1.40, 1.57]

def flux(latlist):
        
    all_flux = []
    
    for lat in latlist: #in radians
    
        sintheta = np.sin(lat)
        costheta = np.cos(lat)
        tantheta = np.tan(lat)
        
        flux = []
    
        for day in range(0, 365, 1):
        
            #2*3.1415/365 to convert day to position around sun
            #delta remains in rad to feed into my trig functions
            
            delta = -0.409*np.cos((2*3.1415/365)*(day + 10))
            
            if -tantheta*np.tan(delta) >1:
                big_h = np.arccos(1)
            elif -tantheta*np.tan(delta) <-1:
                big_h = np.arccos(-1)
            else:
                big_h = np.arccos(-tantheta * np.tan(delta))
            
            solar = (q/3.1415) * (((149597870.7*1000)/a)**2) * (big_h*sintheta*np.sin(delta) + costheta*np.cos(delta)*np.sin(big_h))
            flux.append(solar)
            
        #plt.plot(np.arange(0, 365, 1), flux, linestyle = '--', label = f'$latitude = {lat}$')
        #print("flux at", lat,":", flux)
        all_flux.append(flux)
        
    return(all_flux)

#print(flux(latlist))
