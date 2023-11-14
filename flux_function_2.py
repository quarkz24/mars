# -*- coding: utf-8 -*-
"""
Created on Tue Nov 14 09:34:31 2023

@author: josep
"""

import numpy as np
import matplotlib.pyplot as plt
#import time

##constants (can get rid of these when its a functioning function)

#a = 149597870.7*1000
#q = 1360


#latlist = [0, 0.17, 0.35, 0.52, 0.70, 0.87, 1.05, 1.22, 1.40, 1.57]
latlist = np.linspace(-1.57, 1.57, 32)

def flux(latlist, years, year_length, semimajor, eccentricity, perihelion):
        
    days = years*year_length
    sol_longitude = np.linspace(0, years*2*np.pi, years*year_length)
    r = (semimajor*(1 - eccentricity**2))/(1 + eccentricity*np.cos(sol_longitude-perihelion)) #in AU, the 282 is earth specific! its perihelion Ls
    solar_constant = 1360
    q = solar_constant/r**2 
    q_list = q.tolist()
    
    all_flux = []
    
    for lat in latlist: #in radians
    
        sintheta = np.sin(lat)
        costheta = np.cos(lat)
        tantheta = np.tan(lat)
        
        flux = []
        for day in range(0, days, 1):
        
            #2*3.1415/365 to convert day to position around sun
            #delta remains in rad to feed into my trig functions
            
            delta = -0.409*np.cos((2*np.pi/365)*(day + 10))
            
            if -tantheta*np.tan(delta) >1:
                big_h = np.arccos(1)
            elif -tantheta*np.tan(delta) <-1:
                big_h = np.arccos(-1)
            else:
                big_h = np.arccos(-tantheta * np.tan(delta))
            
            a = 149597870.7*1000
            solar = (q_list[day]/np.pi) * (((149597870.7*1000)/a)**2) * (big_h*sintheta*np.sin(delta) + costheta*np.cos(delta)*np.sin(big_h))
            flux.append(solar)
            
        #plt.plot(np.arange(0, days, 1), flux, linestyle = '--', label = f'$latitude = {lat}$')
        #print("flux at", lat,":", flux)
        #print(len(flux))
        all_flux.append(flux)
    #print(len(all_flux))
    return(all_flux)

millionkm_to_AU = 0.00668458712
deg_to_rad = np.pi/180

#latlist, years, year_length, semimajor, eccentricity, perihelion
solartest = flux(latlist, 1, 365, 149.5978707*millionkm_to_AU, 0.017, 282.02*deg_to_rad)
print(solartest[0][364])