# -*- coding: utf-8 -*-
"""
Created on Wed Oct  4 16:18:25 2023

@author: josep
"""

##V4 adapted from V3, want to add in loop to iterate over 10 deg between 0 -> 90 and 0 -> -90
##my trig functions work in radians, as does np.sin etc, so latitude in rad (10deg = 0.17rad)
##still getting nan issues
##problem is tan(delta) returning large values (as it mathematically should)

##NORTHERN HEMISPHERE
##PACKAGES AND CONSTANTS AND FUNCTIONS

##import packages

import numpy as np
import matplotlib.pyplot as plt
import time

start = time.time()

##define constants (or things that are constant for now)
##my trig functions work in radians, as does np.sin etc, so latitude in rad (10deg = 0.17rad)

a = 149597870.7*1000
q = 1360

latlist = [0, 0.17, 0.35, 0.52, 0.70, 0.87, 1.05, 1.22, 1.40, 1.56]
for lat in latlist: #in radians

    sintheta = np.sin(lat)
    costheta = np.cos(lat)
    tantheta = np.tan(lat)
    
    flux = []

    for day in range(0, 365, 1):
    
        #2*3.1415/365 to convert day to position around sun
        #delta remains in rad to feed into my trig functions
        if tantheta >1:
            tantheta = 1
        elif tantheta <-1:
            tantheta = -1
        
        delta = -0.409*np.cos((2*3.1415/365)*(day + 10))
        big_h = np.arccos(-tantheta * np.tan(delta))
        solar = (q/3.1415) * (((149597870.7*1000)/a)**2) * (big_h*sintheta*np.sin(delta) + costheta*np.cos(delta)*np.sin(big_h))
        flux.append(solar)
        
    plt.plot(np.arange(0, 365, 1), flux, linestyle = '--', label = f'$latitude = {lat}$')

##TESTING
#print("delta list is:", d)
#print("big H list is:", h)
#print("flux:", flux)

plt.ylabel("Diurnally averaged incident solar flux", fontsize = 10)
plt.xlabel("Days of the year")
plt.title("Diurnally averaged incident solar flux per day of the year, beginning Jan 1 (Northern hemisphere)")

##LEGEND
ax = plt.subplot(111)
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.show()

##SOUTHERN HEMISPHERE

latlist2 = [-0.017, -0.17, -0.35, -0.52, -0.70, -0.87, -1.05, -1.22, -1.40, -1.56]

for lat in latlist2: #in radians


    sintheta = np.sin(lat)
    costheta = np.cos(lat)
    tantheta = np.tan(lat)
    
    flux = []

    for day in range(0, 365, 1):
    
        #2*3.1415/365 to convert day to position around sun
        #delta remains in rad to feed into my trig functions
        if tantheta >1:
            tantheta = 1
        elif tantheta <-1:
            tantheta = -1
            
        delta = -0.409*np.cos((2*3.1415/365)*(day + 10))
        big_h = np.arccos(-tantheta * np.tan(delta))
        solar = (q/3.1415) * (((149597870.7*1000)/a)**2) * (big_h*sintheta*np.sin(delta) + costheta*np.cos(delta)*np.sin(big_h))
        flux.append(solar)
        
    plt.plot(np.arange(0, 365, 1), flux, linestyle = '--', label = f'$latitude = {lat}$')

##TESTING
#print("delta list is:", d)
#print("big H list is:", h)
#print("flux:", flux)

plt.ylabel("Diurnally averaged incident solar flux", fontsize = 10)
plt.xlabel("Days of the year")
plt.title("Diurnally averaged incident solar flux per day of the year, beginning Jan 1 (Southern hemisphere)")

##LEGEND
ax = plt.subplot(111)
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.show()

end = time.time()
print("Time elapsed:", end - start)