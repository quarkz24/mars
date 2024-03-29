# -*- coding: utf-8 -*-
"""
Created on Wed Oct  4 16:18:25 2023

@author: josep
"""

##V4 adapted from V3, want to add in loop to iterate over 10 deg between 0 -> 90 and 0 -> -90
##my trig functions work in radians, as does np.sin etc, so latitude in rad (10deg = 0.17rad)

##NORTHERN HEMISPHERE
##PACKAGES AND CONSTANTS AND FUNCTIONS

##import packages

import numpy as np
import matplotlib.pyplot as plt
import time

start = time.time()

##constants

##earth
a = (149597870.7*1000)
q = 1362 #per second
#delta = -0.409*np.cos((2*3.1415/365)*(day + 10))
days = 365

##mars
a = (149597870.7*1000)*1.524
q = 589
#delta = -0.438*np.cos((2*3.1415/687)*(day)) #present day obliquity is 25.1deg (variable over hundreds of thousands of years)
days = 687



latlist = [0, 0.17, 0.35, 0.52, 0.70, 0.87, 1.05, 1.22, 1.40, 1.57]
for lat in latlist: #in radians

    sintheta = np.sin(lat)
    costheta = np.cos(lat)
    tantheta = np.tan(lat)
    
    flux = []

    for day in range(0, 687, 1):
    
        #2*3.1415/365 to convert day to position around sun
        #delta remains in rad to feed into my trig functions
        
        delta = -0.438*np.cos((2*3.1415/687)*(day))
        
        if -tantheta*np.tan(delta) >1:
            big_h = np.arccos(1)
        elif -tantheta*np.tan(delta) <-1:
            big_h = np.arccos(-1)
        else:
            big_h = np.arccos(-tantheta * np.tan(delta))
        
        solar = (q/3.1415) * (((149597870.7*1000)/a)**2) * (big_h*sintheta*np.sin(delta) + costheta*np.cos(delta)*np.sin(big_h))
        flux.append(solar)
        
    plt.plot(np.arange(0, 687, 1), flux, linestyle = '--', label = f'$latitude = {lat}$')

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

latlist2 = [-0, -0.17, -0.35, -0.52, -0.70, -0.87, -1.05, -1.22, -1.40, -1.57]

for lat in latlist2: #in radians


    sintheta = np.sin(lat)
    costheta = np.cos(lat)
    tantheta = np.tan(lat)
    
    flux = []

    for day in range(0, 687, 1):
    
        #2*3.1415/365 to convert day to position around sun
        #delta remains in rad to feed into my trig functions
            
        delta = -0.438*np.cos((2*3.1415/687)*(day))
        #big_h = np.arccos(-tantheta * np.tan(delta))
        
        if -tantheta*np.tan(delta) >1:
            big_h = np.arccos(1)
        elif -tantheta*np.tan(delta) <-1:
            big_h = np.arccos(-1)
        else:
            big_h = np.arccos(-tantheta * np.tan(delta))
            
        solar = (q/3.1415) * (((149597870.7*1000)/a)**2) * (big_h*sintheta*np.sin(delta) + costheta*np.cos(delta)*np.sin(big_h))
        flux.append(solar)
        
    plt.plot(np.arange(0, 687, 1), flux, linestyle = '--', label = f'$latitude = {lat}$')

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