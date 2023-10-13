# -*- coding: utf-8 -*-
"""
Created on Fri Oct 13 09:05:39 2023

@author: josep
"""

##import packages

import numpy as np
import matplotlib.pyplot as plt
import time
from flux_function import flux

start = time.time()

##constants

q = 1360
d = 0.0001 #diffusivity, placeholder, inaccurate
a = 0.5 #also false
c = 40 * 5.25e9 #heat cap, initially constant, of water planet 
sb = 5.670374419e-8

n = 144
del_t = 1
del_lamb = 1.25 * (3.1415/180) #in degrees->radians
initial_temp = 300
final_time = 365

##arrays+lists

temps = [initial_temp]*n #144 length list of initial temps
lats = np.linspace(-1.57, 1.57, n) #both start and end inclusive
solar = flux(lats) #returns a list of lists, indexable, returns list of flux over a year per lat

print("lats: ",len(lats))
print("temps: ", len(temps))
print("incoming flux: ", len(solar)) #144 lists of flux, for 144 latitudes


##looping
##only differential terms included so far - without outgoing or incoming flux, think this should have no change?

for day in range(0, final_time, del_t): #gives list 0->364
    for y in range(len(lats)-1): #0->143
        temps[y+1] = (2*del_t/c)*(-d*np.tan(lats[y]))*((temps[y+1]-temps[y-1])/2*del_lamb) + d*((temps[y+1]-2*temps[y]+temps[y-1])/del_lamb**2) + temps[day-1]
                                   
    temps[len(lats)-1] = (2*del_t/c)*(-d*(temps[y] - temps[y-1])/(del_lamb**2)) + temps[day-1] # north pole
    temps[0] = (2*del_t/c)*(d*((temps[y+1] - temps[y])/del_lamb**2)) + temps[day-1] #south pole - what is day-1 for first day?
    ##maybe need to put time addition in time loop
    
    print("temps for day", day, ":", temps)
    
    plt.figure(1)
    plt.plot(lats, temps)
    plt.axis([-1.57, 1.57, 0, 700])
    plt.xlabel('Distance (m), day')
    plt.ylabel('Temperature (C)')
    plt.show()
    plt.pause(0.01) 
    
end = time.time()
print("Time elapsed:", end - start)