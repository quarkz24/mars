# -*- coding: utf-8 -*-
"""
Created on Wed Oct 11 14:16:24 2023

@author: josep
"""

##import packages

import numpy as np
import matplotlib.pyplot as plt
import time
from flux_function import flux

start = time.time()

##constants
#diffusion NOT accurate

q = 1360
d = 0.0001 #diffusivity, from 'habitable climates'
a = 149597870.7*1000 
c = 40 * 5.25e9 #heat cap, initially constant, of water planet 
sb = 5.670374419e-8

n = 144
del_t = 1
del_lamb = 1.25 * (3.1415/180) #in degrees->radians
initial_t = 300
final_time = 365

print

##arrays/lists
#180/1.25 = 144 nodes, so 143 midpoints?
latitude_list = np.linspace(-1.57 + del_lamb/2, 1.57-del_lamb/2,n) #goes from -1.559 to 1.559, 144 points, in 1.25deg interval
print("latitudes: ", latitude_list, "length of latitude list: ", len(latitude_list))
temp_list = [initial_t]*len(latitude_list) # 144 long, all 300
print("initial temps: ", temp_list, "length of initial temps list: ", len(temp_list))
#dTdt = np.zeros(n) # list of 144 zeroes
#print("deriv. list: ", dTdt, "length: ", len(dTdt))
times = np.arange(0, final_time, 1) # 0->364, 365 long
print("times: ", times, "length: ", len(times))
temp_list = [initial_t]*len(latitude_list) #144 300s
print("temperatures list: ", temp_list, "length: ", len(temp_list))

#s(1-a), I, not included
for day in range(1, len(times)): #1->364 - only 364 loops?
    
    for y in range(0, n-1): #1->142, since n-1 is 143 - surely need 143? #note i have changed this - guy in vid starts at 1
        
        tir = 0.79*((temp_list[y]/273)**3)
        outgoing = (sb * temp_list[y]**4)/(1+ (3/4)*tir)
        
        temp_list[y] = (2*del_t/c)*(d*((temp_list[y+1] - 2*temp_list[y] + temp_list[y-1])/del_lamb**2) - d*np.tan(latitude_list[y])*((temp_list[y+1] - temp_list[y-1])/2*del_lamb) - outgoing) + temp_list[y-1] #wondering whether del_lamb rly is 1.25, given the shortened latitude list
    
    temp_list[0] = (2*del_t/c)*(d*(temp_list[0+1] - temp_list[0])/del_lamb**2 - outgoing) + temp_list[0-1] #for south pole - what do with temp_list[-1]? boundary condition?
    temp_list[n-1] = (2*del_t/c)*(d*(-(temp_list[n-1] - temp_list[n-2])/del_lamb**2) - outgoing) + temp_list[n-2] #for north pole

    plt.figure(1)
    plt.plot(latitude_list, temp_list)
    plt.axis([-1.57, 1.57, 0, 700])
    plt.xlabel('Distance (m), day')
    plt.ylabel('Temperature (C)')
    plt.show()
    plt.pause(0.01) 
    

end = time.time()
print("Time elapsed:", end - start)
