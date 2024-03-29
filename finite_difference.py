# -*- coding: utf-8 -*-
"""
Created on Tue Oct 10 14:04:17 2023

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
d = 5.394*100 #diffusivity, from 'habitable climates'
a = 149597870.7*1000 
c = 40 * 5.25e9 #heat cap, initially constant, of water planet 
del_t = 1
del_lamb = 2 * (3.1415/180) #in degrees->radians
initial_t = 350
sb = 5.670374419e-8

##arrays/lists

latitude_list = np.linspace(-1.57, 1.57, 90) #in radians
temp_list = [initial_t]*len(latitude_list) #this is the temperature at each of the 45 points on the Northern hemisphere

#print("temps is:", temp_list)
##for loops - finite difference method

#possibly add in: if first run through (if i = 0, first lat point, the equator) then taken the i-1 element to be 300 - 
#but this only works for first day - need to find way to take the -2 degrees point, which is also evolving
#maybe within second for loop, before new_temp executes, need line to work out I and A (both temp dependent)
#think i can define flux outside of loop? since it isn't temp dependent, only latitude
#need to index into solar list

solar = flux(latitude_list)

for day in range(0, 365, 1):
    for i in range(len(latitude_list)-1): #does i start at zero?
        
        #IR cooling function
        tir = 0.79*((temp_list[i]/273)**3)
        outgoing = (sb * temp_list[i]**4)/(1+ (3/4)*tir)
        
        #albedo
        albedo = 0.525 - 0.245*np.tanh((temp_list[i]-268)/5)
        
        #finite difference method
        new_temp = (2*del_t/c)*(solar[i][day]*(1-albedo) + d*((temp_list[i+1] - 2*temp_list[i] + temp_list[i-1])/del_lamb) - d*np.tan(latitude_list[i])*((temp_list[i+1] - temp_list[i-1])/2*del_lamb)) + temp_list[day-1] - outgoing
        temp_list[i] = new_temp
   
    print("temps for day", day, ":", temp_list)
    
    #plotting
    plt.figure(1)
    plt.plot(latitude_list, temp_list)
    plt.axis([-1.57, 1.57, 0, 700])
    plt.xlabel('Distance (m), day ', )
    plt.ylabel('Temperature (C)')
    plt.show()
    plt.pause(0.01) 
#print(solar)
print()
        


end = time.time()
print("Time elapsed:", end - start)
