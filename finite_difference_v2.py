# -*- coding: utf-8 -*-
"""
Created on Wed Oct 11 11:20:56 2023

@author: josep
"""

import numpy as np
import matplotlib.pyplot as plt
import time
from flux_function import flux

start = time.time()

##constants

l = 3.1415
lat = 1.5708
n = 144 #l/n = 1.25 degrees sep.

q = 1360
d = 5.394*100 #diffusivity, from 'habitable climates'
a = 149597870.7*1000 
c = 40 * 5.25e9 #heat cap, initially constant, of water planet 
sb = 5.670374419e-8

del_t = 1
del_lamb = l/n #in degrees->radians
initial_temp = 300 #don't need different temps for poles
final_time = 365 #days, 1yr for now

latitude_list = np.linspace(-1.57, 1.57, 90) # this prints a fine list, but maybe need to be in middle of columns
temp_list = [initial_temp]*len(latitude_list)
time_list = np.arange(0, final_time)
#print(latitude_list, len(latitude_list)) #this works
#print(temp_list, len(temp_list)) #this works

dTdt = np.empty(n)
solar = flux(latitude_list) #despite runtime warning, this does work, creates zero values as it should
#print(type(solar[0]))

for day in range(1, len(time_list)):
    for i in range(1, n-1):
        
        #IR cooling function
        tir = 0.79*((temp_list[i]/273)**3)
        outgoing = (sb * temp_list[i]**4)/(1+ (3/4)*tir)
        #print(type(outgoing))
        #albedo
        albedo = 0.525 - 0.245*np.tanh((temp_list[i]-268)/5)
        #print(type(albedo))
        dTdt[i] = solar[i][day]*(1-albedo) + d*((temp_list[i+1] - 2*temp_list[i] + temp_list[i-1])/del_lamb**2) - d*np.tan(latitude_list[i])*((temp_list[i+1]-temp_list[i-1])/2*del_lamb) - outgoing
    
    dTdt[1.57] = solar*(1-albedo) - ((temp_list[i] - temp_list[i-1])/del_lamb**2) - outgoing #dT/dlamb is zero, second deriv is back diff
    dTdt[-1.57] = solar*(1-albedo) + ((temp_list[i+1] - temp_list[i])/del_lamb**2) - outgoing #dT/dlamb is zero, second deriv is forward diff
    
    temp_list = temp_list + dTdt*del_t
    
#remember to divide through by C
end = time.time()
print("Time elapsed:", end - start)


