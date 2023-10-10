# -*- coding: utf-8 -*-
"""
Created on Tue Oct 10 14:04:17 2023

@author: josep
"""

##import packages

import numpy as np
import matplotlib.pyplot as plt
import time

start = time.time()

##constants

q = 1360
#d = #diffusivity
a = 149597870.7*1000 
#c = #heat cap
#solar = 
#albedo = 
del_t = 1
del_lamb = 2 * (3.1415/180) #in degrees->radians
initial_t = 300
#outgoing = ? #maybe take this out for now

##arrays/lists

latitude_list = np.linspace(0, 1.57, 45) #in radians
#times_list = #is this necessary? last program just used range(0, 365, 1)
temp_list = [300]*len(latitude_list) #this is the temperature at each of the 45 points on the Northern hemisphere

##for loops - finite difference method

#possibly add in: if first run through (if i = 0, first lat point, the equator) then taken the i-1 element to be 300 - 
#but this only works for first day - need to find way to take the -2 degrees point, which is also evolving
for day in range(0, 365, 1):
    for i in range(len(latitude_list)):
        new_temp = (2*del_t/c)*(solar*(1-albedo) + d*((temp_list[i+1] - 2*temp_list[i] + temp_list[i-1])/del_lamb) - d*np.tan(latitude_list[i])*((temp_list[i+1] - temp_list[i-1])/2*del_lamb)) + temp_list[i-1] - outgoing
        temp_list[i] = new_temp
        
##plotting
