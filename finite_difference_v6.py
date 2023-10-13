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
a = 149597870.7*1000 
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

##looping
##only differential terms included so far

for day in range(0, final_time, del_t): #gives list 0->364
    for y in range(0, len(lats)): #0->143
        temps[y+1] = (2*del_t/c) *(-d*np.tan(lats[y])*((temps[y+1]-temps[y-1])/2*del_lamb) + d*((temps[y+1]-2*temps[y]+temps[y-1])/del_lamb**2) + temps[day-1]
    temps[0] = (2*del_t/c)*(d*((temps[y+1]-temps[y])/del_lamb**2)) + temps[day-1] #south pole - what is day-1 for first day?
    temps[len(lats)] = (2*del_t/c)*(-d*(temps[y] - temps[y-1])/(del_lamb**2)) + temps[day-1] # north pole
    ##maybe need to put time addition in time loop