# -*- coding: utf-8 -*-
"""
Created on Wed Oct  4 09:20:43 2023

@author: josep
"""

##V3 adapted from V1, using np.cos functions in delta (also could try increasing order of Taylor expansion)

##PACKAGES AND CONSTANTS AND FUNCTIONS

##import packages

import numpy as np
import matplotlib.pyplot as plt

##define constants (or things that are constant for now)
##my trig functions work in radians, as does np.sin etc, so latitude in rad (10deg = 0.17rad)

latitude = 0.1 #in radians
a = 149597870.7*1000

sintheta = np.sin(latitude)
costheta = np.cos(latitude)
tantheta = np.tan(latitude)

q = 1360

##taylor functions for efficiency - use radians

def tan(x):
	return x + (x**3)/6 + (2*x**5)/15
	
def sin(x):
	return x - (x**3)/6 + (x**5)/120
	
def cos(x):
	return 1 - (x**2)/2 + (x**4)/24

##LOOPS

##to be in a loop over values of declination (del) or time (t)
##day 0 to day 365, step = 1 day

flux = []
d = []
h = []

for day in range(0, 365, 1):
    
    #2*3.1415/365 to convert day to position around sun
    #-0.409rad is peak delta, remains in rad to feed into my trig functions
    
    delta = -0.409*np.cos((2*3.1415/365)*(day + 10)) # why is delta returning as anything larger than 0.409?
    big_h = np.arccos(-tantheta * np.tan(delta))
    solar = (q/3.1415) * (((149597870.7*1000)/a)**2) * (big_h*sintheta*np.sin(delta) + costheta*np.cos(delta)*np.sin(big_h))
    flux.append(solar)
    d.append(delta)
    h.append(big_h)

print("delta list is:", d)
print("big H list is:", h)
#print("flux:", flux)

##PLOTTING

plt.plot(np.arange(0, 365, 1), flux, linestyle = '--')

plt.ylabel("Diurnally averaged incident solar flux", fontsize = 10)
plt.xlabel("Days of the year")
plt.title("Diurnally averaged incident solar flux per day of the year, beginning Jan 1")

#plt.plot(np.arange(0, 365, 1), d)