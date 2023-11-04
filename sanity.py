# -*- coding: utf-8 -*-
"""
Created on Wed Nov  1 12:00:10 2023

@author: josep
"""
"""
import numpy as np
import matplotlib.pyplot as plt

temps = np.linspace(240, 350, 110)
print(temps)
y = np.zeros(len(temps))
print(y)
cl = 5.25e6

def heatcap(capacity):
    if temps[i] >= 273:
        c = cl*40*0.7 + cl*0.3
    elif 263 < temps[i] < 273:
        c = 9.2*cl
    else:
        c = 2*cl #used to be 2
    return c

for i in range(len(temps)):
    val = heatcap(temps[i]) 
    y[i] = val
    print(y)

print(y)
plt.plot(temps, y)
"""
"""
##general
second_order_term = d*(temps[index+1] - 2*temps[index] + temps[index-1])/(del_lamb**2)
first_order_term = -d*np.tan(temps[index])*((temps[index+1] - temps[index-1])/(2*del_lamb))
rad_term = solar[index][day]*(1-a) - ir

temp_diff[index] = (second_order_term + first_order_term + rad_term)*(del_t/heatcap(cl))

##south pole, ie the first index, ie -90deg, ie index = 0
second_order_term = d*((temps[index+1] - temps[index])/del_lamb**2)
rad_term = solar[index][day]*(1-a) - ir

temp_diff[index] = (second_order_term + rad_term)*(del_t/heatcap(cl))

##north pole, ie last index, ie +90deg, ie n-1 term
second_order_term = -d*((temps[index] - temps[index-1])/del_lamb**2)
rad_term = solar[index][day]*(1-a) - ir

temp_diff[index] = (second_order_term + rad_term)*(del_t/heatcap(cl))
"""

