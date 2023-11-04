# -*- coding: utf-8 -*-
"""
Created on Sat Nov  4 11:54:12 2023

@author: josep
"""

##import packages

import numpy as np
import matplotlib.pyplot as plt
import time
import math
from flux_function import flux

start = time.time()

def constants():
    d_t_s = 24*60*60
    q = 1362 #W m^-2 = J s^-1 m^-2
    d = (5.394e-1) #diffusivity, J m^-2 s^-1 K^-1
    cl = (5.25e6) #land heat cap, J m^-2 K^-1
    sb = 5.670374419e-8 #stefan-boltzmann, J s^-1 m^-2 K^-4
    n = 60
    del_t = 1
    del_secs = del_t*d_t_s
    del_lamb = np.pi/(n-1) #in degrees->radians
    initial_temp = 300
    years = 100
    days = 365*years
    return q, d, cl, sb, n, del_t, del_secs, del_lamb, initial_temp, years, days

q, d, cl, sb, n, del_t, del_secs, del_lamb, initial_temp, years, days = constants()

def list_of_lists():
    temp_and_time = [[0 for i in range(n)] for j in range(days)] #list of temps per latitude, per timestep
    temp_and_time[0] = [initial_temp for i in range (n)] #sets first temp at all lats to 300K
    lats = np.linspace(-np.pi/2, np.pi/2, n) #array of all n latitudes equally spaced
    solar = flux(lats, days)
    return temp_and_time, lats, solar

temp_and_time, lats, solar = list_of_lists()

def heatcap(capacity): #applying proportional lands to other temps is bad
    fraction = 1 - np.exp((temp_and_time[day-1][i] - 273)/10)
    if temp_and_time[day-1][i] >= 273:
        c = cl*40*0.7 + cl*0.3
    elif 263 < temp_and_time[day-1][i] < 273:
        c = fraction*9.2*cl + (1-fraction)*cl
    else:
        c = 0.7*2*cl + 0.3*cl #used to be 2
    return c

def outgoing_2(index):
    tir = 0.79*((temps[index]/273)**3)
    infrared = ((sb*(temps[index])**4)/(1+0.75*tir))
    albedo = 0.525 - 0.245*np.tanh((temps[index]-268)/5)
    return infrared, albedo

print(list_of_lists())

for latitude in range(len(lats)):
    print(lats[latitude])
    
print(lats[-1])

for day in range(1, days, del_t):
    for i in range(len(lats)): #n latitudes
        if i != 0 and i != n-1:
            tir = 0.79*((temp_and_time[day-1][i]/273)**3)
            ir = ((sb*(temp_and_time[day-1][i])**4)/(1+0.75*tir))
            a = 0.525 - 0.245*np.tanh((temp_and_time[day-1][i]-268)/5)
            
            first_order = -d*np.tan(lats[i])*(temp_and_time[day-1][i+1] - temp_and_time[day-1][i-1])/(2*del_lamb)
            second_order = d*(temp_and_time[day-1][i+1] - 2*temp_and_time[day-1][i] + temp_and_time[day-1][i-1])/(del_lamb**2)
            rad_terms = (solar[i][day-1])*(1-a) - ir
            
            change = (first_order + second_order + rad_terms)*(del_secs/heatcap(cl))
        if i == 0: #negative 90deg
            tir = 0.79*((temp_and_time[day-1][i]/273)**3)
            ir = ((sb*(temp_and_time[day-1][i])**4)/(1+0.75*tir))
            a = 0.525 - 0.245*np.tanh((temp_and_time[day-1][i]-268)/5)
            
            second_order = d*(temp_and_time[day-1][0+1] - temp_and_time[day-1][0])/(del_lamb**2)
            rad_terms = (solar[0][day-1])*(1-a) - ir
            
            change = (second_order + rad_terms)*(del_secs/heatcap(cl))
        if i == n-1: #positive 90deg
            tir = 0.79*((temp_and_time[day-1][-1]/273)**3)
            ir = ((sb*(temp_and_time[day-1][-1])**4)/(1+0.75*tir))
            a = 0.525 - 0.245*np.tanh((temp_and_time[day-1][-1]-268)/5)
            
            second_order = -d*(temp_and_time[day-1][-1] - temp_and_time[day-1][-1-1])/(del_lamb**2)
            rad_terms = (solar[-1][day-1])*(1-a) - ir
            
            change = (second_order + rad_terms)*(del_secs/heatcap(cl))
        
        temp_and_time[day][i] = temp_and_time[day-1][i] + change
        
    if day % 365 == 0:
        plt.plot(lats, temp_and_time[day])
        plt.plot(lats, n*[273.15], color = 'black') #zero degrees celcius
        plt.plot(lats, 302.3 - 45.3*(np.sin(lats)**2), label = "N&C fit", color = 'red') #coakley model
        plt.xlabel(f'latitude, $year = {day/365}$')
        plt.ylabel('Temps')
        plt.show()
        plt.pause(0.5)

end = time.time()
print("Time elapsed:", end - start)