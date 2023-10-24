# -*- coding: utf-8 -*-
"""
Created on Tue Oct 24 13:12:08 2023

@author: josep
"""

##import packages

import numpy as np
import matplotlib.pyplot as plt
import time
from flux_function import flux
import matplotlib.pylab as plt
import seaborn as sns

start = time.time()

##constants

q = 1362
d = 5.394e2 #diffusivity
c = (5.25e9) #heat cap, initially constant, of dry planet
sb = 5.670374419e-8 #stefan-boltzmann constant
a = 0.3

n = 144
del_t = 1
del_lamb = 3.1415/n #in degrees->radians
initial_temp = 300
years = 40
final_time = 365*years #years * days

second_to_day = 24*60*60

##arrays+lists

temps = [initial_temp]*n #144 length list of initial temps
lats = np.linspace(-1.57, 1.57, n) #both start and end inclusive
solar = flux(lats, final_time) #returns a list of lists, indexable, returns list of flux over a year per lat
temp_diff = np.empty(n) # empty list
heat = [[0] for i in range(years)]
cumu = np.zeros(n)
av_list = []

##looping

for day in range(0, final_time, del_t): #gives list 0->364, 365 entries

    for y in range(len(lats)-1): #0->143, 144 entries
    
        #ir cooling function
        tir = 0.79*((temps[y]/273)**3)
        ir = ((sb*(temps[y])**4)/(1+0.75*tir))*second_to_day
        
        #albedo
        a = 0.525 - 0.245*np.tanh((temps[y]-268)/5)
        #print("albedo at ", lats[y], "is ", a)
        
        temp_diff[y] = (2*del_t/c)*(second_to_day*solar[y][day]*(1-a) - d*np.tan(lats[y])*((temps[y+1]-temps[y-1])/2*del_lamb) + d*((temps[y+1]-2*temps[y] + temps[y-1])/(del_lamb**2)) - ir)

    #ir cooling function
    tir = 0.79*((temps[0]/273)**3)
    ir = ((sb*(temps[0])**4)/(1+0.75*tir))*second_to_day
    
    #albedo
    a = 0.525 - 0.245*np.tanh((temps[0]-268)/5)
    #print("albedo at ", lats[y], "is ", a)
    temp_diff[0] = (2*del_t/c)*(second_to_day*solar[y][day]*(1-a) - d*((temps[y+1] - temps[y])/del_lamb**2) - ir) #south pole - what is day-1 for first day?
    
    #ir cooling function
    tir = 0.79*((temps[len(lats)-1]/273)**3)
    ir = ((sb*(temps[len(lats)-1])**4)/(1+0.75*tir))*second_to_day
    
    #albedo
    a = 0.525 - 0.245*np.tanh((temps[len(lats)-1]-268)/5)
    #print("albedo at ", lats[y], "is ", a)
    temp_diff[len(lats)-1] = (2*del_t/c)*(second_to_day*solar[y][day]*(1-a) - d*(temps[y] - temps[y-1])/(del_lamb**2) - ir) # north pole
    
    temps = temp_diff + temps
    #cumu = [temps[i] + cumu[i] for i in range(len(temps))]
    cumu = cumu + temps
    #print("temps on day ", day, temps)
    #print("cumulative at ", day, "is ", cumu)
    
    if day % 365 == 0:
        
        #heat.append(temps)
        year = int(day/365)
        #print("year is ", year)
        heat[year] = temps
        av_list.append([val/365 for val in cumu])
        cumu = np.zeros(n)
        print("year ", year, " length ", len(av_list))
        print("av ", av_list[year])
        
        fit = 302.3 - 45.3*(np.sin(lats)**2) #coakley model
        plt.figure(1)
        #plt.plot(lats, temps)
        plt.plot(lats, av_list[year])
        plt.plot(lats, n*[273.15]) #zero degrees celcius
        #plt.plot(lats, fit)
        plt.axis([-1.57, 1.57, 200, 350])
        plt.xlabel(f'latitude, $year = {day/365}$')
        plt.ylabel('Temps')
        plt.show()
        plt.pause(0.01) 

plt.figure(figsize=(10,10))
#y_axis_labels = np.linspace(-90, 90, 144)
heat_map = sns.heatmap(np.array(heat).T, linewidth = 0 , annot = False, cbar_kws={'label': 'Temperature (K)'}, yticklabels = 6)
plt.title(f"Temperatures per year as model settles over ${years}$ years, at each latitude")
plt.xlabel("Year")
plt.ylabel("Latitude")
plt.show()

end = time.time()
print("Time elapsed:", end - start)