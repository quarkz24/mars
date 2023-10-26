# -*- coding: utf-8 -*-
"""
Created on Wed Oct 25 16:03:52 2023

@author: josep
"""

##import packages

import numpy as np
import matplotlib.pyplot as plt
import time
from flux_function import flux
import seaborn as sns

start = time.time()

##constants

second_to_day = 24*60*60
q = 1362 #bolometric solar flux at 1.0AU
d = (5.394e-1)*second_to_day #diffusivity
#cl = (5.25e6) #heat capacity of land
cl = 5.25e7
sb = 5.670374419e-8 #stefan-boltzmann constant
n = 144 #no. of steps
del_t = 1 #timestep in days
del_lamb = 3.1415/n #in degrees->radians
initial_temp = 300 #initial uniform temperature
years = 200
final_time = 365*years #years * days

##arrays+lists

temps = [initial_temp]*n #144 length list of initial temps
lats = np.linspace(-1.57, 1.57, n) #both start and end inclusive
solar = flux(lats, final_time) #returns a list of lists, indexable, returns list of flux over a year per lat
temp_diff = np.empty(n) # empty list
heat = [[0] for i in range(years)]
cumu = np.zeros(n)
av_list = []

##functions

def heatcap(capacity):
    if temps[y] >= 273:
        c = cl*40*0.7 + cl*0.3
    elif 263 < temps[y] < 273:
        c = 9.2*cl
    else:
        c = 2*cl
    return c

def temp_difference(index):
    if index != 0 and index != -1: #most latitudes
        tir = 0.79*((temps[index]/273)**3)  
        ir = ((sb*(temps[index])**4)/(1+0.75*tir))*second_to_day
        a = 0.525 - 0.245*np.tanh((temps[index]-268)/5)
        temp_diff[index] = (del_t/heatcap(cl))*(second_to_day*solar[y][day]*(1-a) - d*np.tan(lats[index])*((temps[index+1]-temps[index-1])/2*del_lamb) + d*((temps[index+1]-2*temps[index] + temps[index-1])/(del_lamb**2)) - ir)
        
    elif index == 0: #south pole
        tir = 0.79*((temps[0]/273)**3)
        ir = ((sb*(temps[0])**4)/(1+0.75*tir))*second_to_day
        a = 0.525 - 0.245*np.tanh((temps[0]-268)/5)
        temp_diff[0] = (del_t/heatcap(cl))*(second_to_day*solar[0][day]*(1-a) + d*((temps[0+1] - temps[0])/del_lamb**2) - ir)
        
    else: #north pole
        tir = 0.79*((temps[-1]/273)**3)
        ir = ((sb*(temps[-1])**4)/(1+0.75*tir))*second_to_day
        a = 0.525 - 0.245*np.tanh((temps[-1]-268)/5)
        temp_diff[-1] = (del_t/heatcap(cl))*(second_to_day*solar[-1][day]*(1-a) - d*(temps[-1] - temps[-1-1])/(del_lamb**2) - ir)
        
    return temp_diff
        
##looping

for day in range(0, final_time, del_t): #gives list 0->364, 365 entries
    for y in range(len(lats)-1): #0->143, 144 entries
        
        change = temp_difference(y)
    
    temps = change + temps
    cumu = cumu + temps
    
    if day % 365 == 0:
        
        year = int(day/365)
        heat[year] = temps
        av_list.append([val/365 for val in cumu])
        cumu = np.zeros(n) #clear cumulative list for next year
        plt.figure(1)
        if year > 9:
            add = [sum(i) for i in zip(av_list[-1], av_list[-2], av_list[-3], av_list[-4], av_list[-5], av_list[-6], av_list[-7], av_list[-8], av_list[-9], av_list[-10])]
            div = [val/10 for val in add]
            plt.plot(lats, div, label = "10 year average", color = 'blue', linestyle = '--') #ten yr average
        
        fit = 302.3 - 45.3*(np.sin(lats)**2) #coakley model
        plt.plot(lats, temps, label = "finite diff method", color = 'blue') #temps per year
        plt.plot(lats, n*[273.15], color = 'black') #zero degrees celcius
        plt.plot(lats, fit, label = "N&C fit", color = 'red') #coakley model
        plt.axis([-1.57, 1.57, 200, 350])
        plt.xlabel(f'latitude, $year = {day/365}$')
        plt.ylabel('Temps')
        plt.legend()
        plt.show()
        plt.pause(0.01) 

plt.figure(figsize=(10,10))
heat_map = sns.heatmap(np.array(heat).T, linewidth = 0 , annot = False, cbar_kws={'label': 'Temperature (K)'}, yticklabels = 6)
plt.title(f"Temperatures per year as model settles over ${years}$ years, at each latitude")
plt.xlabel("Year")
plt.ylabel("Latitude")
plt.show()

end = time.time()
print("Time elapsed:", end - start)