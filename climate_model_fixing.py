# -*- coding: utf-8 -*-
"""
Created on Mon Oct 23 09:13:52 2023

@author: josep
"""

##finite difference model including ir cooling function, albedo, and solar incident radiation
##import packages

import numpy as np
import matplotlib.pyplot as plt
import time
from flux_function import flux
import matplotlib.pylab as plt
import seaborn as sns

start = time.time()

##constants

q = 1360
d = 5.394e2 #diffusivity
c = (5.25e9) #heat cap, initially constant, of dry planet
sb = 5.670374419e-8 #stefan-boltzmann constant

n = 144
del_t = 1
del_lamb = 3.1415/n #in degrees->radians
initial_temp = 300
years = 2
final_time = 365*years #years * days

second_to_day = 24*60*60

##arrays+lists

temps = [initial_temp]*n #144 length list of initial temps
lats = np.linspace(-1.57, 1.57, n) #both start and end inclusive
solar = flux(lats, final_time) #returns a list of lists, indexable, returns list of flux over a year per lat
temp_diff = np.empty(n) # empty list
heat = [[0] for i in range(years)]
#little_list = np.empty(n) #not sure if range(temps) will work here
little_list = [0]*n
big_list = [[0] for i in range(years)]

print("little list: ", little_list)
print("big list: ", big_list)

##looping

for day in range(0, final_time, del_t): #gives list 0->364, 365 entries
    print("day is: ", day)
    for y in range(len(lats)-1): #0->143, 144 entries
    
        #ir cooling function
        tir = 0.79*((temps[y]/273)**3)
        ir = ((sb*(temps[y])**4)/(1+0.75*tir))*second_to_day
        
        #albedo
        a = 0.525 - 0.245*np.tanh((temps[y]-268)/5)
        
        temp_diff[y] = (del_t/c)*(second_to_day*solar[y][day]*(1-a) - d*np.tan(lats[y])*((temps[y+1]-temps[y-1])/2*del_lamb) + d*((temps[y+1]-2*temps[y] + temps[y-1])/(del_lamb**2)) - ir)
        #print("lat: ", y, ", tempdiff: ", temp_diff[y])
        
    #ir, albedo for south
    tir = 0.79*((temps[0]/273)**3)
    ir = ((sb*(temps[0])**4)/(1+0.75*tir))*second_to_day
    a = 0.525 - 0.245*np.tanh((temps[0]-268)/5)
    
    #changed it to +d, and y=0
    temp_diff[0] = (del_t/c)*(second_to_day*solar[0][day]*(1-a) + d*((temps[0+1] - temps[0])/del_lamb**2) - ir) #south pole - what is day-1 for first day?
    #print("lat: south pole", 0, ", tempdiff: ", temp_diff[0])
    
    #ir, albedo for north
    tir = 0.79*((temps[len(lats)-1]/273)**3)
    ir = ((sb*(temps[len(lats)-1])**4)/(1+0.75*tir))*second_to_day
    a = 0.525 - 0.245*np.tanh((temps[len(lats)-1]-268)/5)
    
    temp_diff[len(lats)-1] = (del_t/c)*(second_to_day*solar[len(lats)-1][day]*(1-a) - d*(temps[len(lats)-1] - temps[len(lats)-1-1])/(del_lamb**2) - ir) # north pole
    #print("lat: north pole", len(lats)-1, ", tempdiff: ", temp_diff[len(lats)-1])
    
    temps = temp_diff + temps
    little_list = little_list + temps
    
    if day % 365 == 0:
        
        print("on day ", day, "cumulative list is ", little_list)
        plt.figure(1)
        #heat.append(temps)
        year = int(day/365)
        big_list[year] = [val/365 for val in little_list] #takes average per lat for the year, put it in first list of the 2d list
        print("list of lists indexed at year: ", year, "is ", big_list[year])
        #add up all lists from element -1 to element -10
        #FIVE YEAR AVERAGE WHILE I TEST IT
        if day > 365*6:
            average_total = [big_list[-1][i]+big_list[-2][i]+big_list[-3][i]+big_list[-4][i]+big_list[-5][i] for i in range(len(little_list))]
            ten_year_average = [val/5 for val in average_total] #CHANGE THE NUMBER BACK TO 10 WHEN WORKING
            plt.plot(ten_year_average, temps)
        #print("year is ", year)
        heat[year] = temps
        
        fit = 302.3 - 45.3*(np.sin(lats)**2) #coakley model
        #plt.figure(1)
        plt.plot(lats, temps)
        plt.plot(lats, n*[273.15]) #zero degrees celcius
        plt.plot(lats, fit)
        plt.axis([-1.6, 1.6, 200, 350])
        plt.xlabel(f'latitude, $year = {day/365}$')
        plt.ylabel('Temps')
        plt.show()
        plt.pause(0.001) 
        little_list = np.empty(n) # to clear list at end of use

plt.figure(figsize=(10,10))
#y_axis_labels = np.linspace(-90, 90, 144)
heat_map = sns.heatmap(np.array(heat).T, linewidth = 0 , annot = False, cbar_kws={'label': 'Temperature (K)'}, yticklabels = 6)
plt.title(f"Temperatures per year as model settles over ${years}$ years, at each latitude")
plt.xlabel("Year")
plt.ylabel("Latitude")
plt.show()

end = time.time()
print("Time elapsed:", end - start)