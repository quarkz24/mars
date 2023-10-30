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

start = time.time()

##constants
def constants():
    second_to_day = 24*60*60
    q = 1362
    d = (5.394e-1)*second_to_day #diffusivity
    #cl = (5.25e6) #heat cap, initially constant, of dry planet
    cl = 5.25e7
    sb = 5.670374419e-8 #stefan-boltzmann constant
    n = 144
    del_t = 1
    del_lamb = 3.1415/n #in degrees->radians
    initial_temp = 300
    years = 100
    partial_years = 5*years #to plot changing seasons
    final_time = 365*years #years * days
    return second_to_day, q, d, cl, sb, n, del_t, del_lamb, initial_temp, years, partial_years, final_time

second_to_day, q, d, cl, sb, n, del_t, del_lamb, initial_temp, years, partial_years, final_time = constants()

def list_of_lists():
    temps = [initial_temp]*n #144 length list of initial temps
    prev_temps = [initial_temp]*n
    lats = np.linspace(-1.57, 1.57, n) #both start and end inclusive
    solar = flux(lats, final_time) #returns a list of lists, indexable, returns list of flux over a year per lat
    temp_diff = np.empty(n) # empty list
    heat = [[0] for i in range(partial_years)]
    cumu = np.zeros(n)
    true_diff = np.zeros(years)
    percent_diff = np.zeros(years)
    self_diff = np.zeros(years)
    av_list = []
    coakley_fit = 302.3 - 45.3*(np.sin(lats)**2)
    return temps, prev_temps, lats, solar, temp_diff, heat, cumu, true_diff, percent_diff, self_diff, av_list, coakley_fit

temps, prev_temps, lats, solar, temp_diff, heat, cumu, true_diff, percent_diff, self_diff, av_list, coakley_fit = list_of_lists()

##functions

def heatcap(capacity):
    if temps[y] >= 273:
        c = cl*40*0.7 + cl*0.3
    elif 263 < temps[y] < 273:
        c = 9.2*cl
    else:
        c = 2*cl
    return c

def heatmap(array):
    plt.figure(figsize=(10,10))
    plt.imshow(np.array(array).T, extent=[0, years, 90, -90], aspect = 'auto', cmap = 'plasma')
    plt.colorbar(label = "Temperature, K", orientation = "vertical")
    plt.title(f"Temperatures per year as model settles over ${years}$ years, at each latitude")
    plt.xlabel("Year")
    plt.ylabel("Latitude")
    plt.show()

def outgoing_1(index):
    infrared = ((sb*(temps[y])**4)/(1.75))*second_to_day
    albedo = 0.5 - 0.2*np.tanh((temps[index]-268)/5)
    return infrared, albedo

def outgoing_2(index):
    tir = 0.79*((temps[y]/273)**3)
    infrared = ((sb*(temps[y])**4)/(1+0.75*tir))*second_to_day
    albedo = 0.525 - 0.245*np.tanh((temps[y]-268)/5)
    return infrared, albedo

def outgoing_3(index):
    const_a = 2.033e4
    const_b = 2.094
    infrared = const_a + const_b*temps[index]
    albedo = 0.475 - 0.225*np.tanh((temps[y]-268)/5)
    return infrared, albedo

def convergence_test_1a(diff_list):
    plt.figure(figsize=(10,10))
    plt.plot(np.arange(0, years), diff_list)
    plt.title("Convergence test 1a: difference to N&C fit over years")
    plt.xlabel("Year")
    plt.ylabel("Mean difference, K")
    plt.show()
    
def convergence_test_1b(diff_list):
    plt.figure(figsize=(10,10))
    plt.plot(np.arange(0, years), diff_list)
    plt.title("Convergence test 1b: percentage difference to N&C fit over years")
    plt.xlabel("Year")
    plt.ylabel("Mean difference, %")
    plt.show()

def convergence_test_2(diff_list):
    plt.figure(figsize=(10,10))
    plt.plot(np.arange(0, years), diff_list)
    plt.title("Convergence test 2: $\delta$ = $|T_{year} - T_{year-1}|$")
    plt.xlabel("Year")
    plt.ylabel("$\delta$, K")
    plt.show()
    
def temp_difference(index): 
    if index != 0 and index != n-1: #most latitudes
        ir, a = outgoing_2(y)
        temp_diff[index] = (del_t/heatcap(cl))*(second_to_day*solar[y][day]*(1-a) - d*np.tan(lats[index])*((temps[index+1]-temps[index-1])/2*del_lamb) + d*((temps[index+1]-2*temps[index] + temps[index-1])/(del_lamb**2)) - ir)
        
    elif index == 0: #south pole
        ir, a = outgoing_2(0)
        temp_diff[0] = (del_t/heatcap(cl))*(second_to_day*solar[0][day]*(1-a) + d*((temps[0+1] - temps[0])/del_lamb**2) - ir)
        
    elif index == n-1: #north pole
        ir, a = outgoing_2(n-1)
        temp_diff[-1] = (del_t/heatcap(cl))*(second_to_day*solar[-1][day]*(1-a) - d*(temps[-1] - temps[-1-1])/(del_lamb**2) - ir)
        
    return temp_diff
        
##looping

for day in range(0, final_time, del_t): #gives list 0->364, 365 entries
    for y in range(len(lats)): #0->143, 144 entries
        
        change = temp_difference(y)
    
    temps = change + temps
    cumu = cumu + temps
    
    if day % 73 == 0:
        x = int(day/73)
        heat[x] = temps
    
    if day % 365 == 0:
        year = int(day/365)
        
        true_diff[year] = np.average(abs(coakley_fit - temps))
        percent_diff[year] = np.average(abs((coakley_fit[y]-temps[y])/coakley_fit[y])*100)
        self_diff[year] = np.average(abs(temps - prev_temps))
        prev_temps = temps
        
        av_list.append([val/365 for val in cumu])
        cumu = np.zeros(n) #clear cumulative list for next year
        
        plt.figure(1)
        if year > 9:
            add = [sum(i) for i in zip(av_list[-1], av_list[-2], av_list[-3], av_list[-4], av_list[-5], av_list[-6], av_list[-7], av_list[-8], av_list[-9], av_list[-10])]
            div = [val/10 for val in add]
            plt.plot(lats, div, label = "10 year average", color = 'blue', linestyle = '--') #ten yr average

        plt.plot(lats, temps, label = "finite diff method", color = 'blue') #temps per year
        plt.plot(lats, n*[273.15], color = 'black') #zero degrees celcius
        plt.plot(lats, coakley_fit, label = "N&C fit", color = 'red') #coakley model
        plt.axis([-1.57, 1.57, 200, 350])
        plt.xlabel(f'latitude, $year = {day/365}$')
        plt.ylabel('Temps')
        plt.legend()
        plt.show()
        #plt.pause(0.01) 

heatmap(heat)
convergence_test_1a(true_diff)
convergence_test_1b(percent_diff)
convergence_test_2(self_diff)

end = time.time()
print("Time elapsed:", end - start)