# -*- coding: utf-8 -*-
"""
Created on Fri Nov 17 12:12:26 2023

@author: josep
"""

##import packages

import numpy as np
import matplotlib.pyplot as plt
import time
import math
from flux_function_2 import flux

start = time.time()

##constants
def constants():
    #constants
    q = 1362 #W m^-2 = J s^-1 m^-2
    d = (5.394e-1) #diffusivity, J m^-2 s^-1 K^-1
    #d = 0.01
    cl = (5.25e6) #land heat cap, J m^-2 K^-1
    sb = 5.670374419e-8 #stefan-boltzmann, J s^-1 m^-2 K^-4
    
    #earth constants
    semimaj_earth = 149597870.7*1000
    eccentricity_earth = 0.017
    axial_tilt_earth = 0.409 #in rad
    perihelion_earth = 282.02 #in deg, converted to rad in flux_function
    perihelion_earth = 102.04
    year_length_earth = 365
    
    #mars constants
    semimaj_mars = 227956000*1000*0.7
    #semimaj_mars = 149597870.7*1000
    eccentricity_mars = 0.0935
    axial_tilt_mars = 0.438
    #axial_tilt_mars = np.pi/2
    perihelion_mars = 251
    perihelion_mars = 71
    year_length_mars = 687
    
    #conversions - are either of these global anyway?
    deg_to_rad = np.pi/180
    au = 149597870.7*1000
    
    #variables
    n = 32
    n = 31
    del_t = int((24*60*60))
    del_lamb = np.pi/(n-1) #in degrees->radians
    initial_temp = 300
    years = 20
    
    del_t_mars = int(del_t + (39*60) + 35)
    final_time_mars = year_length_mars*years
    seconds_mars = final_time_mars*del_t_mars
    partial_years_mars = 3*years
    partial_years_mars = 229*years
    
    #functionals
    partial_years = 5*years #to plot changing seasons
    final_time = year_length_earth*years #years * days
    seconds = final_time*(24*60*60)

    return q, d, cl, sb, semimaj_earth, eccentricity_earth, axial_tilt_earth, perihelion_earth, year_length_earth, semimaj_mars, eccentricity_mars, axial_tilt_mars, perihelion_mars, year_length_mars, deg_to_rad, au, n, del_t, del_lamb, initial_temp, years, del_t_mars, final_time_mars, seconds_mars, partial_years_mars, partial_years, final_time, seconds 

q, d, cl, sb, semimaj_earth, eccentricity_earth, axial_tilt_earth, perihelion_earth, year_length_earth, semimaj_mars, eccentricity_mars, axial_tilt_mars, perihelion_mars, year_length_mars, deg_to_rad, au, n, del_t, del_lamb, initial_temp, years, del_t_mars, final_time_mars, seconds_mars, partial_years_mars, partial_years, final_time, seconds = constants()

def list_of_lists():
    #constant lists
    temps = [initial_temp]*n 
    prev_temps = [initial_temp]*n
    lats = np.linspace(-1.57, 1.57, n) #both start and end inclusive
    temp_diff = np.zeros(n) # empty list
    heat = [[0] for i in range(partial_years_mars)]
    cumu = np.zeros(n)
    true_diff = np.zeros(years)
    percent_diff = np.zeros(years)
    self_diff = np.zeros(years)
    av_list = []
    albedo_list = []
    
    #planet dependent
    solar = flux(lats, years, year_length_mars, semimaj_mars, eccentricity_mars, perihelion_mars, axial_tilt_mars)
    #solar = flux(lats, years, year_length_mars, semimaj_earth, eccentricity_earth, perihelion_earth, axial_tilt_earth)
        
    Ls = np.linspace(0, 2*np.pi, year_length_earth)
    coakley_fit = 302.3 - 45.3*(np.sin(lats)**2)
    return temps, prev_temps, lats, temp_diff, heat, cumu, true_diff, percent_diff, self_diff, av_list, albedo_list, solar, Ls, coakley_fit

temps, prev_temps, lats, temp_diff, heat, cumu, true_diff, percent_diff, self_diff, av_list, albedo_list, solar, Ls, coakley_fit = list_of_lists()
#print(len(solar[0])) #6870
#print(solar[0][years*year_length_mars])

##functions
def diffusion():
    p0 = 1013 #millibar
    d0 = 0.58
    cp0 = 1000
    m0 = 28 #g per mol
    omega0 = 7.27e-5 #rad per second
    
    p = 6.518 #millibar
    cp = 1016
    m = 43.3 #dont know this yet
    omega = 7.09e-5
    
    d = d0*(p/p0)*(cp/cp0)*((m0/m)**2)*((omega0/omega)**2)
    #d=0.02
    #d=0.05
    #d=0.1
    #d=0.2
    return d

d = diffusion()

def heatcap(capacity): #applying proportional lands to other temps is bad
    fraction = 1 - np.exp((temps[y] - 273)/10) #fraction of ocean which is ice
    if temps[y] >= 273:
        c = cl*40*0.7 + cl*0.3
        c = cl
    elif 263 < temps[y] < 273:
        c = 0.7*(fraction*9.2*cl + (1-fraction)*40*cl) + 0.3*cl
        c = cl
    else:
        c = 0.7*2*cl + 0.3*cl #used to be 2
        c = 2*cl
    return c

def heatcap(capacity): 

    if temps[y] >= 273:
        c = cl
    elif temps[y] < 273:
        #c = 0.7*2*cl + 0.3*cl
        c=2*cl
    return c

#def heatcap(capacity): 
#
 #   if y == 0 or y == 1 or y == n-1 or y ==n:
  #      c = 2*cl
   # else:
    #    c = cl
    #return c 

#def heatcap(capacity):
    #c = cl*40
    #return c

def heatmap(array):
    plt.figure(figsize=(10,10))
    plt.imshow(np.array(array).T, extent=[0, years, 90, -90], aspect = 'auto', cmap = 'plasma')
    plt.colorbar(label = "Temperature, K", orientation = "vertical")
    plt.title(f"Temperatures per year as model settles over ${years}$ years, at each latitude")
    plt.xlabel("Year")
    plt.ylabel("Latitude")
    plt.show()
    
    plt.figure(figsize=(10,10))
    plt.imshow(np.array(array[-25:]).T, extent=[years - 5, years, 90, -90], aspect = 'auto', cmap = 'plasma')
    plt.imshow(np.array(array[-1145:]).T, extent=[years - 5, years, 90, -90], aspect = 'auto', cmap = 'plasma')
    plt.colorbar(label = "Temperature, K", orientation = "vertical")
    plt.title(f"Latest 5 years of temperature in climate model, evolved over ${years}$ years")
    plt.xlabel("Year")
    plt.ylabel("Latitude")
    plt.show()

def outgoing_1(index):
    infrared = ((sb*(temps[index])**4)/(1.75))
    albedo = 0.5 - 0.2*np.tanh((temps[index]-268)/5)
    return infrared, albedo

def outgoing_2(index):
    tir = 0.79*((temps[index]/273)**3)
    infrared = ((sb*(temps[index])**4)/(1+0.75*tir))
    albedo = 0.525 - 0.245*np.tanh((temps[index]-268)/5)
    return infrared, albedo

def outgoing_3(index):
    const_a = 2.033e4
    const_b = 2.094
    infrared = const_a + const_b*temps[index]
    albedo = 0.475 - 0.225*np.tanh((temps[index]-268)/5)
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
        ir, a = outgoing_2(index)
        second_order_term = d*(temps[index+1] - 2*temps[index] + temps[index-1])/(del_lamb**2)
        first_order_term = -d*np.tan(lats[index])*((temps[index+1] - temps[index-1])/(2*del_lamb))
        rad_term = solar[index][day]*(1-a) - ir
        
        #if day % 687 == 0:
            #print("year: ", day/687)
            #print("latitude: ", lats[index])
            #print("rad terms: ", rad_term)
            #print("albedo: ", a)
            #print("ir outgoing: ", ir)
            #print("ir incoming: ", solar[index][day])
            #print(" ")
            
        temp_diff[index] = (second_order_term + first_order_term + rad_term)*(del_t/heatcap(cl))
        
    elif index == 0: #south pole
        ir, a = outgoing_2(0)
        second_order_term = d*((temps[0+1] - temps[0])/del_lamb**2)
        rad_term = solar[0][day]*(1-a) - ir

        temp_diff[0] = (second_order_term + rad_term)*(del_t/heatcap(cl))
        #print("success for ", day)

    elif index == n-1: #north pole
        ir, a = outgoing_2(n-1)
        second_order_term = -d*((temps[-1] - temps[-1-1])/del_lamb**2)
        rad_term = solar[n-1][day]*(1-a) - ir

        temp_diff[-1] = (second_order_term + rad_term)*(del_t/heatcap(cl))

    return temp_diff
        

##looping

for s in range(0, seconds_mars, del_t_mars): #list 0->364, 365 entries
    day = math.floor(s/(24*60*60))
    day = math.floor(s/88775)
    
    for y in range(len(lats)): #0->143, 144 entries
            
        change = temp_difference(y)
    
    temps = change + temps
    cumu = cumu + temps
    
    #if np.isnan(temps[0]) == False:
     #   print("day is: ", day, "temps are: ", temps)

    if day % 3 == 0:
        x = int(day/3)
        heat[x] = temps
    
    if day % 687 == 0:
        year = int(day/687)
        #true_diff[year] = np.average(abs(coakley_fit - temps))
        #percent_diff[year] = np.average(abs((coakley_fit[y]-temps[y])/coakley_fit[y])*100)
        self_diff[year] = np.average(abs(temps - prev_temps))
        prev_temps = temps
        
        av_list.append([val/687 for val in cumu])
        cumu = np.zeros(n) #clear cumulative list for next year
        
        plt.figure(1)
        if year > 9:
            add = [sum(i) for i in zip(av_list[-1], av_list[-2], av_list[-3], av_list[-4], av_list[-5], av_list[-6], av_list[-7], av_list[-8], av_list[-9], av_list[-10])]
            div = [val/10 for val in add]
            plt.plot(lats, div, label = "10 year average", color = 'blue', linestyle = '--') #ten yr average

        plt.plot(lats, temps, label = "finite diff method", color = 'blue') #temps per year
        plt.plot(lats, n*[273.15], color = 'black') #zero degrees celcius
        plt.plot(lats, coakley_fit, label = "N&C fit", color = 'red') #coakley model
        plt.axis([-1.57, 1.57, 125, 225])
        plt.xlabel(f'latitude, $year = {day/687}$')
        plt.ylabel('Temps / K')
        plt.legend()
        plt.show()
        plt.pause(0.01) 
   
#print(heat)
#print(len(heat))
#print(len(heat[0]))
#print(len(heat[148]))
#print(len(heat[149]))
heatmap(heat)
#convergence_test_1a(true_diff)
#convergence_test_1b(percent_diff)
convergence_test_2(self_diff)

end = time.time()
print("Time elapsed:", end - start)