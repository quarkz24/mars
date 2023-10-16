# -*- coding: utf-8 -*-
"""
Created on Mon Oct 16 09:30:59 2023

@author: josep
"""

##trying to do temps[days-1] outside of lat loop
##this just breaks it and the only temps evolving are the ones at the poles - whichare the ones i didnt change
##import packages

import numpy as np
import matplotlib.pyplot as plt
import time
from flux_function import flux

start = time.time()

##constants

q = 1360*24*60*60
d = 0.1 #diffusivity, placeholder, inaccurate
a = 0.1 #also false
c = (5.25e9)*40 #heat cap, initially constant, of water planet 
sb = 5.670374419e-8

n = 144
del_t = 1
del_lamb = 1.25 * (3.1415/180) #in degrees->radians
initial_temp = 300
final_time = 365*40

##arrays+lists

temps = [initial_temp]*n #144 length list of initial temps
lats = np.linspace(-1.57, 1.57, n) #both start and end inclusive
solar = flux(lats) #returns a list of lists, indexable, returns list of flux over a year per lat
temp_diff = np.empty(n)
##do i need a new list? a different method - feed in 300s, populate new list
print("lats: ",len(lats))
print("temps: ", len(temps))
print("incoming flux: ", len(solar)) #144 lists of flux, for 144 latitudes
print("incoming flux in south pole winter: ", solar[0][180]) #just verifying im getting an answer i expect

##looping

for day in range(0, final_time, del_t): #gives list 0->364, 365 entries
    for y in range(len(lats)-1): #0->143, 144 entries
        #ir cooling function
        tir = 0.79*((temps[y]/273)**3)
        ir = (sb*temps[y])/(1+0.75*tir)
        
        #x = temps[day-1]
        #x = solar[y][day]
        temp_diff[y] = (2*del_t/c)*(solar[y][day]*(1-a) - d*np.tan(lats[y])*((temps[y+1]-temps[y-1])/2*del_lamb) + d*((temps[y+1]-2*temps[y] + temps[y-1])/(del_lamb**2)) - ir)
    
    temps[len(lats)-1] = (2*del_t/c)*(solar[y][day]*(1-a) - d*(temps[y] - temps[y-1])/(del_lamb**2)) # north pole
    temps[0] = (2*del_t/c)*(solar[y][day]*(1-a) - d*((temps[y+1] - temps[y])/del_lamb**2)) #south pole - what is day-1 for first day?
    
    temps = temp_diff + temps
    
    print("temps for day", day, ":", temps)
    
    #plt.figure(1)
    #plt.plot(lats, temps)
    #plt.axis([-1.57, 1.57, 0, 700])
    #plt.xlabel(f'latitude, $day = {day}$') #label = f'$latitude = {lat}$'
    #plt.ylabel('Temps')
    #plt.show()
    #plt.pause(0.01) 
    
end = time.time()
print("Time elapsed:", end - start)