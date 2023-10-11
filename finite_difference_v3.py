# -*- coding: utf-8 -*-
"""
Created on Wed Oct 11 12:48:07 2023

@author: josep
"""

import numpy as np
import matplotlib.pyplot as plt
from flux_function import flux

L = 0.1  
n = 142   
T0 = 0  
T1s = 40 
T2s = 20 
#dx = L/n

q = 1360
#d = 5.394*100 #diffusivity, from 'habitable climates'
d = 0.005394
a = 149597870.7*1000 
c = 40 * 5.25e9 #heat cap, initially constant, of water planet 
dt = 1
dx = 0.0218 #in degrees->radians
initial_t = 300
sb = 5.670374419e-8
albedo = 0.5 #actually a function

t_final = 365  #Tiempo final en segundos
dt = 1

x = np.linspace(-1.57+dx/2, 1.57-dx/2, n)
T = np.ones(n)*T0
dTdt = np.empty(n)
t = np.arange(0, t_final, dt)

solar = flux(x)

for j in range(1,len(t)):
    plt.clf()
    for i in range(1, n-1): # remember, i took out outgoing flux
    
        #dTdt[i] = (2*dt/c)*(d*((T[i+1] - 2*T[i] + T[i-1])/dx) - d*np.tan(x[i])*((T[i+1] - T[i-1])/2*dx)) + T[i-1] 
        #dTdt[i] = (2*dt/c)*(solar[i][j]*(1-albedo) + d*((T[i+1] - 2*T[i] + T[i-1])/dx) - d*np.tan(x[i])*((T[i+1] - T[i-1])/2*dx)) + T[i-1]
        
        dTdt[i] = d*(-(T[i]-T[i-1])/dx**2+(T[i+1]-T[i])/dx**2)
        
    dTdt[n-1] = solar[i][j]*(1-albedo) - ((T[i] - T[i-1])/dx**2) #dT/dlamb is zero, second deriv is back diff
    dTdt[0] = solar[i][j]*(1-albedo) + ((T[i+1] - T[i])/dx**2) #dT/dlamb is zero, second deriv is forward diff
    
    #dTdt[0] = d*(-(T[0]-T1s)/dx**2+(T[1]-T[0])/dx**2)
    #dTdt[n-1] = d*(-(T[n-1]-T[n-2])/dx**2+(T2s-T[n-1])/dx**2)

    T = T + dTdt*dt

    plt.figure(1)
    plt.plot(x,T)
    plt.axis([-1.57, 1.57, 0, 700])
    plt.xlabel('Distance (m), day ', )
    plt.ylabel('Temperature (C)')
    plt.show()
    plt.pause(0.01) 