# -*- coding: utf-8 -*-
"""
Created on Tue Oct 17 11:20:20 2023

@author: josep
"""

import numpy as np
import matplotlib.pyplot as plt
from flux_function import flux


L = np.pi  #Grosor de la pared en metros
n = 145   #Numero de nodos utilizados   
T0 = 300  #Temperatura inicial
dx = L/n
d = 0.0001 #Difusividad termica K/(Rho*C_p) 
t_final = 365  #Final temp in days - originally in seconds
dt = 1
c = (5.25e9)
sb = 5.670374419e-8


x = np.linspace(-np.pi + dx/2, np.pi-dx/2, n)
T = np.ones(n)*T0
dTdt = np.empty(n)

print("dx in radians: ", dx, ". Initial temp array of length", len(T), ": ", T)

t = np.arange(0, t_final, dt)

print("days to evaluate at: ", t)

solar = flux(x, t_final)


for j in range(1,len(t)):
    
    plt.clf()

    for i in range(1, n-1):

        #ir cooling function
        tir = 0.79*((T[i]/273)**3)
        ir = (sb*T[i])/(1+0.75*tir)
        
        #albedo
        a = 0.525 - 0.245*np.tanh((T[i]-268)/5)
        
        #dTdt[i] = alpha*(-(T[i]-T[i-1])/dx**2+(T[i+1]-T[i])/dx**2)
        dTdt[i] = (2*dt/c)*(solar[i][j]*(1-a) - d*np.tan(x[i])*((T[i+1]-T[i-1])/2*dx) + d*((T[i+1]-2*T[i] + T[i-1])/(dx**2)) - ir)
        print("temp change at latitude ", i, " on day ", j, " is ", dTdt[i])
    #dTdt[0] = alpha*(-(T[0]-T1s)/dx**2+(T[1]-T[0])/dx**2)
    #dTdt[n-1] = alpha*(-(T[n-1]-T[n-2])/dx**2+(T2s-T[n-1])/dx**2)
    
    dTdt[n-1] = (2*dt/c)*(solar[i][j]*(1-a) - d*(T[i] - T[i-1])/(dx**2) - ir) # north pole
    dTdt[0] = (2*dt/c)*(solar[i][j]*(1-a) - d*((T[i+1] - T[i])/dx**2) - ir) #south pole

    T = T + dTdt*dt
    plt.figure(1)
    plt.plot(x,T)
    plt.axis([-1.6, 1.6, 0, 400])
    plt.xlabel('Latitude')
    plt.ylabel('Temperature (K)')
    plt.show()
    plt.pause(0.01)