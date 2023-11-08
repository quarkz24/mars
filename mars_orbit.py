# -*- coding: utf-8 -*-
"""
Created on Tue Nov  7 14:36:49 2023

@author: josep
"""

import numpy as np
import matplotlib.pyplot as plt

##constants
semimaj_mars = 227.956*0.00668458712
semimin_mars = 226.957*0.00668458712
eccentricity_mars = 0.0935
q_0 = 1360

semimaj_earth = 149.5978707*0.00668458712
semimin_earth = 149.579110*0.00668458712
eccentricity_earth = 0.017

##lists
#x = np.linspace(-300, 300, 100)
t = np.linspace(0, 2*np.pi, 100) # for parametrisation of elliptical orbit
Ls = np.linspace(0, 2*np.pi, 100) # in rad
#Ls = np.linspace(0, 4*np.pi, 200) # in rad

## MARS ############################################################
y_vals = []
q_vals = []
y = tuple(zip(semimaj_mars*np.cos(t), semimin_mars*np.sin(t)))
focus_x = eccentricity_mars*semimaj_mars
xs = [val[0] for val in y]
ys = [val[1] for val in y]

##plotting
ax = plt.subplot(111)
ax.set_aspect('equal')
plt.plot(xs, ys, linewidth = 1, color = "red")
plt.xlabel("Distance from centre of Mars' orbit along major axis (AU)", fontsize = 8)
plt.ylabel("Distance from centre of Mars' orbit along minor axis (AU)", fontsize = 8)
#plt.title("Mars' orbit centred at (0, 0)")
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.plot(xs, [0]*len(xs), linewidth = 1, color = "black") # major axis
plt.plot([0]*len(ys), ys, linewidth = 1, color = "black") # minor axis
plt.scatter(focus_x, 0, marker = "x", color = "orange") # the sun
plt.show()
    
r = (semimaj_mars*(1 - eccentricity_mars**2))/(1 + eccentricity_mars*np.cos(Ls)) #needs to be in AU, currently in 10^6 k
q = q_0/r**2 
mean_q_mars = sum(q)/len(r)
print("mean q mars: ", sum(q)/len(r))

ax = plt.subplot(111)
plt.plot(Ls, q, linewidth = 1)
plt.plot(Ls, [mean_q_mars]*len(r), linewidth = 1, color = 'black')
plt.xlabel("Areocentric (Mars-centric) longitude $L_{s}$ in radians", fontsize = 10)
plt.ylabel("Solar irradiation $q (Wm^{-2})$", fontsize = 10)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.show()

ax = plt.subplot(111)
plt.plot(r, q, linewidth = 1)
#plt.plot()
plt.xlabel("Mars-Sun distance (AU)", fontsize = 10)
plt.ylabel("Solar irradiation $q (Wm^{-2})$", fontsize = 10)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.show()

########################################################
## EARTH ###############################################

y_vals = []
q_vals = []
y = tuple(zip(semimaj_earth*np.cos(t), semimin_earth*np.sin(t)))
focus_x = eccentricity_earth*semimaj_earth
xs = [val[0] for val in y]
ys = [val[1] for val in y]

##plotting
ax = plt.subplot(111)
ax.set_aspect('equal')
plt.plot(xs, ys, linewidth = 1, color = "blue")
plt.xlabel("Distance from centre of Earths orbit along major axis (AU)", fontsize = 8)
plt.ylabel("Distance from centre of Earths orbit along minor axis (AU)", fontsize = 8)
#plt.title("Mars' orbit centred at (0, 0)")
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.plot(xs, [0]*len(xs), linewidth = 1, color = "black") # major axis
plt.plot([0]*len(ys), ys, linewidth = 1, color = "black") # minor axis
plt.scatter(focus_x, 0, marker = "x", color = "orange") # the sun
plt.show()
    
r = (semimaj_earth*(1 - eccentricity_earth**2))/(1 + eccentricity_earth*np.cos(Ls)) #needs to be in AU, currently in 10^6 k
q = q_0/r**2 
mean_q_earth = sum(q)/len(r)
print("mean q earth: ", sum(q)/len(r))

ax = plt.subplot(111)
plt.plot(Ls, q, linewidth = 1)
plt.plot(Ls, [mean_q_earth]*len(r), linewidth = 1, color = 'black')
plt.xlabel("Solar longitude $L_{s}$ in radians", fontsize = 10)
plt.ylabel("Solar irradiation $q (Wm^{-2})$", fontsize = 10)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.show()

ax = plt.subplot(111)
plt.plot(r, q, linewidth = 1)
#plt.plot()
plt.xlabel("Earth-Sun distance (AU)", fontsize = 10)
plt.ylabel("Solar irradiation $q (Wm^{-2})$", fontsize = 10)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.show()

######################################
y_vals_mars = []
y_mars = tuple(zip(semimaj_mars*np.cos(t), semimin_mars*np.sin(t)))
focus_x_mars = eccentricity_mars*semimaj_mars
xs_mars = [val[0] for val in y_mars]
ys_mars = [val[1] for val in y_mars]

y_vals_earth = []
y_earth = tuple(zip(semimaj_earth*np.cos(t), semimin_earth*np.sin(t)))
focus_x_earth = eccentricity_earth*semimaj_earth
xs_earth = [val[0] for val in y_earth]
ys_earth = [val[1] for val in y_earth]

ax = plt.subplot(111)
ax.set_aspect('equal')
plt.plot(xs_mars, ys_mars, linewidth = 1, color = "red")
plt.plot(xs_earth, ys_earth, linewidth = 1, color = "blue")
plt.xlabel("Distance from centre of orbit along major axis (AU)", fontsize = 8)
plt.ylabel("Distance from centre of orbit along minor axis (AU)", fontsize = 8)
#plt.title("Mars' orbit centred at (0, 0)")
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.plot(xs, [0]*len(xs), linewidth = 1, color = "black") # major axis
plt.plot([0]*len(ys), ys, linewidth = 1, color = "black") # minor axis
plt.scatter(focus_x_mars, 0, marker = "x", color = "orange") # the sun
plt.scatter(focus_x_earth, 0, marker = "x", color = "black")
plt.show()



