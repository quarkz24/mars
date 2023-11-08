# -*- coding: utf-8 -*-
"""
Created on Tue Nov  7 11:52:08 2023

@author: josep
"""

import numpy as np
import matplotlib.pyplot as plt
from flux_function import flux

au = 149597870.7*1000
axial_tilt_earth = 0.409
q_earth = 1360
year_length_earth = 365

mars_a = 1.524*au
q_mars = 589
axial_tilt_mars = 0.438
year_length_mars = 687

days = 687
lats = np.linspace(0, 1.57, 10)

flux_earth = flux(lats, 365, q_earth, au, axial_tilt_earth, year_length_earth)
flux_mars = flux(lats, 687, q_mars, mars_a, axial_tilt_mars, year_length_mars)

for latitude in range(len(flux_earth)):
    plt.plot(np.arange(0, year_length_earth, 1), flux_earth[latitude], linestyle = '--', label = f'$latitude = {lats}$')
    plt.ylabel("Diurnally averaged incident solar flux", fontsize = 10)
    plt.xlabel("Days of the year")
    plt.title("Diurnally averaged incident solar flux per day of the year (Earth)")

plt.show()
for latitude in range(len(flux_mars)):
    plt.plot(np.arange(0, year_length_mars, 1), flux_mars[latitude], linestyle = '--', label = f'$latitude = {lats}$')
    plt.ylabel("Diurnally averaged incident solar flux", fontsize = 10)
    plt.xlabel("Days of the year")
    plt.title("Diurnally averaged incident solar flux per day of the year (Mars)")

#flux_mars = flux(lats, days, q_mars, mars_a, axial_tilt_mars, year_length_mars)
#plt.plot(np.arange(0, 687, 1), flux_mars, linestyle = '--', label = f'$latitude = {lats}$')
#plt.show()