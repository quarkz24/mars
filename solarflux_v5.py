
# coding: utf-8

# In[3]:


##V5 adapted from V4, reattempting with more Taylor terms to see if that fixes the strange shape they originally made...
##...while also being more efficient (less trig evaluations)
##this doesn't work lol

##NORTHERN HEMISPHERE
##PACKAGES AND CONSTANTS AND FUNCTIONS

##import packages

import numpy as np
import matplotlib.pyplot as plt
import time

start = time.time()

##define constants (or things that are constant for now)
##my trig functions work in radians, as does np.sin etc, so latitude in rad (10deg = 0.17rad)

a = 149597870.7*1000
q = 1360

##taylor functions for efficiency - use radians

def tan(x):
	return x + (x**3)/6 + (2*x**5)/15
	
def sin(x):
	return x - (x**3)/6 + (x**5)/120 - x**7/5040
	
def cos(x):
	return 1 - (x**2)/2 + (x**4)/24 - x**6/720

latlist = [0.017, 0.17, 0.35, 0.52, 0.70, 0.87, 1.05, 1.22, 1.40, 1.56]
for lat in latlist: #in radians

    sintheta = np.sin(lat)
    costheta = np.cos(lat)
    tantheta = np.tan(lat)
    
    flux = []

    for day in range(0, 365, 1):
    
        #2*3.1415/365 to convert day to position around sun
        #-0.409rad is peak delta, remains in rad to feed into my trig functions
    
        delta = -0.409*cos((2*3.1415/365)*(day + 10)) # why is delta returning as anything larger than 0.409?
        big_h = np.arccos(-tantheta * tan(delta))
        solar = (q/3.1415) * (((149597870.7*1000)/a)**2) * (big_h*sintheta*sin(delta) + costheta*cos(delta)*sin(big_h))
        flux.append(solar)
        
    plt.plot(np.arange(0, 365, 1), flux, linestyle = '--')

#print("delta list is:", d)
#print("big H list is:", h)
#print("flux:", flux)

plt.ylabel("Diurnally averaged incident solar flux", fontsize = 10)
plt.xlabel("Days of the year")
plt.title("Diurnally averaged incident solar flux per day of the year, beginning Jan 1 (Northern hemisphere)")

end = time.time()
print("Time elapsed:", end - start)


# In[ ]:


##V4 adapted from V3, want to add in loop to iterate over 10 deg between 0 -> 90 and 0 -> -90

##SOUTHERN HEMISPHERE
##PACKAGES AND CONSTANTS AND FUNCTIONS

##import packages

import numpy as np
import matplotlib.pyplot as plt
import time

start = time.time()

##define constants (or things that are constant for now)
##my trig functions work in radians, as does np.sin etc, so latitude in rad (10deg = 0.17rad)

a = 149597870.7*1000
q = 1360

latlist2 = [-0.017, -0.17, -0.35, -0.52, -0.70, -0.87, -1.05, -1.22, -1.40, -1.56]

for lat in latlist2: #in radians

    sintheta = np.sin(lat)
    costheta = np.cos(lat)
    tantheta = np.tan(lat)
    
    flux = []

    for day in range(0, 365, 1):
    
        #2*3.1415/365 to convert day to position around sun
        #-0.409rad is peak delta, remains in rad to feed into my trig functions
    
        delta = -0.409*np.cos((2*3.1415/365)*(day + 10)) # why is delta returning as anything larger than 0.409?
        big_h = np.arccos(-tantheta * np.tan(delta))
        solar = (q/3.1415) * (((149597870.7*1000)/a)**2) * (big_h*sintheta*np.sin(delta) + costheta*np.cos(delta)*np.sin(big_h))
        flux.append(solar)
        
    plt.plot(np.arange(0, 365, 1), flux, linestyle = '--')

#print("delta list is:", d)
#print("big H list is:", h)
#print("flux:", flux)

plt.ylabel("Diurnally averaged incident solar flux", fontsize = 10)
plt.xlabel("Days of the year")
plt.title("Diurnally averaged incident solar flux per day of the year, beginning Jan 1 (Southern hemisphere)")

end = time.time()
print("Time elapsed:", end - start)

