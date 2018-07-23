import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon

# read in the data, which has header Dist Mv CL Typ LTef logg Age Mass B-I U-B V-I V-K V  mux  muy  Vr UU  VV  WW [Fe/H] l   b   Av Mbol Mbol
color, vimag, vmag, v = np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/MWFG_stars.txt', usecols=(8, 10, 12, 15), unpack=True)

# making a CMD
def magnitudes(filter1, filter2):  # I band
	return filter2 - filter1

mag = magnitudes(vimag, vmag)

#sort the points into 4 age groups
MWFG_MS_v=[]
MWFG_AGy_v=[]
MWFG_AGo_v=[]
MWFG_RG_v=[]

# MS
for i in range(len(color)):
	if color[i]<=1:
		MWFG_MS_v.append(v[i])

#polygon creates a shape that corresponds to the area of the CMD that contains young AGB stars
polygon= Polygon([(3.5,18), (8,18), (2.7, 20.4), (3.5,20.5), (8,20.5)])
for i in range(len(color)):
	point= Point(color[i], mag[i])
	if polygon.contains(point)==True:
		MWFG_AGy_v.append(v[i])

#polygon creates a shape that corresponds to the area of the CMD that contains old AGB stars
polygon= Polygon([(2.7,20.4), (7.5, 23), (8, 20.5)])
for i in range(len(color)):
	point= Point(color[i], mag[i])
	if polygon.contains(point)==True:
		MWFG_AGo_v.append(v[i])

#polygon creates a shape that corresponds to the area of the CMD that contains RGB stars
polygon= Polygon([(2,23), (2.7, 20.4), (7.5,23)])
for i in range(len(color)):
	point= Point(color[i], mag[i])
	if polygon.contains(point)==True:
		MWFG_RG_v.append(v[i])

# plotting
plt.hist(MWFG_MS_v, color='b')
plt.xlabel('Radial Velocity (km/s)')
plt.savefig('/Users/amandaquirk/Desktop/MWFG_MS_velocities.png')
plt.close()

plt.hist(MWFG_AGy_v, color='m')
plt.xlabel('Radial Velocity (km/s)')
plt.savefig('/Users/amandaquirk/Desktop/MWFG_AGy_velocities.png')
plt.close()

plt.hist(MWFG_AGo_v, color='m')
plt.xlabel('Radial Velocity (km/s)')
plt.savefig('/Users/amandaquirk/Desktop/MWFG_AGo_velocities.png')
plt.close()

plt.hist(MWFG_RG_v, color='r')
plt.xlabel('Radial Velocity (km/s)')
plt.savefig('/Users/amandaquirk/Desktop/MWFG_RG_velocities.png')
plt.close()



