from astropy import units as u
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt 
import numpy as np 
from shapely.geometry import Point 
from shapely.geometry.polygon import Polygon
from matplotlib.ticker import MaxNLocator

#read in the data, which has header Dist Mv CL Typ LTef logg Age Mass B-I U-B V-I V-K V  mux  muy  Vr UU  VV  WW [Fe/H] l   b   Av Mbol Mbol
color, vimag, vmag, v=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/MWFG_stars.txt', usecols=(8,10,12,15,), unpack=True)

#making a CMD
def magnitudes(filter1, filter2): #I band
	return filter2 - filter1

mag=magnitudes(vimag, vmag)

#brings in M31 data to assign coordinates to MWFG
MS_xi, MS_eta=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/MS_smoothed_chemin.txt', usecols=(0,1), unpack=True)
AGy_xi, AGy_eta=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/AGy_smoothed_chemin.txt', usecols=(0,1), unpack=True)
AGo_xi, AGo_eta=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/AGo_smoothed_chemin.txt', usecols=(0,1), unpack=True)
RG_xi, RG_eta=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/RG_smoothed_chemin.txt', usecols=(0,1), unpack=True)

#below creates a CMD- this is just a visual check
# plt.scatter(color, mag, color='b', alpha=0.4)
# plt.xlabel('F475W - F814W (mag)')
# plt.xlim(-1,8)
# plt.ylim(18,23)
# plt.gca().invert_yaxis()
# plt.ylabel('F814W (mag)')
# plt.savefig('/Users/amandaquirk/Desktop/CDM_MWFG.png')
# plt.close()

#sort the points into 4 age groups
MWFG_MS_color=[]
MWFG_MS_mag=[]
MWFG_MS_v=[]
MWFG_AGy_color=[]
MWFG_AGy_mag=[]
MWFG_AGy_v=[]
MWFG_AGo_color=[]
MWFG_AGo_mag=[]
MWFG_AGo_v=[]
MWFG_RG_color=[]
MWFG_RG_mag=[]
MWFG_RG_v=[]

#MS
for i in range(len(color)):
	if color[i]<=1:
		MWFG_MS_color.append(color[i])
		MWFG_MS_mag.append(mag[i])
		MWFG_MS_v.append(v[i])

#polygon creates a shape that corresponds to the area of the CMD that contains young AGB stars
polygon= Polygon([(3.5,18), (8,18), (2.7, 20.4), (3.5,20.5), (8,20.5)])
for i in range(len(color)):
	point= Point(color[i], mag[i])
	if polygon.contains(point)==True:
		MWFG_AGy_color.append(color[i])
		MWFG_AGy_mag.append(mag[i])
		MWFG_AGy_v.append(v[i])

#polygon creates a shape that corresponds to the area of the CMD that contains old AGB stars
polygon= Polygon([(2.7,20.4), (7.5, 23), (8, 20.5)])
for i in range(len(color)):
	point= Point(color[i], mag[i])
	if polygon.contains(point)==True:
		MWFG_AGo_color.append(color[i])
		MWFG_AGo_mag.append(mag[i])
		MWFG_AGo_v.append(v[i])

#polygon creates a shape that corresponds to the area of the CMD that contains RGB stars
polygon= Polygon([(2,23), (2.7, 20.4), (7.5,23)])
for i in range(len(color)):
	point= Point(color[i], mag[i])
	if polygon.contains(point)==True:
		MWFG_RG_color.append(color[i])
		MWFG_RG_mag.append(mag[i])
		MWFG_RG_v.append(v[i])


# print('Length of MS={}, AGy={}, AGo={}, RG={}. Initial number of stars={}'.format(len(MWFG_MS_v), len(MWFG_AGy_v), len(MWFG_AGo_v), len(MWFG_RG_v), len(v)))

# #below does a visual check that the age group appears on the correct CMD
# plt.scatter(MWFG_MS_color, MWFG_MS_mag, color='b', alpha=0.4)
# plt.xlabel('F475W - F814W (mag)')
# plt.xlim(-1,8)
# plt.ylim(18,23)
# plt.gca().invert_yaxis()
# plt.ylabel('F814W (mag)')
# plt.savefig('/Users/amandaquirk/Desktop/CDM_MWFG_MS.png')
# plt.close()

# plt.scatter(MWFG_AGy_color, MWFG_AGy_mag, color='b', alpha=0.4)
# plt.xlabel('F475W - F814W (mag)')
# plt.xlim(-1,8)
# plt.ylim(18,23)
# plt.gca().invert_yaxis()
# plt.ylabel('F814W (mag)')
# plt.savefig('/Users/amandaquirk/Desktop/CDM_MWFG_AGy.png')
# plt.close()

# plt.scatter(MWFG_AGo_color, MWFG_AGo_mag, color='b', alpha=0.4)
# plt.xlabel('F475W - F814W (mag)')
# plt.xlim(-1,8)
# plt.ylim(18,23)
# plt.gca().invert_yaxis()
# plt.ylabel('F814W (mag)')
# plt.savefig('/Users/amandaquirk/Desktop/CDM_MWFG_AGo.png')
# plt.close()

# plt.scatter(MWFG_RG_color, MWFG_RG_mag, color='b', alpha=0.4)
# plt.xlabel('F475W - F814W (mag)')
# plt.xlim(-1,8)
# plt.ylim(18,23)
# plt.gca().invert_yaxis()
# plt.ylabel('F814W (mag)')
# plt.savefig('/Users/amandaquirk/Desktop/CDM_MWFG_RG.png')
# plt.close()

#assigning a random coordinate to the MWFG stars
MWFG_MS_xi=[]
MWFG_MS_eta=[]
MWFG_AGy_xi=[]
MWFG_AGy_eta=[]
MWFG_AGo_xi=[]
MWFG_AGo_eta=[]
MWFG_RG_xi=[]
MWFG_RG_eta=[]

for i in range(len(MWFG_MS_mag)):
	N=np.random.random_integers(0, len(MS_xi)-1)
	MWFG_MS_xi.append(MS_xi[N])
	MWFG_MS_eta.append(MS_eta[N])

for i in range(len(MWFG_AGy_mag)):
	N=np.random.random_integers(0, len(AGy_xi)-1)
	MWFG_AGy_xi.append(AGy_xi[N])
	MWFG_AGy_eta.append(AGy_eta[N])

for i in range(len(MWFG_AGo_mag)):
	N=np.random.random_integers(0, len(AGo_xi)-1)
	MWFG_AGo_xi.append(AGo_xi[N])
	MWFG_AGo_eta.append(AGo_eta[N])

for i in range(len(MWFG_RG_mag)):
	N=np.random.random_integers(0, len(RG_xi)-1)
	MWFG_RG_xi.append(RG_xi[N])
	MWFG_RG_eta.append(RG_eta[N])

#plotting
# plt.figure(figsize=(4, 6))
# plt.scatter(MWFG_MS_xi, MWFG_MS_eta, s=1, c=MWFG_MS_v, cmap='rainbow',vmin=-300,vmax=0)
# plt.gca().invert_xaxis()
# plt.xlabel('xi (kpc)')
# plt.ylabel('eta (kpc)')
# clb=plt.colorbar()
# clb.set_label('v (km/s)')
# plt.savefig('/Users/amandaquirk/Desktop/MWFG_MS_individual.png')
# plt.close()

# plt.figure(figsize=(4, 6))
# plt.scatter(MWFG_AGy_xi, MWFG_AGy_eta, s=1, c=MWFG_AGy_v, cmap='rainbow',vmin=-300,vmax=0)
# plt.gca().invert_xaxis()
# plt.xlabel('xi (kpc)')
# plt.ylabel('eta (kpc)')
# clb=plt.colorbar()
# clb.set_label('v (km/s)')
# plt.savefig('/Users/amandaquirk/Desktop/MWFG_AGy_individual.png')
# plt.close()

# plt.figure(figsize=(4, 6))
# plt.scatter(MWFG_AGo_xi, MWFG_AGo_eta, s=1, c=MWFG_AGo_v, cmap='rainbow',vmin=-300,vmax=0)
# plt.gca().invert_xaxis()
# plt.xlabel('xi (kpc)')
# plt.ylabel('eta (kpc)')
# clb=plt.colorbar()
# clb.set_label('v (km/s)')
# plt.savefig('/Users/amandaquirk/Desktop/MWFG_AGo_individual.png')
# plt.close()

# plt.figure(figsize=(4, 6))
# plt.scatter(MWFG_RG_xi, MWFG_RG_eta, s=1, c=MWFG_RG_v, cmap='rainbow',vmin=-300,vmax=0)
# plt.gca().invert_xaxis()
# plt.xlabel('xi (kpc)')
# plt.ylabel('eta (kpc)')
# clb=plt.colorbar()
# clb.set_label('v (km/s)')
# plt.savefig('/Users/amandaquirk/Desktop/MWFG_RG_individual.png')
# plt.close()

#smoothing the velocities-- not using errors
# MWFG_MS_smoothedv=[] #will contain the averaged velocity of the smoothed circle centered on each point, if a good center
# MWFG_MS_xi_goodcenter=[] #will contain x coordinate of star that is a good center
# MWFG_MS_eta_goodcenter=[] #will contain y coordinate of star that is a good center
# for i in range(len(MWFG_MS_xi)):
# 	#picks out the first point
# 	c1=SkyCoord(MWFG_MS_xi[i]/13.67, MWFG_MS_eta[i]/13.67, unit=(u.deg,u.deg)) 
# 	#will contain all of the velocities that are within the circle- the average of this will be added to the smoothed v array
# 	velocities=[]
# 	positions=[]
# 	for j in range(len(MWFG_MS_xi)):
# 		#below will go through all points to determine if the points are within the smoothing circle
# 		c2=SkyCoord(MWFG_MS_xi[j]/13.67, MWFG_MS_eta[j]/13.67, unit=(u.deg,u.deg))
# 		sep=c1.separation(c2)
# 		#below adds the velocity to the array containing the smoothed circle data
# 		if sep.arcsecond<200:
# 			velocities.append(MWFG_MS_v[j])
# 	#below adds the mean of all of the velocities to the points within the circle to the smoothed array
# 	#if star[i] has fewer than 15 neighbors, it cannot be used as a center but can be used as a neighbor
# 	if len(velocities)>15:
# 		MWFG_MS_smoothedv.append(np.mean(velocities))
# 		MWFG_MS_xi_goodcenter.append(MWFG_MS_xi[i])
# 		MWFG_MS_eta_goodcenter.append(MWFG_MS_eta[i])

# MWFG_AGy_smoothedv=[] #will contain the averaged velocity of the smoothed circle centered on each point, if a good center
# MWFG_AGy_xi_goodcenter=[] #will contain x coordinate of star that is a good center
# MWFG_AGy_eta_goodcenter=[] #will contain y coordinate of star that is a good center
# for i in range(len(MWFG_AGy_xi)):
# 	#picks out the first point
# 	c1=SkyCoord(MWFG_AGy_xi[i]/13.67, MWFG_AGy_eta[i]/13.67, unit=(u.deg,u.deg)) 
# 	#will contain all of the velocities that are within the circle- the average of this will be added to the smoothed v array
# 	velocities=[]
# 	positions=[]
# 	for j in range(len(MWFG_AGy_xi)):
# 		#below will go through all points to determine if the points are within the smoothing circle
# 		c2=SkyCoord(MWFG_AGy_xi[j]/13.67, MWFG_AGy_eta[j]/13.67, unit=(u.deg,u.deg))
# 		sep=c1.separation(c2)
# 		#below adds the velocity to the array containing the smoothed circle data
# 		if sep.arcsecond<200:
# 			velocities.append(MWFG_AGy_v[j])
# 	#below adds the mean of all of the velocities to the points within the circle to the smoothed array
# 	#if star[i] has fewer than 15 neighbors, it cannot be used as a center but can be used as a neighbor
# 	if len(velocities)>15:
# 		MWFG_AGy_smoothedv.append(np.mean(velocities))
# 		MWFG_AGy_xi_goodcenter.append(MWFG_AGy_xi[i])
# 		MWFG_AGy_eta_goodcenter.append(MWFG_AGy_eta[i])

# MWFG_AGo_smoothedv=[] #will contain the averaged velocity of the smoothed circle centered on each point, if a good center
# MWFG_AGo_xi_goodcenter=[] #will contain x coordinate of star that is a good center
# MWFG_AGo_eta_goodcenter=[] #will contain y coordinate of star that is a good center
# for i in range(len(MWFG_AGo_xi)):
# 	#picks out the first point
# 	c1=SkyCoord(MWFG_AGo_xi[i]/13.67, MWFG_AGo_eta[i]/13.67, unit=(u.deg,u.deg)) 
# 	#will contain all of the velocities that are within the circle- the average of this will be added to the smoothed v array
# 	velocities=[]
# 	positions=[]
# 	for j in range(len(MWFG_AGo_xi)):
# 		#below will go through all points to determine if the points are within the smoothing circle
# 		c2=SkyCoord(MWFG_AGo_xi[j]/13.67, MWFG_AGo_eta[j]/13.67, unit=(u.deg,u.deg))
# 		sep=c1.separation(c2)
# 		#below adds the velocity to the array containing the smoothed circle data
# 		if sep.arcsecond<200:
# 			velocities.append(MWFG_AGo_v[j])
# 	#below adds the mean of all of the velocities to the points within the circle to the smoothed array
# 	#if star[i] has fewer than 15 neighbors, it cannot be used as a center but can be used as a neighbor
# 	if len(velocities)>15:
# 		MWFG_AGo_smoothedv.append(np.mean(velocities))
# 		MWFG_AGo_xi_goodcenter.append(MWFG_AGo_xi[i])
# 		MWFG_AGo_eta_goodcenter.append(MWFG_AGo_eta[i])

# MWFG_RG_smoothedv=[] #will contain the averaged velocity of the smoothed circle centered on each point, if a good center
# MWFG_RG_xi_goodcenter=[] #will contain x coordinate of star that is a good center
# MWFG_RG_eta_goodcenter=[] #will contain y coordinate of star that is a good center
# for i in range(len(MWFG_RG_xi)):
# 	#picks out the first point
# 	c1=SkyCoord(MWFG_RG_xi[i]/13.67, MWFG_RG_eta[i]/13.67, unit=(u.deg,u.deg)) 
# 	#will contain all of the velocities that are within the circle- the average of this will be added to the smoothed v array
# 	velocities=[]
# 	positions=[]
# 	for j in range(len(MWFG_RG_xi)):
# 		#below will go through all points to determine if the points are within the smoothing circle
# 		c2=SkyCoord(MWFG_RG_xi[j]/13.67, MWFG_RG_eta[j]/13.67, unit=(u.deg,u.deg))
# 		sep=c1.separation(c2)
# 		#below adds the velocity to the array containing the smoothed circle data
# 		if sep.arcsecond<200:
# 			velocities.append(MWFG_RG_v[j])
# 	#below adds the mean of all of the velocities to the points within the circle to the smoothed array
# 	#if star[i] has fewer than 15 neighbors, it cannot be used as a center but can be used as a neighbor
# 	if len(velocities)>15:
# 		MWFG_RG_smoothedv.append(np.mean(velocities))
# 		MWFG_RG_xi_goodcenter.append(MWFG_RG_xi[i])
# 		MWFG_RG_eta_goodcenter.append(MWFG_RG_eta[i])

# #plotting
# plt.figure(figsize=(4, 6))
# plt.scatter(MWFG_MS_xi_goodcenter, MWFG_MS_eta_goodcenter, s=1, c=MWFG_MS_smoothedv, cmap='rainbow',vmin=-300,vmax=0)
# plt.gca().invert_xaxis()
# plt.xlabel('xi (kpc)')
# plt.ylabel('eta (kpc)')
# clb=plt.colorbar()
# clb.set_label('v (km/s)')
# plt.savefig('/Users/amandaquirk/Desktop/MWFG_MS_smoothed.png')
# plt.close()

# plt.figure(figsize=(4, 6))
# plt.scatter(MWFG_AGy_xi_goodcenter, MWFG_AGy_eta_goodcenter, s=1, c=MWFG_AGy_smoothedv, cmap='rainbow',vmin=-300,vmax=0)
# plt.gca().invert_xaxis()
# plt.xlabel('xi (kpc)')
# plt.ylabel('eta (kpc)')
# clb=plt.colorbar()
# clb.set_label('v (km/s)')
# plt.savefig('/Users/amandaquirk/Desktop/MWFG_AGy_smoothed.png')
# plt.close()

# plt.figure(figsize=(4, 6))
# plt.scatter(MWFG_AGo_xi_goodcenter, MWFG_AGo_eta_goodcenter, s=1, c=MWFG_AGo_smoothedv, cmap='rainbow',vmin=-300,vmax=0)
# plt.gca().invert_xaxis()
# plt.xlabel('xi (kpc)')
# plt.ylabel('eta (kpc)')
# clb=plt.colorbar()
# clb.set_label('v (km/s)')
# plt.savefig('/Users/amandaquirk/Desktop/MWFG_AGo_smoothed.png')
# plt.close()

# plt.figure(figsize=(4, 6))
# plt.scatter(MWFG_RG_xi_goodcenter, MWFG_RG_eta_goodcenter, s=1, c=MWFG_RG_smoothedv, cmap='rainbow',vmin=-300,vmax=0)
# plt.gca().invert_xaxis()
# plt.xlabel('xi (kpc)')
# plt.ylabel('eta (kpc)')
# clb=plt.colorbar()
# clb.set_label('v (km/s)')
# plt.savefig('/Users/amandaquirk/Desktop/MWFG_RG_smoothed.png')
# plt.close()

#rotation curve stuff
#shifted coordinates
def x(xi, eta):
	xi_deg=xi/13.67
	eta_deg=eta/13.67
	sine=np.sin(float(37*np.pi) / 180)
	cosine=np.cos(float(37*np.pi) / 180)
	x=(xi_deg*cosine)-(eta_deg*sine)
	return x 

def y(xi, eta):
	xi_deg=xi/13.67
	eta_deg=eta/13.67
	sine=np.sin(float(37*np.pi) / 180)
	cosine=np.cos(float(37*np.pi) / 180)
	y=(eta_deg*cosine)+(xi_deg*sine)
	return y

# MWFG_MS_x_smoothed=x(np.array(MWFG_MS_xi_goodcenter), np.array(MWFG_MS_eta_goodcenter))
# MWFG_MS_y_smoothed=y(np.array(MWFG_MS_xi_goodcenter), np.array(MWFG_MS_eta_goodcenter))

# MWFG_AGy_x_smoothed=x(np.array(MWFG_AGy_xi_goodcenter), np.array(MWFG_AGy_eta_goodcenter))
# MWFG_AGy_y_smoothed=y(np.array(MWFG_AGy_xi_goodcenter), np.array(MWFG_AGy_eta_goodcenter))

# MWFG_AGo_x_smoothed=x(np.array(MWFG_AGo_xi_goodcenter), np.array(MWFG_AGo_eta_goodcenter))
# MWFG_AGo_y_smoothed=y(np.array(MWFG_AGo_xi_goodcenter), np.array(MWFG_AGo_eta_goodcenter))

# MWFG_RG_x_smoothed=x(np.array(MWFG_RG_xi_goodcenter), np.array(MWFG_RG_eta_goodcenter))
# MWFG_RG_y_smoothed=y(np.array(MWFG_RG_xi_goodcenter), np.array(MWFG_RG_eta_goodcenter))

MWFG_MS_x_individual=x(np.array(MWFG_MS_xi), np.array(MWFG_MS_eta))
MWFG_MS_y_individual=y(np.array(MWFG_MS_xi), np.array(MWFG_MS_eta))

MWFG_AGy_x_individual=x(np.array(MWFG_AGy_xi), np.array(MWFG_AGy_eta))
MWFG_AGy_y_individual=y(np.array(MWFG_AGy_xi), np.array(MWFG_AGy_eta))

MWFG_AGo_x_individual=x(np.array(MWFG_AGo_xi), np.array(MWFG_AGo_eta))
MWFG_AGo_y_individual=y(np.array(MWFG_AGo_xi), np.array(MWFG_AGo_eta))

MWFG_RG_x_individual=x(np.array(MWFG_RG_xi), np.array(MWFG_RG_eta))
MWFG_RG_y_individual=y(np.array(MWFG_RG_xi), np.array(MWFG_RG_eta))

#deprojected distance
def distance(x, y):
	inclination_factor=np.cos(float(77*np.pi) / 180)**2
	ang_distance_sq=[(a**2)+(float(b**2)/inclination_factor) for a,b in zip(y,x)]
	ang_dist=[np.sqrt(a) for a in ang_distance_sq]
	dist=[a*13.67 for a in ang_dist]
	return dist

# MWFG_MS_r_smoothed=distance(MWFG_MS_x_smoothed, MWFG_MS_y_smoothed)
# MWFG_AGy_r_smoothed=distance(MWFG_AGy_x_smoothed, MWFG_AGy_y_smoothed)
# MWFG_AGo_r_smoothed=distance(MWFG_AGo_x_smoothed, MWFG_AGo_y_smoothed)
# MWFG_RG_r_smoothed=distance(MWFG_RG_x_smoothed, MWFG_RG_y_smoothed)

MWFG_MS_r_individual=distance(MWFG_MS_x_individual, MWFG_MS_y_individual)
MWFG_AGy_r_individual=distance(MWFG_AGy_x_individual, MWFG_AGy_y_individual)
MWFG_AGo_r_individual=distance(MWFG_AGo_x_individual, MWFG_AGo_y_individual)
MWFG_RG_r_individual=distance(MWFG_RG_x_individual, MWFG_RG_y_individual)

#the function below will assign each star and gas cloud a position angle using their xi and eta coordinates
def PA(x,y): #defined as west of north as positioned in maps (inverted x axis)
	if x>0:
		rad=np.arctan(float(y)/x)
		deg=90-(float(rad*180)/np.pi)
	else:
		rad=np.arctan(float(y)/x)
		deg=270-(float(rad*180)/np.pi)
	return deg + 37

# MWFG_MS_PA_smoothed=[]
# for i in range(len(MWFG_MS_r_smoothed)):
# 	MWFG_MS_PA_smoothed.append(PA(MWFG_MS_x_smoothed[i], MWFG_MS_y_smoothed[i]))
# MWFG_AGy_PA_smoothed=[]
# for i in range(len(MWFG_AGy_r_smoothed)):
# 	MWFG_AGy_PA_smoothed.append(PA(MWFG_AGy_x_smoothed[i], MWFG_AGy_y_smoothed[i]))
# MWFG_AGo_PA_smoothed=[]
# for i in range(len(MWFG_AGo_r_smoothed)):
# 	MWFG_AGo_PA_smoothed.append(PA(MWFG_AGo_x_smoothed[i], MWFG_AGo_y_smoothed[i]))
# MWFG_RG_PA_smoothed=[]
# for i in range(len(MWFG_RG_r_smoothed)):
# 	MWFG_RG_PA_smoothed.append(PA(MWFG_RG_x_smoothed[i], MWFG_RG_y_smoothed[i]))

MWFG_MS_PA_individual=[]
for i in range(len(MWFG_MS_r_individual)):
	MWFG_MS_PA_individual.append(PA(MWFG_MS_x_individual[i], MWFG_MS_y_individual[i]))
MWFG_AGy_PA_individual=[]
for i in range(len(MWFG_AGy_r_individual)):
	MWFG_AGy_PA_individual.append(PA(MWFG_AGy_x_individual[i], MWFG_AGy_y_individual[i]))
MWFG_AGo_PA_individual=[]
for i in range(len(MWFG_AGo_r_individual)):
	MWFG_AGo_PA_individual.append(PA(MWFG_AGo_x_individual[i], MWFG_AGo_y_individual[i]))
MWFG_RG_PA_individual=[]
for i in range(len(MWFG_RG_r_individual)):
	MWFG_RG_PA_individual.append(PA(MWFG_RG_x_individual[i], MWFG_RG_y_individual[i]))

#now we need to assign each star a PA_ring and i_ring using HI data. we define these values by seeing which HI ring a star is cloest to and assigning it the corresponding PA_ring and i_ring value (HI_PA and HI_i respectively). we also will pair a HI velocity to each star to later calculate the asymmetric drift
#reads in the HI data
HI_r, HI_PA, HI_i, HI_v=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/HI_PA_i_vrot.txt', unpack=True)

#below defines a function to determine which HI ring a star is cloest to
def find_nearest_ring(radius):
	idx=[]
	for i in range(len(HI_r)):
		idx.append(np.abs(HI_r[i]-radius))
		n=np.argmin(idx)
	return n #index of the radius closest to the star's radius

#assing a PA  and i to each star and HI component
# MWFG_MS_PA_ring_smoothed=[]
# MWFG_MS_i_smoothed=[]
# for j in range(len(MWFG_MS_r_smoothed)):
# 	N=find_nearest_ring(MWFG_MS_r_smoothed[j])
# 	MWFG_MS_PA_ring_smoothed.append(HI_PA[N])
# 	MWFG_MS_i_smoothed.append(HI_i[N])

# MWFG_AGy_PA_ring_smoothed=[]
# MWFG_AGy_i_smoothed=[]
# for j in range(len(MWFG_AGy_r_smoothed)):
# 	N=find_nearest_ring(MWFG_AGy_r_smoothed[j])
# 	MWFG_AGy_PA_ring_smoothed.append(HI_PA[N])
# 	MWFG_AGy_i_smoothed.append(HI_i[N])

# MWFG_AGo_PA_ring_smoothed=[]
# MWFG_AGo_i_smoothed=[]
# for j in range(len(MWFG_AGo_r_smoothed)):
# 	N=find_nearest_ring(MWFG_AGo_r_smoothed[j])
# 	MWFG_AGo_PA_ring_smoothed.append(HI_PA[N])
# 	MWFG_AGo_i_smoothed.append(HI_i[N])

# MWFG_RG_PA_ring_smoothed=[]
# MWFG_RG_i_smoothed=[]
# for j in range(len(MWFG_RG_r_smoothed)):
# 	N=find_nearest_ring(MWFG_RG_r_smoothed[j])
# 	MWFG_RG_PA_ring_smoothed.append(HI_PA[N])
# 	MWFG_RG_i_smoothed.append(HI_i[N])

MWFG_MS_PA_ring_individual=[]
MWFG_MS_i_individual=[]
for j in range(len(MWFG_MS_r_individual)):
	N=find_nearest_ring(MWFG_MS_r_individual[j])
	MWFG_MS_PA_ring_individual.append(HI_PA[N])
	MWFG_MS_i_individual.append(HI_i[N])

MWFG_AGy_PA_ring_individual=[]
MWFG_AGy_i_individual=[]
for j in range(len(MWFG_AGy_r_individual)):
	N=find_nearest_ring(MWFG_AGy_r_individual[j])
	MWFG_AGy_PA_ring_individual.append(HI_PA[N])
	MWFG_AGy_i_individual.append(HI_i[N])

MWFG_AGo_PA_ring_individual=[]
MWFG_AGo_i_individual=[]
for j in range(len(MWFG_AGo_r_individual)):
	N=find_nearest_ring(MWFG_AGo_r_individual[j])
	MWFG_AGo_PA_ring_individual.append(HI_PA[N])
	MWFG_AGo_i_individual.append(HI_i[N])

MWFG_RG_PA_ring_individual=[]
MWFG_RG_i_individual=[]
for j in range(len(MWFG_RG_r_individual)):
	N=find_nearest_ring(MWFG_RG_r_individual[j])
	MWFG_RG_PA_ring_individual.append(HI_PA[N])
	MWFG_RG_i_individual.append(HI_i[N])

#tilted ring rotational velocity
def Vrot_tilted_ring(v,PA_ring,PA_star, i_ring): 
	vsys= -300 #km/s, as defined in Claire's thesis
	A=[float(a-vsys)/np.sin(float(b*np.pi) / 180) for a, b in zip(v, i_ring)]
	B= [(np.tan(float((a-b)*np.pi) / 180))**2 for a,b in zip(PA_ring, PA_star)]
	C= [(np.cos(float(a*np.pi) / 180))**2 for a in i_ring]
	rotation_velocity= [a*np.sqrt(1+(float(b)/c)) for a,b,c in zip(A,B,C)]
	positive=[np.absolute(a) for a in rotation_velocity]
	return positive

# MWFG_MS_vrot_tilt_smoothed=Vrot_tilted_ring(MWFG_MS_smoothedv,MWFG_MS_PA_ring_smoothed,MWFG_MS_PA_smoothed, MWFG_MS_i_smoothed)
# MWFG_AGy_vrot_tilt_smoothed=Vrot_tilted_ring(MWFG_AGy_smoothedv,MWFG_AGy_PA_ring_smoothed,MWFG_AGy_PA_smoothed, MWFG_AGy_i_smoothed)
# MWFG_AGo_vrot_tilt_smoothed=Vrot_tilted_ring(MWFG_AGo_smoothedv,MWFG_AGo_PA_ring_smoothed,MWFG_AGo_PA_smoothed, MWFG_AGo_i_smoothed)
# MWFG_RG_vrot_tilt_smoothed=Vrot_tilted_ring(MWFG_RG_smoothedv,MWFG_RG_PA_ring_smoothed,MWFG_RG_PA_smoothed, MWFG_RG_i_smoothed)

MWFG_MS_vrot_tilt_individual=Vrot_tilted_ring(MWFG_MS_v,MWFG_MS_PA_ring_individual,MWFG_MS_PA_individual, MWFG_MS_i_individual)
MWFG_AGy_vrot_tilt_individual=Vrot_tilted_ring(MWFG_AGy_v,MWFG_AGy_PA_ring_individual,MWFG_AGy_PA_individual, MWFG_AGy_i_individual)
MWFG_AGo_vrot_tilt_individual=Vrot_tilted_ring(MWFG_AGo_v,MWFG_AGo_PA_ring_individual,MWFG_AGo_PA_individual, MWFG_AGo_i_individual)
MWFG_RG_vrot_tilt_individual=Vrot_tilted_ring(MWFG_RG_v,MWFG_RG_PA_ring_individual,MWFG_RG_PA_individual, MWFG_RG_i_individual)

#plotting
# plt.scatter(MWFG_MS_r_smoothed, MWFG_MS_vrot_tilt_smoothed, s=1, c='b')
# plt.xlim(4,20)
# plt.ylim(10,300)
# plt.xlabel('r [kpc]')
# plt.ylabel('rotational v [km/s]')
# plt.savefig('/Users/amandaquirk/Desktop/MWFG_MS_vrot_tilt_smoothed.png')
# plt.close()

# plt.scatter(MWFG_AGy_r_smoothed, MWFG_AGy_vrot_tilt_smoothed, s=1, c='m')
# plt.xlim(4,20)
# plt.ylim(10,300)
# plt.xlabel('r [kpc]')
# plt.ylabel('rotational v [km/s]')
# plt.savefig('/Users/amandaquirk/Desktop/MWFG_AGy_vrot_tilt_smoothed.png')
# plt.close()

# plt.scatter(MWFG_AGo_r_smoothed, MWFG_AGo_vrot_tilt_smoothed, s=1, c='m')
# plt.xlim(4,20)
# plt.ylim(10,300)
# plt.xlabel('r [kpc]')
# plt.ylabel('rotational v [km/s]')
# plt.savefig('/Users/amandaquirk/Desktop/MWFG_AGo_vrot_tilt_smoothed.png')
# plt.close()

# plt.scatter(MWFG_RG_r_smoothed, MWFG_RG_vrot_tilt_smoothed, s=1, c='r')
# plt.xlim(4,20)
# plt.ylim(10,300)
# plt.xlabel('r [kpc]')
# plt.ylabel('rotational v [km/s]')
# plt.savefig('/Users/amandaquirk/Desktop/MWFG_RG_vrot_tilt_smoothed.png')
# plt.close()

from matplotlib import rc 

rc('font', family = 'serif')
f, axes= plt.subplots(4,1, sharey=False, sharex=True, figsize=(4,9.5))
axes[0].scatter(MWFG_MS_r_individual, MWFG_MS_vrot_tilt_individual, s=4, c='b', alpha=0.4)
axes[1].scatter(MWFG_AGy_r_individual, MWFG_AGy_vrot_tilt_individual, s=2, c='m', alpha=0.4)
axes[2].scatter(MWFG_AGo_r_individual, MWFG_AGo_vrot_tilt_individual, s=2, c='k', alpha=0.4)
axes[3].scatter(MWFG_RG_r_individual, MWFG_RG_vrot_tilt_individual, s=2, c='r', alpha=0.4)
axes[0].annotate('MS', xy=(19,115), horizontalalignment='right', fontsize=12)
axes[1].annotate('young AGB', xy=(19,115), horizontalalignment='right', fontsize=12)
axes[2].annotate('older AGB', xy=(19,115), horizontalalignment='right', fontsize=12)
axes[3].annotate('RGB', xy=(19,115), horizontalalignment='right', fontsize=12)

for ax in axes:
	ax.set_xlim(4, 20)
	#ax.set_ylabel(r'$\rm v_{\rm rot} \ (km\ s^{-1})$', fontsize=13)
	ax.set_ylim(100,300)
	ax.tick_params(axis='x',which='both',bottom='on',top='off', direction='out')
	ax.tick_params(axis='x',which='both',top='on', direction='in')
	ax.tick_params(axis='y',which='both',left='on',top='off', direction='out')
	ax.tick_params(axis='y',which='both',right='on', direction='in')
	ax.tick_params(which='both', width=2)
	ax.tick_params(which='major', length=7)
	ax.tick_params(which='minor', length=4)
	ax.tick_params(labelsize=12) 
	ax.minorticks_on()
	for axis in ['top','bottom','left','right']:
	        ax.spines[axis].set_linewidth(2)
axes[3].set_xlabel(r'$\rm Radial\ Distance:\ \itr \ \rm(kpc)$', fontsize=13)
nbins = len(axes[0].get_yticklabels())-1
axes[3].yaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='upper'))
axes[1].yaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='upper'))
axes[2].yaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='upper'))
f.subplots_adjust(left=0.17)
f.text(0.008, 0.5, r'$\rm Rotation\ Velocity:\ \itv_{\rm rot}\ \rm(km\ s^{-1})$', va='center', rotation='vertical', fontsize=13)
plt.subplots_adjust(wspace=0, hspace=0)
plt.savefig('/Users/amandaquirk/Desktop/MFWG_rotation_curves.pdf', bbox_inches='tight') 

# plt.scatter(MWFG_MS_r_individual, MWFG_MS_vrot_tilt_individual, s=8, c='b')
# plt.xlim(4,20)
# plt.ylim(10,300)
# plt.xlabel('r [kpc]')
# plt.ylabel('rotational v [km/s]')
# plt.savefig('/Users/amandaquirk/Desktop/MWFG_MS_vrot_tilt_individual.png')
# plt.close()

# plt.scatter(MWFG_AGy_r_individual, MWFG_AGy_vrot_tilt_individual, s=1, c='m')
# plt.xlim(4,20)
# plt.ylim(10,300)
# plt.xlabel('r [kpc]')
# plt.ylabel('rotational v [km/s]')
# plt.savefig('/Users/amandaquirk/Desktop/MWFG_AGy_vrot_tilt_individual.png')
# plt.close()

# plt.scatter(MWFG_AGo_r_individual, MWFG_AGo_vrot_tilt_individual, s=1, c='m')
# plt.xlim(4,20)
# plt.ylim(10,300)
# plt.xlabel('r [kpc]')
# plt.ylabel('rotational v [km/s]')
# plt.savefig('/Users/amandaquirk/Desktop/MWFG_AGo_vrot_tilt_individual.png')
# plt.close()

# plt.scatter(MWFG_RG_r_individual, MWFG_RG_vrot_tilt_individual, s=1, c='r')
# plt.xlim(4,20)
# plt.ylim(10,300)
# plt.xlabel('r [kpc]')
# plt.ylabel('rotational v [km/s]')
# plt.savefig('/Users/amandaquirk/Desktop/MWFG_RG_vrot_tilt_individual.png')
# plt.close()

#files
# file=open('/Users/amandaquirk/Desktop/MWFG_MS_vrot_tilt_smoothed.txt', 'w')
# file.write('#r (kpc), v_rot tr model\n')
# for i in range(len(MWFG_MS_vrot_tilt_smoothed)):
# 	file.write('{} {}\n'.format(MWFG_MS_r_smoothed[i],MWFG_MS_vrot_tilt_smoothed[i]))
# file.close()

# file=open('/Users/amandaquirk/Desktop/MWFG_AGy_vrot_tilt_smoothed.txt', 'w')
# file.write('#r (kpc), v_rot tr model\n')
# for i in range(len(MWFG_AGy_vrot_tilt_smoothed)):
# 	file.write('{} {}\n'.format(MWFG_AGy_r_smoothed[i],MWFG_AGy_vrot_tilt_smoothed[i]))
# file.close()

# file=open('/Users/amandaquirk/Desktop/MWFG_AGo_vrot_tilt_smoothed.txt', 'w')
# file.write('#r (kpc), v_rot tr model\n')
# for i in range(len(MWFG_AGo_vrot_tilt_smoothed)):
# 	file.write('{} {}\n'.format(MWFG_AGo_r_smoothed[i],MWFG_AGo_vrot_tilt_smoothed[i]))
# file.close()

# file=open('/Users/amandaquirk/Desktop/MWFG_RG_vrot_tilt_smoothed.txt', 'w')
# file.write('#r (kpc), v_rot tr model\n')
# for i in range(len(MWFG_RG_vrot_tilt_smoothed)):
# 	file.write('{} {}\n'.format(MWFG_RG_r_smoothed[i],MWFG_RG_vrot_tilt_smoothed[i]))
# file.close()

# file=open('/Users/amandaquirk/Desktop/MWFG_MS_vrot_tilt_individual.txt', 'w')
# file.write('#r (kpc), v_rot tr model\n')
# for i in range(len(MWFG_MS_vrot_tilt_individual)):
# 	file.write('{} {}\n'.format(MWFG_MS_r_individual[i],MWFG_MS_vrot_tilt_individual[i]))
# file.close()

# file=open('/Users/amandaquirk/Desktop/MWFG_AGy_vrot_tilt_individual.txt', 'w')
# file.write('#r (kpc), v_rot tr model\n')
# for i in range(len(MWFG_AGy_vrot_tilt_individual)):
# 	file.write('{} {}\n'.format(MWFG_AGy_r_individual[i],MWFG_AGy_vrot_tilt_individual[i]))
# file.close()

# file=open('/Users/amandaquirk/Desktop/MWFG_AGo_vrot_tilt_individual.txt', 'w')
# file.write('#r (kpc), v_rot tr model\n')
# for i in range(len(MWFG_AGo_vrot_tilt_individual)):
# 	file.write('{} {}\n'.format(MWFG_AGo_r_individual[i],MWFG_AGo_vrot_tilt_individual[i]))
# file.close()

# file=open('/Users/amandaquirk/Desktop/MWFG_RG_vrot_tilt_individual.txt', 'w')
# file.write('#r (kpc), v_rot tr model\n')
# for i in range(len(MWFG_RG_vrot_tilt_individual)):
# 	file.write('{} {}\n'.format(MWFG_RG_r_individual[i],MWFG_RG_vrot_tilt_individual[i]))
# file.close()
