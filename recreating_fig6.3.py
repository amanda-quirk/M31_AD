from astropy.io import fits 
from astropy.table import Table 
from astropy import units as u
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt 
import numpy as np 
from shapely.geometry import Point 
from shapely.geometry.polygon import Polygon

#import the data file
hdulist= fits.open('/Users/amandaquirk/Documents/AsymmetricDrift/Data/subMasterSPLASH.fits', memmap=True)

#puts data into an accessible table
data=Table(hdulist[1].data)

#creates arrays of the parameters than I need
zqual=data['ZQUAL'] #quality of data point
L=data['LIKELIHOOD'] #if negative, MWFG star
RA=data['RA']
Dec=data['DEC']
F814W=data['F814W']
F475W=data['F475W']
z=data['Z'] #redshift

#these arrays will contain data from stars with good zqual values and that are not MWFGstars
z_good=[]
RA_good=[]
Dec_good=[]
F8W=[]
F4W=[]
for i in range(len(zqual)):
	#below selects the candidate stars
        if zqual[i]>=3 and L[i]>0 or np.isnan(L[i])==True:
                z_good.append(z[i])
                RA_good.append(RA[i])
                Dec_good.append(Dec[i])
                F8W.append(F814W[i])
                F4W.append(F475W[i])

#the next step is to separate the stars into age groups: MS, RGB, young AGB, and old AGB based on a color magnitude diagram. first i create a CDM

#below will contain the color of every good star
color=[]
for i in range(len(F8W)):
	color.append(F4W[i]-F8W[i])

#print(len(color))

#below creates a CMD- this is just a visual check
#plt.scatter(color, F8W, color='b', alpha=0.4)
#plt.xlabel('F475W - F814W (mag)')
#plt.xlim(-1,8)
#plt.ylim(18,23)
#plt.gca().invert_yaxis()
#plt.ylabel('F814W (mag)')
#plt.savefig('/Users/amandaquirk/Desktop/CDM_M31.png')

#separate into age groups based on color and magnitude. 

#below will contain position and redshift values for MS stars
MS_z=[]
MS_RA=[]
MS_Dec=[]
MS_color=[]
MS_F8W=[]
for i in range(len(color)):
	if color[i]<=1:
                MS_z.append(z_good[i])
                MS_RA.append(RA_good[i])
                MS_Dec.append(Dec_good[i])
                #MS_color.append(color[i])
                #MS_F8W.append(F8W[i])

#print(len(MS_z), len(MS_RA), len(MS_Dec))

#below tests to see if the star is in the RGB
RG_z=[]
RG_RA=[]
RG_Dec=[]
RG_color=[]
RG_F8W=[]
#polygon creates a shape that corresponds to the area of the CMD that contains RGB stars
polygon= Polygon([(2,23), (2.7, 20.4), (7.5,23)])
for i in range(len(color)):
        point= Point(color[i], F8W[i])
        if polygon.contains(point)==True:
                RG_z.append(z_good[i])
                RG_RA.append(RA_good[i])
                RG_Dec.append(Dec_good[i])
                #RG_color.append(color[i])
                #RG_F8W.append(F8W[i])

#below tests to see if the star is in the old AGB
AGo_z=[]
AGo_RA=[]
AGo_Dec=[]
AGo_color=[]
AGo_F8W=[]
#polygon creates a shape that corresponds to the area of the CMD that contains old AGB stars
polygon= Polygon([(2.7,20.4), (7.5, 23), (8, 20.5)])
for i in range(len(color)):
        point= Point(color[i], F8W[i])
        if polygon.contains(point)==True:
                AGo_z.append(z_good[i])
                AGo_RA.append(RA_good[i])
                AGo_Dec.append(Dec_good[i])
                #AGo_color.append(color[i])
                #AGo_F8W.append(F8W[i])

#below tests to see if the star is in the young AGB
AGy_z=[]
AGy_RA=[]
AGy_Dec=[]
AGy_color=[]
AGy_F8W=[]
#polygon creates a shape that corresponds to the area of the CMD that contains young AGB stars
polygon= Polygon([(3.5,18), (8,18), (2.7, 20.4), (3.5,20.5), (8,20.5)])
for i in range(len(color)):
        point= Point(color[i], F8W[i])
        if polygon.contains(point)==True:
                AGy_z.append(z_good[i])
                AGy_RA.append(RA_good[i])
                AGy_Dec.append(Dec_good[i])
                #AGy_color.append(color[i])
                #AGy_F8W.append(F8W[i])

#below does a visual check that the age group appears on the correct CMD
#plt.scatter(AGy_color, AGy_F8W, color='b', alpha=0.4)
#plt.xlabel('F475W - F814W (mag)')
#plt.xlim(-1,8)
#plt.ylim(18,23)
#plt.gca().invert_yaxis()
#plt.ylabel('F814W (mag)')
#plt.savefig('/Users/amandaquirk/Desktop/CDM_AGB_young.png')

#now I will convert the RA and Dec of each star candidate into eta and xi that are centered on M31
MS_xi=[]
MS_eta=[]

m31 = SkyCoord(10.6847083*u.deg, 41.26875*u.deg, frame='icrs')
for i in range(len(MS_RA)):
	c=SkyCoord(MS_RA[i], MS_Dec[i], frame='icrs', unit=(u.hourangle,u.deg))
	#below puts the point in reference to M31
	c_inm31=c.transform_to(m31.skyoffset_frame())
	xi, eta=c_inm31.lon, c_inm31.lat
	MS_xi.append(xi.degree)
	MS_eta.append(eta.degree)

#below converts redshift into line of sight velocity (km/s)
MS_v=[]
for i in range(len(MS_RA)):
	MS_v.append(MS_z[i]*(3*10**5))

#below writes data to a file
#file=open('/Users/amandaquirk/Documents/AsymmetricDrift/Data/MS_individualv.txt', 'w')
#file.write('#eta (deg), xi (deg), v(km/s)\n')
#for i in range(len(MS_v)):
#	file.write('{} {} {}\n'.format(MS_eta[i],MS_xi[i], MS_v[i]))
#file.close()

#below is a visual test that the correct stars were chosen
#plt.scatter(MS_xi, MS_eta, s=1, c=MS_v)
#plt.xlim(-.5,1)
#plt.gca().invert_xaxis()
#plt.ylim(-1,1.2)
#plt.xlabel('xi (deg)')
#plt.ylabel('eta (deg)')
#clb=plt.colorbar()
#clb.set_label('v (km/s)')
#plt.savefig('/Users/amandaquirk/Desktop/MS_map1.png')
#plt.close()

RG_xi=[]
RG_eta=[]

m31 = SkyCoord(10.6847083*u.deg, 41.26875*u.deg)#, frame='icrs')
for i in range(len(RG_RA)):
        c=SkyCoord(RG_RA[i], RG_Dec[i], unit=(u.hourangle,u.deg))#, frame='icrs', unit=(u.hourangle,u.deg))
        #below puts the point in reference to M31
        c_inm31=c.transform_to(m31.skyoffset_frame())
        xi, eta=c_inm31.lon, c_inm31.lat
        RG_xi.append(xi.degree)
        RG_eta.append(eta.degree)

#print(min(RG_xi), max(RG_xi), min(RG_eta), max(RG_eta))

#below converts redshift into line of sight velocity (km/s)
RG_v=[]
for i in range(len(RG_RA)):
        RG_v.append(RG_z[i]*(3*10**5))

#below writes data to a file
#file=open('/Users/amandaquirk/Documents/AsymmetricDrift/Data/RG_individualv.txt', 'w')
#file.write('#eta (deg), xi (deg), v(km/s)\n')
#for i in range(len(RG_v)):
#	file.write('{} {} {}\n'.format(RG_eta[i],RG_xi[i], RG_v[i]))
#file.close()

#below is a visual test that the correct stars were chosen
#plt.scatter(RG_xi, RG_eta, s=1, c=RG_v)
#plt.xlim(-.5,1)
#plt.gca().invert_xaxis()
#plt.ylim(-1,1.2)
#plt.xlabel('xi (deg)')
#plt.ylabel('eta (deg)')
#clb=plt.colorbar()
#clb.set_label('v (km/s)')
#plt.savefig('/Users/amandaquirk/Desktop/RG_map1.png')
#plt.close()

AGo_xi=[]
AGo_eta=[]

m31 = SkyCoord(10.6847083*u.deg, 41.26875*u.deg, frame='icrs')
for i in range(len(AGo_RA)):
        c=SkyCoord(AGo_RA[i], AGo_Dec[i], frame='icrs', unit=(u.hourangle,u.deg))
        #below puts the point in reference to M31
        c_inm31=c.transform_to(m31.skyoffset_frame())
        xi, eta=c_inm31.lon, c_inm31.lat
        AGo_xi.append(xi.degree)
        AGo_eta.append(eta.degree)

#below converts redshift into line of sight velocity (km/s)
AGo_v=[]
for i in range(len(AGo_RA)):
        AGo_v.append(AGo_z[i]*(3*10**5))

#below writes data to a file
#file=open('/Users/amandaquirk/Documents/AsymmetricDrift/Data/AGo_individualv.txt', 'w')
#file.write('#eta (deg), xi (deg), v(km/s)\n')
#for i in range(len(AGo_v)):
#	file.write('{} {} {}\n'.format(AGo_eta[i],AGo_xi[i], AGo_v[i]))
#file.close()

#below is a visual test that the correct stars were chosen
#plt.scatter(AGo_xi, AGo_eta, s=1, c=AGo_v)
#plt.xlim(-.5,1)
#plt.gca().invert_xaxis()
#plt.ylim(-1,1.2)
#plt.xlabel('xi (deg)')
#plt.ylabel('eta (deg)')
#clb=plt.colorbar()
#clb.set_label('v (km/s)')
#plt.savefig('/Users/amandaquirk/Desktop/AGo_map1.png')
#plt.close()


AGy_v=[]
AGy_RAgood=[]
AGy_Decgood=[]
for i in range(len(AGy_RA)):
        #below converts redshift into line of sight velocity (km/s)
	vel=(AGy_z[i]*(3*10**5))
	if vel<200: #eliminates outlier velocity-- this is done completely arbitrarily
                AGy_v.append(vel)
                AGy_RAgood.append(AGy_RA[i])
                AGy_Decgood.append(AGy_Dec[i])

AGy_xi=[]
AGy_eta=[]

m31 = SkyCoord(10.6847083*u.deg, 41.26875*u.deg, frame='icrs')
for i in range(len(AGy_RAgood)):
        c=SkyCoord(AGy_RAgood[i], AGy_Decgood[i], frame='icrs', unit=(u.hourangle,u.deg))
        #below puts the point in reference to M31
        c_inm31=c.transform_to(m31.skyoffset_frame())
        xi, eta=c_inm31.lon, c_inm31.lat
        AGy_xi.append(xi.degree)
        AGy_eta.append(eta.degree)

#below writes data to a file
#file=open('/Users/amandaquirk/Documents/AsymmetricDrift/Data/AGy_individualv.txt', 'w')
#file.write('#eta (deg), xi (deg), v(km/s)\n')
#for i in range(len(AGy_v)):
#	file.write('{} {} {}\n'.format(AGy_eta[i],AGy_xi[i], AGy_v[i]))
#file.close()

#below is a visual test that the correct stars were chosen	
#plt.scatter(AGy_xi, AGy_eta, s=1, c=AGy_v)
#plt.xlim(-.5,1)
#plt.gca().invert_xaxis()
#plt.ylim(-1,1.2)
#plt.xlabel('xi (deg)')
#plt.ylabel('eta (deg)')
#clb=plt.colorbar()
#clb.set_label('v (km/s)')
#plt.savefig('/Users/amandaquirk/Desktop/AGy_map1.png')
#plt.close()

#i used below to examine the outliers in the AGy velocities
#plt.hist(AGy_v)
#print(len(AGy_v))
#plt.show()
#print(min(AGy_v), max(AGy_v), np.mean(AGy_v))	

#now we need need to create a plot showing the average velocity. to do this, i will use a smooth circle to create an average velocity for each new point- MS
MS_smoothedv=[] #will contain the averaged velocity of the smoothed circle centered on each point, if a good center
MS_xi_goodcenter=[] #will contain x coordinate of star that is a good center
MS_eta_goodcenter=[] #will contain y coordinate of star that is a good center
MS_dispersion=[] #will contain the LOS dipersion
for i in range(len(MS_xi)):
	#picks out the first point
	c1=SkyCoord(MS_xi[i], MS_eta[i], unit=(u.deg,u.deg)) 
	#will contain all of the velocities that are within the circle- the average of this will be added to the smoothed v array
	velocities=[]
	positions=[]
	for j in range(len(MS_xi)):
		#below will go through all points to determine if the points are within the smoothing circle
		c2=SkyCoord(MS_xi[j], MS_eta[j], unit=(u.deg,u.deg))
		sep=c1.separation(c2)
		#below adds the velocity to the array containing the smoothed circle data
		if sep.arcsecond<200:
			velocities.append(MS_v[j])
	#below adds the mean of all of the velocities to the points within the circle to the smoothed array
	#if star[i] has fewer than 15 neighbors, it cannot be used as a center but can be used as a neighbor
	if len(velocities)>15:
		MS_smoothedv.append(np.mean(velocities))
		MS_dispersion.append(np.std(velocities))
		MS_xi_goodcenter.append(MS_xi[i])
		MS_eta_goodcenter.append(MS_eta[i])

#write dispersions to a file
# file=open('/Users/amandaquirk/Documents/AsymmetricDrift/Data/MS_dispersion.txt', 'w')
# file.write('#sigma_LOS (km/s)\n')
# for i in range(len(MS_xi_goodcenter)):
# 	file.write('{}\n'.format(MS_dispersion[i]))
# file.close()

#below writes data to a file
#file=open('/Users/amandaquirk/Documents/AsymmetricDrift/Data/MS_smoothed_goodcenters.txt', 'w')
#file.write('#eta (deg), xi (deg), average v(km/s)\n')
#for i in range(len(MS_xi_goodcenter)):
#	file.write('{} {} {}\n'.format(MS_eta_goodcenter[i],MS_xi_goodcenter[i], MS_smoothedv[i]))
#file.close()

#below creates a 3D plot showing the physical location of the stars and their LOS velocities
#plt.scatter(MS_xi_goodcenter, MS_eta_goodcenter, s=1, c=MS_smoothedv, cmap='rainbow',vmin=-300,vmax=0)
#plt.xlim(-.5,1)
#plt.gca().invert_xaxis()
#plt.ylim(-1,1.2)
#plt.xlabel('xi (deg)')
#plt.ylabel('eta (deg)')
#clb=plt.colorbar()
#clb.set_label('v (km/s)')
#plt.savefig('/Users/amandaquirk/Documents/AsymmetricDrift/Data/Plots/MS_avg_map_fixed.png')
#plt.close()

#RGB
RG_smoothedv=[] #will contain the averaged velocity of the smoothed circle centered on each point, if a good center
RG_xi_goodcenter=[] #will contain x coordinate of star that is a good center
RG_eta_goodcenter=[] #will contain y coordinate of star that is a good center
RG_dispersion=[]

for i in range(len(RG_xi)):
        #picks out the first point
	c1=SkyCoord(RG_xi[i], RG_eta[i], unit=(u.deg,u.deg)) 
        #will contain all of the velocities that are within the circle- the average of this will be added to the smoothed v array
	velocities=[]
	positions=[]
	for j in range(len(RG_xi)):
                #below will go through all points to determine if the points are within the smoothing circle
		c2=SkyCoord(RG_xi[j], RG_eta[j], unit=(u.deg,u.deg)) 
		sep=c1.separation(c2)
                #below adds the velocity to the array containing the smoothed circle data
		if sep.arcsecond<200:
			velocities.append(RG_v[j])
        #below adds the mean of all of the velocities to the points within the circle to the smoothed array
        #if star[i] has fewer than 15 neighbors, it cannot be used as a center but can be used as a neighbor
	if len(velocities)>15:
		RG_smoothedv.append(np.mean(velocities))
		RG_xi_goodcenter.append(RG_xi[i])
		RG_eta_goodcenter.append(RG_eta[i])
		RG_dispersion.append(np.std(velocities))

#write dispersions to a file
file=open('/Users/amandaquirk/Documents/AsymmetricDrift/Data/RG_dispersion.txt', 'w')
file.write('#sigma_LOS (km/s)\n')
for i in range(len(RG_xi_goodcenter)):
	file.write('{}\n'.format(RG_dispersion[i]))
file.close()

#below writes data to a file
#file=open('/Users/amandaquirk/Documents/AsymmetricDrift/Data/RG_smoothed_goodcenters.txt', 'w')
#file.write('#eta (deg), xi (deg), average v(km/s)\n')
#for i in range(len(RG_xi_goodcenter)):
#	file.write('{} {} {}\n'.format(RG_eta_goodcenter[i],RG_xi_goodcenter[i], RG_smoothedv[i]))
#file.close()

#below creates a 3D plot showing the physical location of the stars and their LOS velocities
#plt.scatter(RG_xi_goodcenter, RG_eta_goodcenter, s=1, c=RG_smoothedv, cmap='rainbow',vmin=-300,vmax=0)
#plt.xlim(-.5,1)
#plt.gca().invert_xaxis()
#plt.ylim(-1,1.2)
#plt.xlabel('xi (deg)')
#plt.ylabel('eta (deg)')
#clb=plt.colorbar()
#clb.set_label('v (km/s)')
#plt.savefig('/Users/amandaquirk/Documents/AsymmetricDrift/Data/Plots/RG_avg_map_fixed.png')
#plt.close()

#AGB old
AGo_smoothedv=[] #will contain the averaged velocity of the smoothed circle centered on each point, if a good center
AGo_xi_goodcenter=[] #will contain x coordinate of star that is a good center
AGo_eta_goodcenter=[] #will contain y coordinate of star that is a good center
AGo_dispersion=[]
for i in range(len(AGo_xi)):
        #picks out the first point
	c1=SkyCoord(AGo_xi[i], AGo_eta[i], unit=(u.deg,u.deg)) 
        #will contain all of the velocities that are within the circle- the average of this will be added to the smoothed v array
	velocities=[]
	positions=[]
	for j in range(len(AGo_xi)):
                #below will go through all points to determine if the points are within the smoothing circle
		c2=SkyCoord(AGo_xi[j], AGo_eta[j], unit=(u.deg,u.deg)) 
		sep=c1.separation(c2)
                #below adds the velocity to the array containing the smoothed circle data
		if sep.arcsecond<275:
			velocities.append(AGo_v[j])
        #below adds the mean of all of the velocities to the points within the circle to the smoothed array
        #if star[i] has fewer than 15 neighbors, it cannot be used as a center but can be used as a neighbor
	if len(velocities)>15:
		AGo_smoothedv.append(np.mean(velocities))
		AGo_xi_goodcenter.append(AGo_xi[i])
		AGo_eta_goodcenter.append(AGo_eta[i])
		AGo_dispersion.append(np.std(velocities))

#write dispersions to a file
file=open('/Users/amandaquirk/Documents/AsymmetricDrift/Data/AGo_dispersion.txt', 'w')
file.write('#sigma_LOS (km/s)\n')
for i in range(len(AGo_xi_goodcenter)):
	file.write('{}\n'.format(AGo_dispersion[i]))
file.close()

#below writes data to a file
#file=open('/Users/amandaquirk/Documents/AsymmetricDrift/Data/AGo_smoothed_goodcenters.txt', 'w')
#file.write('#eta (deg), xi (deg), average v(km/s)\n')
#for i in range(len(AGo_xi_goodcenter)):
#	file.write('{} {} {} {} {}\n'.format(AGo_eta_goodcenter[i],AGo_xi_goodcenter[i], AGo_smoothedv[i]))
#file.close()

#below creates a 3D plot showing the physical location of the stars and their LOS velocities
#plt.scatter(AGo_xi_goodcenter, AGo_eta_goodcenter, s=1, c=AGo_smoothedv, cmap='rainbow',vmin=-300,vmax=0)
#plt.xlim(-.5,1)
#plt.gca().invert_xaxis()
#plt.ylim(-1,1.2)
#plt.xlabel('xi (deg)')
#plt.ylabel('eta (deg)')
#clb=plt.colorbar()
#clb.set_label('v (km/s)')
#plt.savefig('/Users/amandaquirk/Documents/AsymmetricDrift/Data/Plots/AGo_avg_map_fixed.png')
#plt.close()

#AGB young
AGy_smoothedv=[] #will contain the averaged velocity of the smoothed circle centered on each point, if a good center
AGy_xi_goodcenter=[] #will contain x coordinate of star that is a good center
AGy_eta_goodcenter=[] #will contain y coordinate of star that is a good center
AGy_dispersion=[]
for i in range(len(AGy_xi)):
        #picks out the first point
	c1=SkyCoord(AGy_xi[i], AGy_eta[i], unit=(u.deg,u.deg)) 
        #will contain all of the velocities that are within the circle- the average of this will be added to the smoothed v array
	velocities=[]
	for j in range(len(AGy_xi)):
                #below will go through all points to determine if the points are within the smoothing circle
		c2=SkyCoord(AGy_xi[j], AGy_eta[j], unit=(u.deg,u.deg)) 
		sep=c1.separation(c2)
                #below adds the velocity to the array containing the smoothed circle data
		if sep.arcsecond<275:
			velocities.append(AGy_v[j])
        #below adds the mean of all of the velocities to the points within the circle to the smoothed array
        #if star[i] has fewer than 15 neighbors, it cannot be used as a center but can be used as a neighbor
	if len(velocities)>15:
		AGy_smoothedv.append(np.mean(velocities))
		AGy_xi_goodcenter.append(AGy_xi[i])
		AGy_eta_goodcenter.append(AGy_eta[i])
		AGy_dispersion.append(np.std(velocities))

#write dispersions to a file
file=open('/Users/amandaquirk/Documents/AsymmetricDrift/Data/AGy_dispersion.txt', 'w')
file.write('#sigma_LOS (km/s)\n')
for i in range(len(AGy_xi_goodcenter)):
	file.write('{}\n'.format(AGy_dispersion[i]))
file.close()

#below writes data to a file
#file=open('/Users/amandaquirk/Documents/AsymmetricDrift/Data/AGy_smoothed_goodcenters.txt', 'w')
#file.write('#eta (deg), xi (deg), average v(km/s)\n')
#for i in range(len(AGy_xi_goodcenter)):
#	file.write('{} {} {} {} {}\n'.format(AGy_eta_goodcenter[i],AGy_xi_goodcenter[i], AGy_smoothedv[i]))
#file.close()

#below creates a 3D plot showing the physical location of the stars and their LOS velocities
#plt.scatter(AGy_xi_goodcenter, AGy_eta_goodcenter, s=1, c=AGy_smoothedv, cmap='rainbow',vmin=-300,vmax=0)
#plt.xlim(-.5,1)
#plt.gca().invert_xaxis()
#plt.ylim(-1,1.2)
#plt.xlabel('xi (deg)')
#plt.ylabel('eta (deg)')
#clb=plt.colorbar()
#clb.set_label('v (km/s)')
#plt.savefig('/Users/amandaquirk/Documents/AsymmetricDrift/Data/Plots/AGy_avg_map_fixed.png')
#plt.close()
