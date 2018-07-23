from astropy import units as u
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt 
import numpy as np 
from astropy.io import fits 
from astropy.table import Table 

#read in the data and convert to a fits table
hdulist=fits.open('/Users/amandaquirk/Documents/AsymmetricDrift/Data/submasterSPLASH_quirk_master.fits', memmap=True)

#puts data into an accessible table
data=Table(hdulist[1].data)

#extract the data from the table
RA=data['RA'] #coordinate
Dec=data['Dec'] #coordinate
v_s=data['v_s'] #averaged LOS velocity of star
v_HImain=data['HI_main'] #velocity of HI main component 
v_HIclose=data['HI_close'] #velocity of HI component closest to the star
n=data['N'] #number of peaks in HI velocity spectrum: 1-5
age=data['Age'] #age tag: 1-4, 1 for MS and 4 for RG
error=data['verr'] #velocity error
index=data['Index in submasterSPLASH']
ID=data['ID']

#convert things to floats
v_s=v_s.astype(np.float)
v_HImain=v_HImain.astype(np.float)
v_HIclose=v_HIclose.astype(np.float)
n=n.astype(np.float)
age=age.astype(np.float)
error=error.astype(np.float)
index=index.astype(np.float)

#function to adjust the errors the way claire does
def adjust_err(err):
	return np.sqrt((3.2 * err)**2 + 1.2)

adjusted_err=adjust_err(error)

#function to calculate the weights
def weights(err):
	return 1 / (err**2)

weight=weights(adjusted_err)

#function to normalize the weights
def normed_weight(w):
	sum_weights=sum(w)
	return w / sum_weights 

#sorting the stars into the age groups
MS_xi=[] #will contain coordinates centered on M31's center
MS_eta=[] #will contain coordinates centered on M31's center
MS_vs=[]
MS_vHImain=[]
MS_vHIclose=[]
MS_n=[]
MS_err=[]
MS_index=[]
MS_ID=[]
MS_weight=[]
AGy_xi=[] #will contain coordinates centered on M31's center
AGy_eta=[] #will contain coordinates centered on M31's center
AGy_vs=[]
AGy_vHImain=[]
AGy_vHIclose=[]
AGy_n=[]
AGy_err=[]
AGy_index=[]
AGy_ID=[]
AGy_weight=[]
AGo_xi=[] #will contain coordinates centered on M31's center
AGo_eta=[] #will contain coordinates centered on M31's center
AGo_vs=[]
AGo_vHImain=[]
AGo_vHIclose=[]
AGo_n=[]
AGo_err=[]
AGo_index=[]
AGo_ID=[]
AGo_weight=[]
RG_xi=[] #will contain coordinates centered on M31's center
RG_eta=[] #will contain coordinates centered on M31's center
RG_vs=[]
RG_vHImain=[]
RG_vHIclose=[]
RG_n=[]
RG_err=[]
RG_index=[]
RG_ID=[]
RG_weight=[]

m31 = SkyCoord(10.6847083*u.deg, 41.26875*u.deg, frame='icrs') #will be used to convert RA and Dec to xi and eta
for i in range(len(age)):
	if age[i]==1:
		MS_vs.append(v_s[i])
		MS_vHImain.append(v_HImain[i])
		MS_vHIclose.append(v_HIclose[i])
		MS_n.append(n[i])
		c=SkyCoord(RA[i], Dec[i], frame='icrs', unit=(u.hourangle,u.deg))
		c_inm31=c.transform_to(m31.skyoffset_frame())
		xi, eta=c_inm31.lon, c_inm31.lat
		MS_xi.append(xi.degree)
		MS_eta.append(eta.degree)
		MS_err.append(adjusted_err[i])
		MS_index.append(index[i])
		MS_ID.append(ID[i])
		MS_weight.append(weight[i]) 

	if age[i]==2:
		AGy_vs.append(v_s[i])
		AGy_vHImain.append(v_HImain[i])
		AGy_vHIclose.append(v_HIclose[i])
		AGy_n.append(n[i])
		c=SkyCoord(RA[i], Dec[i], frame='icrs', unit=(u.hourangle,u.deg))
		c_inm31=c.transform_to(m31.skyoffset_frame())
		xi, eta=c_inm31.lon, c_inm31.lat
		AGy_xi.append(xi.degree)
		AGy_eta.append(eta.degree)
		AGy_err.append(adjusted_err[i])
		AGy_index.append(index[i])
		AGy_ID.append(ID[i])
		AGy_weight.append(weight[i])

	if age[i]==3:
		AGo_vs.append(v_s[i])
		AGo_vHImain.append(v_HImain[i])
		AGo_vHIclose.append(v_HIclose[i])
		AGo_n.append(n[i])
		c=SkyCoord(RA[i], Dec[i], frame='icrs', unit=(u.hourangle,u.deg))
		c_inm31=c.transform_to(m31.skyoffset_frame())
		xi, eta=c_inm31.lon, c_inm31.lat
		AGo_xi.append(xi.degree)
		AGo_eta.append(eta.degree)
		AGo_err.append(adjusted_err[i])
		AGo_index.append(index[i])
		AGo_ID.append(ID[i])
		AGo_weight.append(weight[i])

	if age[i]==4:
		RG_vs.append(v_s[i])
		RG_vHImain.append(v_HImain[i])
		RG_vHIclose.append(v_HIclose[i])
		RG_n.append(n[i])
		c=SkyCoord(RA[i], Dec[i], frame='icrs', unit=(u.hourangle,u.deg))
		c_inm31=c.transform_to(m31.skyoffset_frame())
		xi, eta=c_inm31.lon, c_inm31.lat
		RG_xi.append(xi.degree)
		RG_eta.append(eta.degree)
		RG_err.append(adjusted_err[i])
		RG_index.append(index[i])
		RG_ID.append(ID[i])
		RG_weight.append(weight[i])

print(len(MS_xi) ,len(MS_eta), len(MS_vs), len(MS_err), len(MS_n), len(MS_vHImain), len(MS_vHIclose), len(MS_ID), len(MS_index))

#below writes data to a file
file=open('/Users/amandaquirk/Desktop/MS_individual_chemin.txt', 'w')
file.write('#xi (kpc), eta (kpc), v(km/s), v err, n, HI main, HI close, ID, orginal index\n')
for i in range(len(MS_xi)):
	file.write('{} {} {} {} {} {} {} {} {}\n'.format(MS_xi[i],MS_eta[i], MS_vs[i], MS_err[i], MS_n[i], MS_vHImain[i], MS_vHIclose[i], MS_ID[i], MS_index[i]))
file.close()

file=open('/Users/amandaquirk/Desktop/AGy_individual_chemin.txt', 'w')
file.write('#xi (kpc), eta (kpc), v(km/s), v err, n, HI main, HI close, ID, orginal index\n')
for i in range(len(AGy_xi)):
	file.write('{} {} {} {} {} {} {} {} {}\n'.format(AGy_xi[i],AGy_eta[i], AGy_vs[i], AGy_err[i], AGy_n[i], AGy_vHImain[i], AGy_vHIclose[i], AGy_ID[i], AGy_index[i]))
file.close()

file=open('/Users/amandaquirk/Desktop/AGo_individual_chemin.txt', 'w')
file.write('#xi (kpc), eta (kpc), v(km/s), v err, n, HI main, HI close, ID, orginal index\n')
for i in range(len(AGo_xi)):
	file.write('{} {} {} {} {} {} {} {} {}\n'.format(AGo_xi[i],AGo_eta[i], AGo_vs[i], AGo_err[i], AGo_n[i], AGo_vHImain[i], AGo_vHIclose[i], AGo_ID[i], AGo_index[i]))
file.close()

file=open('/Users/amandaquirk/Desktop/RG_individual_chemin.txt', 'w')
file.write('#xi (kpc), eta (kpc), v(km/s), v err, n, HI main, HI close, ID, orginal index\n')
for i in range(len(RG_xi)):
	file.write('{} {} {} {} {} {} {} {} {}\n'.format(RG_xi[i],RG_eta[i], RG_vs[i], RG_err[i], RG_n[i], RG_vHImain[i], RG_vHIclose[i], RG_ID[i], RG_index[i]))
file.close()

#individual velocity plots
# plt.figure(figsize=(4, 6))
# plt.scatter(MS_xi, MS_eta, s=1, c=MS_vs, cmap='rainbow',vmin=-300,vmax=0)
# plt.gca().invert_xaxis()
# plt.xlabel('xi (deg)')
# plt.ylabel('eta (deg)')
# clb=plt.colorbar()
# clb.set_label('v (km/s)')
# plt.savefig('/Users/amandaquirk/Desktop/MS_individual_chemin.png')
# plt.close()

# plt.figure(figsize=(4, 6))
# plt.scatter(AGy_xi, AGy_eta, s=1, c=AGy_vs, cmap='rainbow',vmin=-300,vmax=0)
# plt.gca().invert_xaxis()
# plt.xlabel('xi (deg)')
# plt.ylabel('eta (deg)')
# clb=plt.colorbar()
# clb.set_label('v (km/s)')
# plt.savefig('/Users/amandaquirk/Desktop/AGy_individual_chemin.png')
# plt.close()

# plt.figure(figsize=(4, 6))
# plt.scatter(AGo_xi, AGo_eta, s=1, c=AGo_vs, cmap='rainbow',vmin=-300,vmax=0)
# plt.gca().invert_xaxis()
# plt.xlabel('xi (deg)')
# plt.ylabel('eta (deg)')
# clb=plt.colorbar()
# clb.set_label('v (km/s)')
# plt.savefig('/Users/amandaquirk/Desktop/AGo_individual_chemin.png')
# plt.close()

# plt.figure(figsize=(4, 6))
# plt.scatter(RG_xi, RG_eta, s=1, c=RG_vs, cmap='rainbow',vmin=-300,vmax=0)
# plt.gca().invert_xaxis()
# plt.xlabel('xi (deg)')
# plt.ylabel('eta (deg)')
# clb=plt.colorbar()
# clb.set_label('v (km/s)')
# plt.savefig('/Users/amandaquirk/Desktop/RG_individual_chemin.png')
# plt.close()

# #function does the weighted meean
# def weighted_mean(data,norm_w):
# 	return sum([a*b for a,b in zip(data, norm_w)])

# #function does the weighted RMSE
# def weighted_rmse(norm_w, data, mean):
# 	differences=[x - mean for x in data]
# 	diff_sq=[d**2 for d in differences]
# 	products=[a*b for a,b in zip(diff_sq, norm_w)]
# 	return np.sqrt(sum(products))

# #MS
# MSs_vs=[]
# MSs_vHImain=[]
# MSs_vHIclose=[]
# MSs_n=[]
# MSs_xi=[]
# MSs_eta=[]
# MSs_err=[]
# MSs_dispersion=[]
# MSs_ID=[]
# MSs_index=[]
# for i in range(len(MS_xi)):
# 	#picks out the first point
# 	c1=SkyCoord(MS_xi[i], MS_eta[i], unit=(u.deg,u.deg)) 
# 	#will contain all of the velocities that are within the circle- the average of this will be added to the smoothed v array
# 	velocities=[]
# 	weights=[]
# 	for j in range(len(MS_xi)):
# 		#below will go through all points to determine if the points are within the smoothing circle
# 		c2=SkyCoord(MS_xi[j], MS_eta[j], unit=(u.deg,u.deg))
# 		sep=c1.separation(c2)
# 		#below adds the velocity to the array containing the smoothed circle data
# 		if sep.arcsecond<200:
# 			velocities.append(MS_vs[j])
# 			weights.append(MS_weight[j])
# 	#below adds the mean of all of the velocities to the points within the circle to the smoothed array
# 	#if star[i] has fewer than 15 neighbors, it cannot be used as a center but can be used as a neighbor
# 	if len(velocities)>15:
# 		normed_weights=normed_weight(weights)
# 		avg=weighted_mean(velocities, normed_weights)
# 		MSs_vs.append(avg)
# 		MSs_xi.append(MS_xi[i]*13.67) #kpc
# 		MSs_eta.append(MS_eta[i]*13.67) #kpc
# 		MSs_vHImain.append(MS_vHImain[i])
# 		MSs_vHIclose.append(MS_vHIclose[i])
# 		MSs_n.append(MS_n[i])
# 		MSs_err.append(MS_err[i])
# 		MSs_dispersion.append(weighted_rmse(normed_weights,velocities,avg))
# 		MSs_ID.append(MS_ID[i])
# 		MSs_index.append(MS_index[i])

# #plotting
# plt.figure(figsize=(4, 6))
# plt.scatter(MSs_xi, MSs_eta, s=1, c=MSs_vs, cmap='rainbow',vmin=-300,vmax=0)
# plt.xlim(2.5,12,5)
# plt.gca().invert_xaxis()
# plt.ylim(-2.5,15)
# plt.xlabel('xi (kpc)')
# plt.ylabel('eta (kpc)')
# clb=plt.colorbar()
# clb.set_label('v (km/s)')
# plt.savefig('/Users/amandaquirk/Desktop/MS_smoothed_chemin.png')
# plt.close()

# plt.figure(figsize=(4, 6))
# plt.scatter(MSs_xi, MSs_eta, s=1, c=MSs_dispersion, cmap='copper',vmin=10,vmax=150)
# plt.xlim(2.5,12,5)
# plt.gca().invert_xaxis()
# plt.ylim(-2.5,15)
# plt.xlabel('xi (kpc)')
# plt.ylabel('eta (kpc)')
# clb=plt.colorbar()
# clb.set_label('v (km/s)')
# plt.savefig('/Users/amandaquirk/Desktop/MS_dispersion_chemin.png')
# plt.close()

# #below writes data to a file
# file=open('/Users/amandaquirk/Documents/AsymmetricDrift/Data/MS_smoothed_chemin.txt', 'w')
# file.write('#xi (kpc), eta (kpc), average v(km/s), v err, var, n, HI main, HI close, ID, orginal index\n')
# for i in range(len(MSs_xi)):
# 	file.write('{} {} {} {} {} {} {} {} {} {}\n'.format(MSs_xi[i],MSs_eta[i], MSs_vs[i], MSs_err[i], MSs_dispersion[i], MSs_n[i], MSs_vHImain[i], MSs_vHIclose[i], MSs_ID[i], MSs_index[i]))
# file.close()

# #AGy
# AGys_vs=[]
# AGys_vHImain=[]
# AGys_vHIclose=[]
# AGys_n=[]
# AGys_xi=[]
# AGys_eta=[]
# AGys_err=[]
# AGys_dispersion=[]
# AGys_ID=[]
# AGys_index=[]
# for i in range(len(AGy_xi)):
# 	#picks out the first point
# 	c1=SkyCoord(AGy_xi[i], AGy_eta[i], unit=(u.deg,u.deg)) 
# 	#will contain all of the velocities that are within the circle- the average of this will be added to the smoothed v array
# 	velocities=[]
# 	weights=[]
# 	for j in range(len(AGy_xi)):
# 		#below will go through all points to determine if the points are within the smoothing circle
# 		c2=SkyCoord(AGy_xi[j], AGy_eta[j], unit=(u.deg,u.deg))
# 		sep=c1.separation(c2)
# 		#below adds the velocity to the array containing the smoothed circle data
# 		if sep.arcsecond<275:
# 			velocities.append(AGy_vs[j])
# 			weights.append(AGy_weight[j])
# 	#below adds the mean of all of the velocities to the points within the circle to the smoothed array
# 	#if star[i] has fewer than 15 neighbors, it cannot be used as a center but can be used as a neighbor
# 	if len(velocities)>15:
# 		normed_weights=normed_weight(weights)
# 		avg=weighted_mean(velocities, normed_weights)
# 		AGys_vs.append(avg)
# 		AGys_xi.append(AGy_xi[i]*13.67) #kpc
# 		AGys_eta.append(AGy_eta[i]*13.67) #kpc
# 		AGys_vHImain.append(AGy_vHImain[i])
# 		AGys_vHIclose.append(AGy_vHIclose[i])
# 		AGys_n.append(AGy_n[i])
# 		AGys_err.append(AGy_err[i])
# 		AGys_ID.append(AGy_ID[i])
# 		AGys_index.append(AGy_index[i])
# 		AGys_dispersion.append(weighted_rmse(normed_weights,velocities,avg))

# #plotting
# plt.figure(figsize=(4, 6))
# plt.scatter(AGys_xi, AGys_eta, s=1, c=AGys_vs, cmap='rainbow',vmin=-300,vmax=0)
# plt.xlim(2.5,12,5)
# plt.gca().invert_xaxis()
# plt.ylim(-2.5,15)
# plt.xlabel('xi (kpc)')
# plt.ylabel('eta (kpc)')
# clb=plt.colorbar()
# clb.set_label('v (km/s)')
# plt.savefig('/Users/amandaquirk/Desktop/AGy_smoothed_chemin.png')
# plt.close()

# plt.figure(figsize=(4, 6))
# plt.scatter(AGys_xi, AGys_eta, s=1, c=AGys_dispersion, cmap='copper',vmin=10,vmax=150)
# plt.xlim(2.5,12,5)
# plt.gca().invert_xaxis()
# plt.ylim(-2.5,15)
# plt.xlabel('xi (kpc)')
# plt.ylabel('eta (kpc)')
# clb=plt.colorbar()
# clb.set_label('v (km/s)')
# plt.savefig('/Users/amandaquirk/Desktop/AGy_dispersion_chemin.png')
# plt.close()

# #below writes data to a file
# file=open('/Users/amandaquirk/Documents/AsymmetricDrift/Data/AGy_smoothed_chemin.txt', 'w')
# file.write('#xi (kpc), eta (kpc), average v(km/s), v err, var, n, HI main, HI close, ID, orginal index\n')
# for i in range(len(AGys_xi)):
# 	file.write('{} {} {} {} {} {} {} {} {} {}\n'.format(AGys_xi[i],AGys_eta[i], AGys_vs[i], AGys_err[i], AGys_dispersion[i], AGys_n[i], AGys_vHImain[i], AGys_vHIclose[i], AGys_ID[i], AGys_index[i]))
# file.close()

# #AGo
# AGos_vs=[]
# AGos_vHImain=[]
# AGos_vHIclose=[]
# AGos_n=[]
# AGos_xi=[]
# AGos_eta=[]
# AGos_err=[]
# AGos_dispersion=[]
# AGos_ID=[]
# AGos_index=[]
# for i in range(len(AGo_xi)):
# 	#picks out the first point
# 	c1=SkyCoord(AGo_xi[i], AGo_eta[i], unit=(u.deg,u.deg)) 
# 	#will contain all of the velocities that are within the circle- the average of this will be added to the smoothed v array
# 	velocities=[]
# 	weights=[]
# 	for j in range(len(AGo_xi)):
# 		#below will go through all points to determine if the points are within the smoothing circle
# 		c2=SkyCoord(AGo_xi[j], AGo_eta[j], unit=(u.deg,u.deg))
# 		sep=c1.separation(c2)
# 		#below adds the velocity to the array containing the smoothed circle data
# 		if sep.arcsecond<275:
# 			velocities.append(AGo_vs[j])
# 			weights.append(AGo_weight[j])
# 	#below adds the mean of all of the velocities to the points within the circle to the smoothed array
# 	#if star[i] has fewer than 15 neighbors, it cannot be used as a center but can be used as a neighbor
# 	if len(velocities)>15:
# 		normed_weights=normed_weight(weights)
# 		avg=weighted_mean(velocities, normed_weights)
# 		AGos_vs.append(avg)
# 		AGos_xi.append(AGo_xi[i]*13.67) #kpc
# 		AGos_eta.append(AGo_eta[i]*13.67) #kpc
# 		AGos_vHImain.append(AGo_vHImain[i])
# 		AGos_vHIclose.append(AGo_vHIclose[i])
# 		AGos_n.append(AGo_n[i])
# 		AGos_err.append(AGo_err[i])
# 		AGos_ID.append(AGo_ID[i])
# 		AGos_index.append(AGo_index[i])
# 		AGos_dispersion.append(weighted_rmse(normed_weights,velocities,avg))

# #plotting
# plt.figure(figsize=(4, 6))
# plt.scatter(AGos_xi, AGos_eta, s=1, c=AGos_vs, cmap='rainbow',vmin=-300,vmax=0)
# plt.xlim(2.5,12,5)
# plt.gca().invert_xaxis()
# plt.ylim(-2.5,15)
# plt.xlabel('xi (kpc)')
# plt.ylabel('eta (kpc)')
# clb=plt.colorbar()
# clb.set_label('v (km/s)')
# plt.savefig('/Users/amandaquirk/Desktop/AGo_smoothed_chemin.png')
# plt.close()

# plt.figure(figsize=(4, 6))
# plt.scatter(AGos_xi, AGos_eta, s=1, c=AGos_dispersion, cmap='copper',vmin=10,vmax=150)
# plt.xlim(2.5,12,5)
# plt.gca().invert_xaxis()
# plt.ylim(-2.5,15)
# plt.xlabel('xi (kpc)')
# plt.ylabel('eta (kpc)')
# clb=plt.colorbar()
# clb.set_label('v (km/s)')
# plt.savefig('/Users/amandaquirk/Desktop/AGo_dispersion_chemin.png')
# plt.close()

# #below writes data to a file
# file=open('/Users/amandaquirk/Documents/AsymmetricDrift/Data/AGo_smoothed_chemin.txt', 'w')
# file.write('#xi (kpc), eta (kpc), average v(km/s), v err, var, n, HI main, HI close, ID, orginal index\n')
# for i in range(len(AGos_xi)):
# 	file.write('{} {} {} {} {} {} {} {} {} {}\n'.format(AGos_xi[i],AGos_eta[i], AGos_vs[i], AGos_err[i], AGos_dispersion[i],AGos_n[i], AGos_vHImain[i], AGos_vHIclose[i], AGos_ID[i], AGos_index[i]))
# file.close()

# #RG
# RGs_vs=[]
# RGs_vHImain=[]
# RGs_vHIclose=[]
# RGs_n=[]
# RGs_xi=[]
# RGs_eta=[]
# RGs_err=[]
# RGs_dispersion=[]
# RGs_ID=[]
# RGs_index=[]
# for i in range(len(RG_xi)):
# 	#picks out the first point
# 	c1=SkyCoord(RG_xi[i], RG_eta[i], unit=(u.deg,u.deg)) 
# 	#will contain all of the velocities that are within the circle- the average of this will be added to the smoothed v array
# 	velocities=[]
# 	weights=[]
# 	for j in range(len(RG_xi)):
# 		#below will go through all points to determine if the points are within the smoothing circle
# 		c2=SkyCoord(RG_xi[j], RG_eta[j], unit=(u.deg,u.deg))
# 		sep=c1.separation(c2)
# 		#below adds the velocity to the array containing the smoothed circle data
# 		if sep.arcsecond<200:
# 			velocities.append(RG_vs[j])
# 			weights.append(RG_weight[j])
# 	#below adds the mean of all of the velocities to the points within the circle to the smoothed array
# 	#if star[i] has fewer than 15 neighbors, it cannot be used as a center but can be used as a neighbor
# 	if len(velocities)>15:
# 		normed_weights=normed_weight(weights)
# 		avg=weighted_mean(velocities, normed_weights)
# 		RGs_vs.append(avg)
# 		RGs_xi.append(RG_xi[i]*13.67) #kpc
# 		RGs_eta.append(RG_eta[i]*13.67)
# 		RGs_vHImain.append(RG_vHImain[i]) #kpc
# 		RGs_vHIclose.append(RG_vHIclose[i])
# 		RGs_n.append(RG_n[i])
# 		RGs_err.append(RG_err[i])
# 		RGs_dispersion.append(weighted_rmse(normed_weights,velocities,avg))
# 		RGs_ID.append(RG_ID[i])
# 		RGs_index.append(RG_index[i])

# #plotting
# plt.figure(figsize=(4, 6))
# plt.scatter(RGs_xi, RGs_eta, s=1, c=RGs_vs, cmap='rainbow',vmin=-300,vmax=0)
# plt.xlim(2.5,12,5)
# plt.gca().invert_xaxis()
# plt.ylim(-2.5,15)
# plt.xlabel('xi (kpc)')
# plt.ylabel('eta (kpc)')
# clb=plt.colorbar()
# clb.set_label('v (km/s)')
# plt.savefig('/Users/amandaquirk/Desktop/RG_smoothed_chemin.png')
# plt.close()

# plt.figure(figsize=(4, 6))
# plt.scatter(RGs_xi, RGs_eta, s=1, c=RGs_dispersion, cmap='copper',vmin=10,vmax=150)
# plt.xlim(2.5,12,5)
# plt.gca().invert_xaxis()
# plt.ylim(-2.5,15)
# plt.xlabel('xi (kpc)')
# plt.ylabel('eta (kpc)')
# clb=plt.colorbar()
# clb.set_label('v (km/s)')
# plt.savefig('/Users/amandaquirk/Desktop/RG_dispersion_chemin.png')
# plt.close()

# #below writes data to a file
# file=open('/Users/amandaquirk/Documents/AsymmetricDrift/Data/RG_smoothed_chemin.txt', 'w')
# file.write('#xi (kpc), eta (kpc), average v(km/s), v err, var, n, HI main, HI close, ID, orginal index\n')
# for i in range(len(RGs_xi)):
# 	file.write('{} {} {} {} {} {} {} {} {} {}\n'.format(RGs_xi[i],RGs_eta[i], RGs_vs[i], RGs_err[i], RGs_dispersion[i], RGs_n[i], RGs_vHImain[i], RGs_vHIclose[i], RGs_ID[i], RGs_index[i]))
# file.close()

#adding together the HI data
# HIs_xi=MSs_xi+AGys_xi+AGos_xi+RGs_xi
# HIs_eta=MSs_eta+AGys_eta+AGos_eta+RGs_eta
# HIs_vmain=MSs_vHImain+AGys_vHImain+AGos_vHImain+RGs_vHImain
# HIs_vclose=MSs_vHIclose+AGys_vHIclose+AGos_vHIclose+RGs_vHIclose
# HIs_n=MSs_n+AGys_n+AGos_n+RGs_n

# #plotting
# plt.figure(figsize=(4, 6))
# plt.scatter(HIs_xi, HIs_eta, s=1, c=HIs_vmain, cmap='rainbow',vmin=-300,vmax=0)
# plt.xlim(-.5,1)
# plt.gca().invert_xaxis()
# plt.ylim(-0.75,1.2)
# plt.xlabel('xi (deg)')
# plt.ylabel('eta (deg)')
# clb=plt.colorbar()
# clb.set_label('v (km/s)')
# plt.savefig('/Users/amandaquirk/Desktop/HImain_smoothed_chemin.png')
# plt.close()

# plt.figure(figsize=(4, 6))
# plt.scatter(HIs_xi, HIs_eta, s=1, c=HIs_vclose, cmap='rainbow',vmin=-300,vmax=0)
# plt.xlim(-.5,1)
# plt.gca().invert_xaxis()
# plt.ylim(-0.75,1.2)
# plt.xlabel('xi (deg)')
# plt.ylabel('eta (deg)')
# clb=plt.colorbar()
# clb.set_label('v (km/s)')
# plt.savefig('/Users/amandaquirk/Desktop/HIclose_smoothed_chemin.png')
# plt.close()

# #below writes data to a file
# file=open('/Users/amandaquirk/Documents/AsymmetricDrift/Data/HI_smoothed_chemin.txt', 'w')
# file.write('#xi (deg), eta (deg), main v(km/s), close v(km/s), n\n')
# for i in range(len(HIs_xi)):
# 	file.write('{} {} {} {} {}\n'.format(HIs_xi[i],HIs_eta[i], HIs_vmain[i], HIs_vclose[i], HIs_n[i]))
# file.close()


