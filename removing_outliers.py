import matplotlib.pyplot as plt 
import numpy as np 
from astropy import units as u
from astropy.coordinates import SkyCoord

#read in data with header #xi (deg), eta (deg), v(km/s), v err, n, HI main, HI close, ID, orginal index
MS_xi, MS_eta, MS_v, MS_err, MS_n, MS_HImain, MS_HIclose, MS_ind=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/MS_individual_chemin.txt', usecols=(0,1,2,3,4,5,6,8,), unpack=True)
AGy_xi, AGy_eta, AGy_v, AGy_err, AGy_n, AGy_HImain, AGy_HIclose, AGy_ind=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/AGy_individual_chemin.txt', usecols=(0,1,2,3,4,5,6,8,), unpack=True)
AGo_xi, AGo_eta, AGo_v, AGo_err, AGo_n, AGo_HImain, AGo_HIclose, AGo_ind=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/AGo_individual_chemin.txt', usecols=(0,1,2,3,4,5,6,8,), unpack=True)
RG_xi, RG_eta, RG_v, RG_err, RG_n, RG_HImain, RG_HIclose, RG_ind=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/RG_individual_chemin.txt', usecols=(0,1,2,3,4,5,6,8,), unpack=True)

#function to calculate the weights
def weights(err):
	return 1 / (err**2)

MS_weight=weights(MS_err)
AGy_weight=weights(AGy_err)
AGo_weight=weights(AGo_err)
RG_weight=weights(RG_err)

#function to normalize the weights
def normed_weight(w):
	sum_weights=sum(w)
	return w / sum_weights 

#function does the weighted meean
def weighted_mean(data,norm_w):
	return sum([a*b for a,b in zip(data, norm_w)])

#function does the weighted RMSE
def weighted_rmse(norm_w, data, mean):
	differences=[x - mean for x in data]
	diff_sq=[d**2 for d in differences]
	products=[a*b for a,b in zip(diff_sq, norm_w)]
	return np.sqrt(sum(products))

#below smoothes and removes outliers
#MS
MSs_vs=[]
MSs_vHImain=[]
MSs_vHIclose=[]
MSs_n=[]
MSs_xi=[]
MSs_eta=[]
MSs_err=[]
MSs_dispersion=[]
MSs_index=[]
for i in range(len(MS_xi)):
	#picks out the first point
	c1=SkyCoord(MS_xi[i], MS_eta[i], unit=(u.deg,u.deg)) 
	#will contain all of the velocities that are within the circle- the average of this will be added to the smoothed v array
	velocities=[]
	weights=[]
	for j in range(len(MS_xi)):
		#below will go through all points to determine if the points are within the smoothing circle
		c2=SkyCoord(MS_xi[j], MS_eta[j], unit=(u.deg,u.deg))
		sep=c1.separation(c2)
		#below adds the velocity to the array containing the smoothed circle data
		if sep.arcsecond<200:
			velocities.append(MS_v[j])
			weights.append(MS_weight[j])
	#need to remove outliers
	vel_no_outliers=[]
	weights_no_outliers=[]
	median=np.median(velocities)
	std=np.std(velocities)
	for k in range(len(velocities)):
		if abs(median-velocities[k]) / std < 3:
			vel_no_outliers.append(velocities[k])
			weights_no_outliers.append(weights[k])
	#below adds the mean of all of the velocities to the points within the circle to the smoothed array
	#if star[i] has fewer than 15 neighbors, it cannot be used as a center but can be used as a neighbor
	if len(vel_no_outliers)>15:
		normed_weights=normed_weight(weights_no_outliers)
		avg=weighted_mean(vel_no_outliers, normed_weights)
		MSs_vs.append(avg)
		MSs_xi.append(MS_xi[i]*13.67) #kpc
		MSs_eta.append(MS_eta[i]*13.67) #kpc
		MSs_vHImain.append(MS_HImain[i])
		MSs_vHIclose.append(MS_HIclose[i])
		MSs_n.append(MS_n[i])
		MSs_err.append(MS_err[i])
		MSs_dispersion.append(weighted_rmse(normed_weights,vel_no_outliers,avg))
		MSs_index.append(MS_ind[i])

#AGy
AGys_vs=[]
AGys_vHImain=[]
AGys_vHIclose=[]
AGys_n=[]
AGys_xi=[]
AGys_eta=[]
AGys_err=[]
AGys_dispersion=[]
AGys_index=[]
for i in range(len(AGy_xi)):
	#picks out the first point
	c1=SkyCoord(AGy_xi[i], AGy_eta[i], unit=(u.deg,u.deg)) 
	#will contain all of the velocities that are within the circle- the average of this will be added to the smoothed v array
	velocities=[]
	weights=[]
	for j in range(len(AGy_xi)):
		#below will go through all points to determine if the points are within the smoothing circle
		c2=SkyCoord(AGy_xi[j], AGy_eta[j], unit=(u.deg,u.deg))
		sep=c1.separation(c2)
		#below adds the velocity to the array containing the smoothed circle data
		if sep.arcsecond<200:
			velocities.append(AGy_v[j])
			weights.append(AGy_weight[j])
	#need to remove outliers
	vel_no_outliers=[]
	weights_no_outliers=[]
	median=np.median(velocities)
	std=np.std(velocities)
	for k in range(len(velocities)):
		if abs(median-velocities[k]) / std < 3:
			vel_no_outliers.append(velocities[k])
			weights_no_outliers.append(weights[k])
	#below adds the mean of all of the velocities to the points within the circle to the smoothed array
	#if star[i] has fewer than 15 neighbors, it cannot be used as a center but can be used as a neighbor
	if len(vel_no_outliers)>15:
		normed_weights=normed_weight(weights_no_outliers)
		avg=weighted_mean(vel_no_outliers, normed_weights)
		AGys_vs.append(avg)
		AGys_xi.append(AGy_xi[i]*13.67) #kpc
		AGys_eta.append(AGy_eta[i]*13.67) #kpc
		AGys_vHImain.append(AGy_HImain[i])
		AGys_vHIclose.append(AGy_HIclose[i])
		AGys_n.append(AGy_n[i])
		AGys_err.append(AGy_err[i])
		AGys_dispersion.append(weighted_rmse(normed_weights,vel_no_outliers,avg))
		AGys_index.append(AGy_ind[i])

#AGo
AGos_vs=[]
AGos_vHImain=[]
AGos_vHIclose=[]
AGos_n=[]
AGos_xi=[]
AGos_eta=[]
AGos_err=[]
AGos_dispersion=[]
AGos_index=[]
for i in range(len(AGo_xi)):
	#picks out the first point
	c1=SkyCoord(AGo_xi[i], AGo_eta[i], unit=(u.deg,u.deg)) 
	#will contain all of the velocities that are within the circle- the average of this will be added to the smoothed v array
	velocities=[]
	weights=[]
	for j in range(len(AGo_xi)):
		#below will go through all points to determine if the points are within the smoothing circle
		c2=SkyCoord(AGo_xi[j], AGo_eta[j], unit=(u.deg,u.deg))
		sep=c1.separation(c2)
		#below adds the velocity to the array containing the smoothed circle data
		if sep.arcsecond<200:
			velocities.append(AGo_v[j])
			weights.append(AGo_weight[j])
	#need to remove outliers
	vel_no_outliers=[]
	weights_no_outliers=[]
	median=np.median(velocities)
	std=np.std(velocities)
	for k in range(len(velocities)):
		if abs(median-velocities[k]) / std < 3:
			vel_no_outliers.append(velocities[k])
			weights_no_outliers.append(weights[k])
	#below adds the mean of all of the velocities to the points within the circle to the smoothed array
	#if star[i] has fewer than 15 neighbors, it cannot be used as a center but can be used as a neighbor
	if len(vel_no_outliers)>15:
		normed_weights=normed_weight(weights_no_outliers)
		avg=weighted_mean(vel_no_outliers, normed_weights)
		AGos_vs.append(avg)
		AGos_xi.append(AGo_xi[i]*13.67) #kpc
		AGos_eta.append(AGo_eta[i]*13.67) #kpc
		AGos_vHImain.append(AGo_HImain[i])
		AGos_vHIclose.append(AGo_HIclose[i])
		AGos_n.append(AGo_n[i])
		AGos_err.append(AGo_err[i])
		AGos_dispersion.append(weighted_rmse(normed_weights,vel_no_outliers,avg))
		AGos_index.append(AGo_ind[i])

#RG
RGs_vs=[]
RGs_vHImain=[]
RGs_vHIclose=[]
RGs_n=[]
RGs_xi=[]
RGs_eta=[]
RGs_err=[]
RGs_dispersion=[]
RGs_index=[]
for i in range(len(RG_xi)):
	#picks out the first point
	c1=SkyCoord(RG_xi[i], RG_eta[i], unit=(u.deg,u.deg)) 
	#will contain all of the velocities that are within the circle- the average of this will be added to the smoothed v array
	velocities=[]
	weights=[]
	for j in range(len(RG_xi)):
		#below will go through all points to determine if the points are within the smoothing circle
		c2=SkyCoord(RG_xi[j], RG_eta[j], unit=(u.deg,u.deg))
		sep=c1.separation(c2)
		#below adds the velocity to the array containing the smoothed circle data
		if sep.arcsecond<200:
			velocities.append(RG_v[j])
			weights.append(RG_weight[j])
	#need to remove outliers
	vel_no_outliers=[]
	weights_no_outliers=[]
	median=np.median(velocities)
	std=np.std(velocities)
	for k in range(len(velocities)):
		if abs(median-velocities[k]) / std < 3:
			vel_no_outliers.append(velocities[k])
			weights_no_outliers.append(weights[k])
	#below adds the mean of all of the velocities to the points within the circle to the smoothed array
	#if star[i] has fewer than 15 neighbors, it cannot be used as a center but can be used as a neighbor
	if len(vel_no_outliers)>15:
		normed_weights=normed_weight(weights_no_outliers)
		avg=weighted_mean(vel_no_outliers, normed_weights)
		RGs_vs.append(avg)
		RGs_xi.append(RG_xi[i]*13.67) #kpc
		RGs_eta.append(RG_eta[i]*13.67) #kpc
		RGs_vHImain.append(RG_HImain[i])
		RGs_vHIclose.append(RG_HIclose[i])
		RGs_n.append(RG_n[i])
		RGs_err.append(RG_err[i])
		RGs_dispersion.append(weighted_rmse(normed_weights,vel_no_outliers,avg))
		RGs_index.append(RG_ind[i])

#plotting
plt.figure(figsize=(4, 6))
plt.scatter(MSs_xi, MSs_eta, s=1, c=MSs_vs, cmap='rainbow',vmin=-300,vmax=0)
plt.xlim(2.5,12,5)
plt.gca().invert_xaxis()
plt.ylim(-2.5,15)
plt.xlabel('xi (kpc)')
plt.ylabel('eta (kpc)')
clb=plt.colorbar()
clb.set_label('v (km/s)')
plt.savefig('/Users/amandaquirk/Desktop/MS_smoothed_chemin_no_outliers.png')
plt.close()

plt.figure(figsize=(4, 6))
plt.scatter(MSs_xi, MSs_eta, s=1, c=MSs_dispersion, cmap='copper',vmin=10,vmax=150)
plt.xlim(2.5,12,5)
plt.gca().invert_xaxis()
plt.ylim(-2.5,15)
plt.xlabel('xi (kpc)')
plt.ylabel('eta (kpc)')
clb=plt.colorbar()
clb.set_label('v (km/s)')
plt.savefig('/Users/amandaquirk/Desktop/MS_dispersion_chemin_no_outliers.png')
plt.close()

#below writes data to a file
file=open('/Users/amandaquirk/Desktop/MS_smoothed_chemin_no_outliers.txt', 'w')
file.write('#xi (kpc), eta (kpc), average v(km/s), adjusted v err, var, n, HI main, HI close, ID, orginal index\n')
for i in range(len(MSs_xi)):
	file.write('{} {} {} {} {} {} {} {} {}\n'.format(MSs_xi[i],MSs_eta[i], MSs_vs[i], MSs_err[i], MSs_dispersion[i], MSs_n[i], MSs_vHImain[i], MSs_vHIclose[i], MSs_index[i]))
file.close()

#plotting
plt.figure(figsize=(4, 6))
plt.scatter(AGys_xi, AGys_eta, s=1, c=AGys_vs, cmap='rainbow',vmin=-300,vmax=0)
plt.xlim(2.5,12,5)
plt.gca().invert_xaxis()
plt.ylim(-2.5,15)
plt.xlabel('xi (kpc)')
plt.ylabel('eta (kpc)')
clb=plt.colorbar()
clb.set_label('v (km/s)')
plt.savefig('/Users/amandaquirk/Desktop/AGy_smoothed_chemin_no_outliers.png')
plt.close()

plt.figure(figsize=(4, 6))
plt.scatter(AGys_xi, AGys_eta, s=1, c=AGys_dispersion, cmap='copper',vmin=10,vmax=150)
plt.xlim(2.5,12,5)
plt.gca().invert_xaxis()
plt.ylim(-2.5,15)
plt.xlabel('xi (kpc)')
plt.ylabel('eta (kpc)')
clb=plt.colorbar()
clb.set_label('v (km/s)')
plt.savefig('/Users/amandaquirk/Desktop/AGy_dispersion_chemin_no_outliers.png')
plt.close()

#below writes data to a file
file=open('/Users/amandaquirk/Desktop/AGy_smoothed_chemin_no_outliers.txt', 'w')
file.write('#xi (kpc), eta (kpc), average v(km/s), adjusted v err, var, n, HI main, HI close, ID, orginal index\n')
for i in range(len(AGys_xi)):
	file.write('{} {} {} {} {} {} {} {} {}\n'.format(AGys_xi[i],AGys_eta[i], AGys_vs[i], AGys_err[i], AGys_dispersion[i], AGys_n[i], AGys_vHImain[i], AGys_vHIclose[i], AGys_index[i]))
file.close()

#plotting
plt.figure(figsize=(4, 6))
plt.scatter(AGos_xi, AGos_eta, s=1, c=AGos_vs, cmap='rainbow',vmin=-300,vmax=0)
plt.xlim(2.5,12,5)
plt.gca().invert_xaxis()
plt.ylim(-2.5,15)
plt.xlabel('xi (kpc)')
plt.ylabel('eta (kpc)')
clb=plt.colorbar()
clb.set_label('v (km/s)')
plt.savefig('/Users/amandaquirk/Desktop/AGo_smoothed_chemin_no_outliers.png')
plt.close()

plt.figure(figsize=(4, 6))
plt.scatter(AGos_xi, AGos_eta, s=1, c=AGos_dispersion, cmap='copper',vmin=10,vmax=150)
plt.xlim(2.5,12,5)
plt.gca().invert_xaxis()
plt.ylim(-2.5,15)
plt.xlabel('xi (kpc)')
plt.ylabel('eta (kpc)')
clb=plt.colorbar()
clb.set_label('v (km/s)')
plt.savefig('/Users/amandaquirk/Desktop/AGo_dispersion_chemin_no_outliers.png')
plt.close()

#below writes data to a file
file=open('/Users/amandaquirk/Desktop/AGo_smoothed_chemin_no_outliers.txt', 'w')
file.write('#xi (kpc), eta (kpc), average v(km/s), adjusted v err, var, n, HI main, HI close, ID, orginal index\n')
for i in range(len(AGos_xi)):
	file.write('{} {} {} {} {} {} {} {} {}\n'.format(AGos_xi[i],AGos_eta[i], AGos_vs[i], AGos_err[i], AGos_dispersion[i],AGos_n[i], AGos_vHImain[i], AGos_vHIclose[i], AGos_index[i]))
file.close()

#plotting
plt.figure(figsize=(4, 6))
plt.scatter(RGs_xi, RGs_eta, s=1, c=RGs_vs, cmap='rainbow',vmin=-300,vmax=0)
plt.xlim(2.5,12,5)
plt.gca().invert_xaxis()
plt.ylim(-2.5,15)
plt.xlabel('xi (kpc)')
plt.ylabel('eta (kpc)')
clb=plt.colorbar()
clb.set_label('v (km/s)')
plt.savefig('/Users/amandaquirk/Desktop/RG_smoothed_chemin_no_outliers.png')
plt.close()

plt.figure(figsize=(4, 6))
plt.scatter(RGs_xi, RGs_eta, s=1, c=RGs_dispersion, cmap='copper',vmin=10,vmax=150)
plt.xlim(2.5,12,5)
plt.gca().invert_xaxis()
plt.ylim(-2.5,15)
plt.xlabel('xi (kpc)')
plt.ylabel('eta (kpc)')
clb=plt.colorbar()
clb.set_label('v (km/s)')
plt.savefig('/Users/amandaquirk/Desktop/RG_dispersion_chemin_no_outliers.png')
plt.close()

#below writes data to a file
file=open('/Users/amandaquirk/Desktop/RG_smoothed_chemin_no_outliers.txt', 'w')
file.write('#xi (kpc), eta (kpc), average v(km/s), adjusted v err, var, n, HI main, HI close, orginal index\n')
for i in range(len(RGs_xi)):
	file.write('{} {} {} {} {} {} {} {} {}\n'.format(RGs_xi[i],RGs_eta[i], RGs_vs[i], RGs_err[i], RGs_dispersion[i], RGs_n[i], RGs_vHImain[i], RGs_vHIclose[i], RGs_index[i]))
file.close()