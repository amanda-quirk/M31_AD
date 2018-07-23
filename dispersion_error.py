import matplotlib.pyplot as plt 
import numpy as np 
from astropy import units as u
from astropy.coordinates import SkyCoord

MS_xi, MS_eta, MS_var=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/MS_smoothed_chemin.txt', usecols=(0,1,4,), unpack=True)
AGy_xi, AGy_eta, AGy_var=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/AGy_smoothed_chemin.txt', usecols=(0,1,4,), unpack=True)
AGo_xi, AGo_eta, AGo_var=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/AGo_smoothed_chemin.txt', usecols=(0,1,4,), unpack=True)
RG_xi, RG_eta, RG_var=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/RG_smoothed_chemin.txt', usecols=(0,1,4,), unpack=True)

MS_disp_err=[]
for i in range(len(MS_xi)):
	#picks out the first point
	c1=SkyCoord(MS_xi[i]/13.67, MS_eta[i]/13.67, unit=(u.deg,u.deg)) 
	#will contain all of the velocities that are within the circle- the average of this will be added to the smoothed v array
	dispersions=[]
	for j in range(len(MS_xi)):
		#below will go through all points to determine if the points are within the smoothing circle
		c2=SkyCoord(MS_xi[j]/13.67, MS_eta[j]/13.67, unit=(u.deg,u.deg))
		sep=c1.separation(c2)
		#below adds the dispersion to the array containing the smoothed circle data
		if sep.arcsecond<200:
			dispersions.append(MS_var[j])
	#below computes the standard error
	std=np.std(dispersions)
	N=len(dispersions)
	MS_disp_err.append(std / np.sqrt(N))

AGy_disp_err=[]
for i in range(len(AGy_xi)):
	#picks out the first point
	c1=SkyCoord(AGy_xi[i]/13.67, AGy_eta[i]/13.67, unit=(u.deg,u.deg)) 
	#will contain all of the velocities that are within the circle- the average of this will be added to the smoothed v array
	dispersions=[]
	for j in range(len(AGy_xi)):
		#below will go through all points to determine if the points are within the smoothing circle
		c2=SkyCoord(AGy_xi[j]/13.67, AGy_eta[j]/13.67, unit=(u.deg,u.deg))
		sep=c1.separation(c2)
		#below adds the dispersion to the array containing the smoothed circle data
		if sep.arcsecond<200:
			dispersions.append(AGy_var[j])
	#below computes the standard error
	std=np.std(dispersions)
	N=len(dispersions)
	AGy_disp_err.append(std / np.sqrt(N))

AGo_disp_err=[]
for i in range(len(AGo_xi)):
	#picks out the first point
	c1=SkyCoord(AGo_xi[i]/13.67, AGo_eta[i]/13.67, unit=(u.deg,u.deg)) 
	#will contain all of the velocities that are within the circle- the average of this will be added to the smoothed v array
	dispersions=[]
	for j in range(len(AGo_xi)):
		#below will go through all points to determine if the points are within the smoothing circle
		c2=SkyCoord(AGo_xi[j]/13.67, AGo_eta[j]/13.67, unit=(u.deg,u.deg))
		sep=c1.separation(c2)
		#below adds the dispersion to the array containing the smoothed circle data
		if sep.arcsecond<200:
			dispersions.append(AGo_var[j])
	#below computes the standard error
	std=np.std(dispersions)
	N=len(dispersions)
	AGo_disp_err.append(std / np.sqrt(N))

RG_disp_err=[]
for i in range(len(RG_xi)):
	#picks out the first point
	c1=SkyCoord(RG_xi[i]/13.67, RG_eta[i]/13.67, unit=(u.deg,u.deg)) 
	#will contain all of the velocities that are within the circle- the average of this will be added to the smoothed v array
	dispersions=[]
	for j in range(len(RG_xi)):
		#below will go through all points to determine if the points are within the smoothing circle
		c2=SkyCoord(RG_xi[j]/13.67, RG_eta[j]/13.67, unit=(u.deg,u.deg))
		sep=c1.separation(c2)
		#below adds the dispersion to the array containing the smoothed circle data
		if sep.arcsecond<200:
			dispersions.append(RG_var[j])
	#below computes the standard error
	std=np.std(dispersions)
	N=len(dispersions)
	RG_disp_err.append(std / np.sqrt(N))

file=open('/Users/amandaquirk/Desktop/MS_dispersion_err.txt', 'w')
file.write('#dispersion error (km/s)\n')
for i in range(len(MS_xi)):
	file.write('{}\n'.format(MS_disp_err[i]))
file.close()

file=open('/Users/amandaquirk/Desktop/AGy_dispersion_err.txt', 'w')
file.write('#dispersion error (km/s)\n')
for i in range(len(AGy_xi)):
	file.write('{}\n'.format(AGy_disp_err[i]))
file.close()

file=open('/Users/amandaquirk/Desktop/AGo_dispersion_err.txt', 'w')
file.write('#dispersion error (km/s)\n')
for i in range(len(AGo_xi)):
	file.write('{}\n'.format(MS_disp_err[i]))
file.close()

file=open('/Users/amandaquirk/Desktop/RG_dispersion_err.txt', 'w')
file.write('#dispersion error (km/s)\n')
for i in range(len(RG_xi)):
	file.write('{}\n'.format(RG_disp_err[i]))
file.close()