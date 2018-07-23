from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import match_coordinates_sky
import matplotlib.pyplot as plt 
import numpy as np 

#importing all of the data
#xi (kpc), eta (kpc), average v(km/s), v err,var, n, HI main, HI close <-- header of data file to be read in 
MS_xi, MS_eta=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/MS_smoothed_chemin.txt', usecols=(0,1,), unpack=True)
AGy_xi, AGy_eta=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/AGy_smoothed_chemin.txt', usecols=(0,1,), unpack=True)
AGo_xi, AGo_eta=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/AGo_smoothed_chemin.txt', usecols=(0,1,), unpack=True)
RG_xi, RG_eta=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/RG_smoothed_chemin.txt', usecols=(0,1,), unpack=True)

#RA, Dec, Av, Av error
dust_RA, dust_Dec, dust_Av, dust_err=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/dust_coords_medAv.txt', unpack=True)

#convert the dust coordinates to xi and eta
def m31_coord(RA, Dec, x=True):
	m31 = SkyCoord(10.6847083*u.deg, 41.26875*u.deg, frame='icrs')
	c=SkyCoord(RA, Dec, frame='icrs', unit=(u.deg,u.deg))
	c_inm31=c.transform_to(m31.skyoffset_frame())
	xi, eta=c_inm31.lon, c_inm31.lat
	if x==True:
		return xi*13.67 
	if x==False:
		return eta*13.67 

dust_xi=m31_coord(dust_RA, dust_Dec, x=True)
dust_eta=m31_coord(dust_RA, dust_Dec, x=False)

# plt.figure(figsize=(4, 6))
# plt.scatter(dust_xi, dust_eta, c=dust_Av)
# plt.scatter(RG_xi, RG_eta, s=2, c='r')
# plt.gca().invert_xaxis()
# plt.xlabel('xi (kpc)')
# plt.ylabel('eta (kpc)')
# plt.show()

#turning xi and eta into a coordinate pair
def coords(xi, eta):
	return SkyCoord(xi/13.67, eta/13.67, unit=(u.deg, u.deg)) #dividing by 13.67 to convert to degrees

dust_coords=coords(dust_xi, dust_eta)
MS_coords=coords(MS_xi, MS_eta)
AGy_coords=coords(AGy_xi, AGy_eta)
AGo_coords=coords(AGo_xi, AGo_eta)
RG_coords=coords(RG_xi, RG_eta)

#function finds the closest dust pixel
def find_closest(star, dust):
	idx,d2d,d3d=star.match_to_catalog_sky(dust)
	return idx

MS_Av=dust_Av[find_closest(MS_coords, dust_coords)]
MS_Av_err=dust_err[find_closest(MS_coords, dust_coords)]
AGy_Av=dust_Av[find_closest(AGy_coords, dust_coords)]
AGy_Av_err=dust_err[find_closest(AGy_coords, dust_coords)]
AGo_Av=dust_Av[find_closest(AGo_coords, dust_coords)]
AGo_Av_err=dust_err[find_closest(AGo_coords, dust_coords)]
RG_Av=dust_Av[find_closest(RG_coords, dust_coords)]
RG_Av_err=dust_err[find_closest(RG_coords, dust_coords)]

np.savetxt('/Users/amandaquirk/Desktop/MS_dust.txt', np.c_[MS_xi, MS_eta, MS_Av, MS_Av_err], fmt='%1.16f', delimiter=' ', header='xi (kpc), eta (kpc), median Av, Av error')
np.savetxt('/Users/amandaquirk/Desktop/AGy_dust.txt', np.c_[AGy_xi, AGy_eta, AGy_Av, AGy_Av_err], fmt='%1.16f', delimiter=' ', header='xi (kpc), eta (kpc), median Av, Av error')
np.savetxt('/Users/amandaquirk/Desktop/AGo_dust.txt', np.c_[AGo_xi, AGo_eta, AGo_Av, AGo_Av_err], fmt='%1.16f', delimiter=' ', header='xi (kpc), eta (kpc), median Av, Av error')
np.savetxt('/Users/amandaquirk/Desktop/RG_dust.txt', np.c_[RG_xi, RG_eta, RG_Av, RG_Av_err], fmt='%1.16f', delimiter=' ', header='xi (kpc), eta (kpc), median Av, Av error')

#function to calculate the weights
# def weights(err):
# 	return 1 / (err**2)

# dust_weight=weights(dust_err)

#function does the weighted meean
# def weighted_mean(data,norm_w):
# 	return sum([a*b for a,b in zip(data, norm_w)])

# #function to normalize the weights
# def normed_weight(w):
# 	sum_weights=sum(w)
# 	return w / sum_weights

# #matching star location with dust location
# AGy_Av=[]
# AGy_Av_err=[]
# for i in range(len(AGy_coords)):
# 	#picks out the first point
# 	c1=AGy_coords[i]
# 	Av=[]
# 	Av_err=[]
# 	weights=[]
# 	for j in range(len(dust_coords)):
# 		#below will go through all points to determine if the points are within the smoothing circle
# 		c2=dust_coords[j]
# 		sep=c1.separation(c2)
# 	if sep.arcsecond<0.01:
# 		Av.append(dust_Av[j])
# 		Av_err.append(dust_err[j])
# 		weights.append(dust_weight[j])
# 	if len(Av)>1:
# 		normed_weights=normed_weight(weights)
# 		AGy_Av.append(weighted_mean(Av, normed_weights))
# 		AGy_Av_err.append(np.median(Av_err))
# 	if len(Av)==1:
# 		AGy_Av.append(Av[0])
# 		AGy_Av_err.append(Av_err[0])

# print(len(AGy_xi))
# print(len(AGy_coords))
# print(len(AGy_Av))
# print(len(AGy_Av_err))





