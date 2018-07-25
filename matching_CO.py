from astropy import units as u
from astropy.io import fits
from astropy.wcs import wcs
from astropy.coordinates import SkyCoord
from astropy.coordinates import match_coordinates_sky
import matplotlib.pyplot as plt 
import numpy as np 

#first step is to get the CO values
hdulist=fits.open('/Users/amandaquirk/Documents/AsymmetricDrift/Data/m31_iram/m31_iram.mom1.fits')

#read in header info for astropy's coordinate converter
w = wcs.WCS(hdulist[0].header, naxis = 2)

#read in Av value
data=hdulist[0].data 

#pixel values of dust map
x_vals=np.linspace(0, 549, 550)
y_vals=np.linspace(0, 649, 650)

CO_v=[]
RA=[]
Dec=[]
for i in range(len(x_vals)):
	for j in range(len(y_vals)):
		x=int(x_vals[i]) #gets x coord of pixel
		y=int(y_vals[j]) #gets y coord of pixel
		CO_v.append(data[y,x]) #adds the velocity vales
		coord=np.array([[x, y]], np.float_) #puts coords into an array
		world=w.wcs_pix2world(coord, 1) #converts the coordinate into world coordinate
		RA.append(world[0][0])
		Dec.append(world[0][1])

RA_good = [a for a, b in zip(RA, CO_v) if np.isnan(b)==False]
Dec_good = [a for a, b in zip(Dec, CO_v) if np.isnan(b)==False]
CO_v_good = [b for b in zip(CO_v) if np.isnan(b)==False]

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

CO_xi=m31_coord(RA_good, Dec_good, x=True)
CO_eta=m31_coord(RA_good, Dec_good, x=False)

#importing the star data
#xi (kpc), eta (kpc), average v(km/s), v err,var, n, HI main, HI close <-- header of data file to be read in 
MS_xi, MS_eta=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/MS_smoothed_chemin.txt', usecols=(0,1,), unpack=True)

#line to see where CO data ends
def y(x):
	return 16.8 - 0.8 * x

x = np.linspace(0,13, 10)

plt.scatter(CO_xi, CO_eta, c=CO_v_good)
im = plt.scatter(MS_xi, MS_eta, c= 'b', alpha=0.4)
plt.plot(x, y(x))
plt.gca().invert_xaxis()
plt.xlabel('xi (kpc)')
plt.ylabel('eta (kpc)')
plt.show()
plt.close()

#truncating the star data to just inside the CO field
CO_field = MS_eta < 16.8 - 0.8 * MS_xi
MS_xi = MS_xi[CO_field]
MS_eta = MS_eta[CO_field]

#turning xi and eta into a coordinate pair
def coords(xi, eta):
	return SkyCoord(xi/13.67, eta/13.67, unit=(u.deg, u.deg)) #dividing by 13.67 to convert to degrees

CO_coords=coords(CO_xi, CO_eta)
MS_coords=coords(MS_xi, MS_eta)

#function finds the closest CO
def find_closest(star, data):
	idx,d2d,d3d=star.match_to_catalog_sky(data)
	return idx

indices = find_closest(MS_coords, CO_coords)

#write oMSinial coordinates of CO
def CO_coords_deg(RA, Dec, x=True):
	c = SkyCoord(RA, Dec, unit=(u.deg, u.deg))
	if x==True:
		return c.ra.deg 
	if x==False:
		return c.dec.deg

CO_RA = CO_coords_deg(RA_good, Dec_good, x=True) 
CO_Dec = CO_coords_deg(RA_good, Dec_good, x=False)


MS_CO_v = []
MS_RA = []
MS_Dec = []
for i in indices:
	MS_CO_v.append(CO_v_good[i])
	MS_RA.append(CO_RA[i])
	MS_Dec.append(CO_Dec[i])

plt.scatter(MS_RA, MS_Dec, c=MS_CO_v)
plt.gca().invert_xaxis()
plt.xlabel('xi (kpc)')
plt.ylabel('eta (kpc)')
plt.colorbar()
plt.show()

np.savetxt('/Users/amandaquirk/Desktop/MS_CO.txt', np.c_[MS_xi, MS_eta, MS_CO_v, MS_RA, MS_Dec], fmt='%1.16f', delimiter=' ', header='xi (kpc), eta (kpc), CO v, RA (deg), Dec (deg)')
