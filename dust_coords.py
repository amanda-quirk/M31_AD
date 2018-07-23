from astropy.wcs import wcs
import numpy as np 
from astropy.io import fits
import matplotlib.pyplot as plt

#read in dust map fits file
hdulist=fits.open('/Users/amandaquirk/Documents/AsymmetricDrift/phat_dust_maps/phat_dalcanton.AV.fits')
hdulist_err=fits.open('/Users/amandaquirk/Documents/AsymmetricDrift/phat_dust_maps/phat_dalcanton.AVerr.fits')

#read in header info for astropy's coordinate converter
w=wcs.WCS(hdulist[0].header)

#read in Av value
data=hdulist[0].data
data_err=hdulist_err[0].data #errors

#pixel values of dust map
x_vals=np.linspace(0,1117,1118)
y_vals=np.linspace(0,1369, 1370)

dust_av=[]
dust_err=[]
RA=[]
Dec=[]
for i in range(len(x_vals)):
	for j in range(len(y_vals)):
		x=int(x_vals[i]) #gets x coord of pixel
		y=int(y_vals[j]) #gets y coord of pixel
		dust_av.append(data[y,x]) #adds dust vales
		dust_err.append(data_err[y,x]) #adds the dust value errors
		coord=np.array([[x, y]], np.float_) #puts coords into an array
		world=w.wcs_pix2world(coord,1) #converts the coordinate into world coordinate
		RA.append(world[0][0])
		Dec.append(world[0][1])

#eliminating chunks of data I know are off M31

dust_av_good=[]
dust_err_good=[]
RA_good=[]
Dec_good=[]

for i in range(len(RA)):
	if dust_err[i]>0:
		dust_av_good.append(dust_av[i])
		dust_err_good.append(dust_err[i])
		RA_good.append(RA[i])
		Dec_good.append(Dec[i])

print(len(dust_av), len(dust_av_good))
plt.scatter(RA_good, Dec_good, c=dust_av_good)
plt.show()
file=open('/Users/amandaquirk/Desktop/dust_coords_medAv.txt', 'w')
file.write('#RA (deg), Dec (deg), median Av (mag), Av error\n')
for i in range(len(dust_av_good)):
	file.write('{} {} {} {}\n'.format(RA_good[i],Dec_good[i], dust_av_good[i], dust_err_good[i]))
file.close()