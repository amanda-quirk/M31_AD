import numpy as np 
from astropy.io import fits
from astropy.table import Table 
import matplotlib.pyplot as plt 

#purpose is to take the error values that are negative and replace them with the mean error of the stars the have a similar magnitude

#read in the data and convert to a fits table
#hdulist=fits.open('/Users/amandaquirk/Documents/AsymmetricDrift/Data1/submasterSPLASH_quirk_master.fits', memmap=True)

hdulist=fits.open('/Users/amandaquirk/Desktop/submasterSPLASH_quirk_master.fits', memmap=True)

#puts data into an accessible table
data=Table(hdulist[1].data)

#extract the data from the table
F160W=data['F160W']
F110W=data['F110W']
age=data['Age'] #age tag: 1-4, 1 for MS and 4 for RG
error=data['verr'] #velocity error

age=age.astype(np.float)
error=error.astype(np.float)
F110W=F110W.astype(np.float)
F160W=F160W.astype(np.float)

F160_nonnan=[]
F110_nonnan=[]
for i in range(len(age)):
	if np.isnan(F160W[i])==False:
		F160_nonnan.append(F160W[i])
	if np.isnan(F110W[i])==False:
		F110_nonnan.append(F110W[i])

# print(len(F110_nonnan), len(F160_nonnan))

# plt.hist(F110_nonnan, 10, (15,25))
# plt.show()

good_data=[]
for i in range(len(age)):
	if error[i]>0:#np.isnan(F160W[i])==False and np.isnan(F110W[i])==False and error[i]>0:
		good_data.append(i)

print(len(good_data))

#first thing I need to do is to bin the magnitudes and find a mean error for every bin
#then i will assign the mean error to the star based on which bin they fit in

bins=np.linspace(15,24,10) #these are the divisions for both colors
#if mag=/= nan and err>0: divide into these bins, remove stars that have both bad error and nan mags
