from astropy.io import fits
from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np 

hdulist= fits.open('/Users/amandaquirk/Desktop/Asymmetric Disk/Data/subMasterSPLASH.fits', memmap=True)

data=Table(hdulist[1].data)
ZSNR=data['ZSNR']
L=data['LIKELIHOOD']
flag=data['flags']

#print(np.nanmean(L))

plt.hist(ZSNR[~np.isnan(ZSNR)], bins=15, norm=True)
plt.xlabel('ZSNR')
plt.savefig('/Users/amandaquirk/Desktop/ZSNRparameter.png')

print(len(L))
print(len(L[~np.isnan(L)]))
plt.hist(L[~np.isnan(L)], bins=15, norm=True)
plt.savefig('/Users/amandaquirk/Desktop/Lparameter.png')
