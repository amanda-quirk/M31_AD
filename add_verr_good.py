import numpy as np
from astropy.io import fits 
from astropy.table import Table 
import matplotlib.pyplot as plt

# hdulistC=fits.open('/Users/amandaquirk/Documents/AsymmetricDrift/Data1/binned_master.fits', memmap=True)
# dataC=Table(hdulistC[1].data)

# verrC=dataC['VERR']
# IDC=dataC['ID']
# v_c=dataC['V']

# c1=fits.Column(name='ID', array=IDC, format='16A')
# c2=fits.Column(name='VERR', array=verrC, format='D')
# c3=fits.Column(name='V', array=v_c, format='D')
# t= fits.BinTableHDU.from_columns([c1, c2, c3])
# t.writeto('/Users/amandaquirk/Documents/AsymmetricDrift/Data1/real_Claire_IDs.fits')

hdulistC1=fits.open('/Users/amandaquirk/Documents/AsymmetricDrift/Data1/real_Claire_IDs.fits', memmap=True)
dataC1=Table(hdulistC1[1].data)

verrC1=dataC1['VERR']
IDC1=dataC1['ID']
v_c1=dataC1['V']

hdulist=fits.open('/Users/amandaquirk/Documents/AsymmetricDrift/Data1/submasterSPLASH_quirk_color_fixed.fits', memmap=True)
data=Table(hdulist[1].data)
ID=data['ID']

good_err=[]
good_v=[]
for i in range(len(ID)):
	n=np.where(ID[i]==IDC1)
	N=n[0]
	if len(N)==1:
		good_err.append(verrC1[N])
		good_v.append(v_c1[N])
	if len(N)==2:
		a=N[0]
		b=N[1]
		errs=[verrC1[a], verrC1[b]]
		vels=[v_c1[a], v_c1[b]]
		good_err.append(np.median(errs))
		good_v.append(np.median(vels))

print(good_v)

# plt.hist(good_err, 10, (0,50))
# plt.show()
