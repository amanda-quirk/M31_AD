import numpy as np 
from astropy.io import fits
from astropy.table import Table 

hdulist=fits.open('/Users/amandaquirk/Documents/AsymmetricDrift/Data/submasterSPLASH.fits', memmap=True)
data=Table(hdulist[1].data)
ID=data['OBJNAME']

RGs_xi, RGs_eta, RGs_vs, RGs_err, RGs_dispersion, RGs_n, RGs_vHImain, RGs_vHIclose, RGs_index=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/RG_smoothed_chemin.txt', usecols=(0,1,2,3,4,5,6,7,9,), unpack=True)

RGs_ID=[]
for i in range(len(RGs_xi)):
	n=RGs_index[i]
	N=int(n)
	RGs_ID.append(ID[N])

file=open('/Users/amandaquirk/Documents/AsymmetricDrift/Data/RG_smoothed_chemin.txt', 'w')
file.write('#xi (kpc), eta (kpc), average v(km/s), v err, var, n, HI main, HI close, ID, orginal index\n')
for i in range(len(RGs_xi)):
	file.write('{} {} {} {} {} {} {} {} {} {}\n'.format(RGs_xi[i],RGs_eta[i], RGs_vs[i], RGs_err[i], RGs_dispersion[i], RGs_n[i], RGs_vHImain[i], RGs_vHIclose[i], RGs_ID[i], int(RGs_index[i])))
file.close()