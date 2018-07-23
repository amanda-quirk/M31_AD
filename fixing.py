import numpy as np
from astropy.io import fits 
from astropy.table import Table 

#read in the data and convert to a fits table
hdulist=fits.open('/Users/amandaquirk/Documents/AsymmetricDrift/Data/submasterSPLASH_quirk_color.fits', memmap=True)

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
F160W_stars=data['F160W']
F110W_stars=data['F110W']

#convert things to floats
# v_s=v_s.astype(np.float)
# v_HImain=v_HImain.astype(np.float)
# v_HIclose=v_HIclose.astype(np.float)
# n=n.astype(np.float)
# age=age.astype(np.float)
# error=error.astype(np.float)
# index=index.astype(np.float)

data1 = Table.read('/Users/amandaquirk/Documents/AsymmetricDrift/Data/velocityhi-coords-amanda.dat', format='ascii')
#extract the data from the table
RA1=data1['RA'] #coordinate
Dec1=data1['DEC'] #coordinate
v_s1=data1['vopt'] #LOS velocity of star
v_HImain1=data1['vHIvfmain'] #velocity of HI main component 
v_HIclose1=data1['vHIclosestvopt'] #velocity of HI component closest to the star
n_comp1=data1['nbcomponentHI'] #number of peaks in HI velocity spectrum: 1-5
Age1=data1['age'] #age tag: 1-4, 1 for MS and 4 for RG

realv=[]
for i in range(len(RA)):
	velocities=[]
	for j in range(len(RA1)):
		if v_HIclose[i]==v_HIclose1[j] and v_HImain[i]==v_HImain1[j] and RA[i]==RA1[j] and Dec[i]==Dec1[j] and n[i]==n_comp1[j] and age[i]==Age1[j]:
			velocities.append(v_s1[j])
	if len(velocities)>0:
		realv.append(np.mean(velocities))

#writes data to new file
c1=fits.Column(name='RA', array=RA, format='11A')
c2=fits.Column(name='Dec', array=Dec, format='11A')
c3=fits.Column(name='v_s', array=realv, format='D')
c4=fits.Column(name='Age', array=age, format='K')
c5=fits.Column(name='HI_main', array=v_HImain, format='D')
c6=fits.Column(name='HI_close', array=v_HIclose, format='D')
c7=fits.Column(name='N', array=n, format='K')
c8=fits.Column(name='verr', array=error, format='11A')
c9=fits.Column(name='F160W', array=F160W_stars, format='11A')
c10=fits.Column(name='F110W', array=F110W_stars, format='11A')
c11=fits.Column(name='ID', array=ID, format='16A')
c12=fits.Column(name='Index in submasterSPLASH', array=index, format='K')
t= fits.BinTableHDU.from_columns([c1, c2, c3,c4,c5,c6,c7,c8,c9,c10,c11,c12])
t.writeto('/Users/amandaquirk/Desktop/submasterSPLASH_quirk_color.fits')