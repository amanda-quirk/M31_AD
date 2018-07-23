from astropy.io import fits 
from astropy.table import Table 
import numpy as np

#read in the data and convert to a fits table
hdulist=fits.open('/Users/amandaquirk/Downloads/submasterSPLASH_quirk_master.fits', memmap=True)

#puts data into an accessible table
data=Table(hdulist[1].data)

#extract the data from the table
RA1=data['RA'] #coordinate
Dec1=data['Dec'] #coordinate
v_HImain1=data['HI_main'] #velocity of HI main component 
v_HIclose1=data['HI_close'] #velocity of HI component closest to the star
n1=data['N'] #number of peaks in HI velocity spectrum: 1-5
Age1=data['Age'] #age tag: 1-4, 1 for MS and 4 for RG
index1=data['Index in submasterSPLASH']
ID1=data['ID']
# F160W=data['F160W']
# F110W=data['F110W']
# verr=data['verr']
# v_s=data['v_s']
F160W1=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data1/F160W.txt')
F110W1=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data1/F110W.txt')
verr1=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data1/errors_quirk.txt')
v_s1=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data1/quirk_v.txt')

#removing stars with placeholder errors
RA=[]
Dec=[]
v_HImain=[]
v_HIclose=[]
n=[]
Age=[]
index=[]
ID=[]
F160W=[]
F110W=[]
verr=[]
v_s=[]
for i in range(len(verr1)):
	if verr1[i]>0:
		RA.append(RA1[i])
		Dec.append(Dec1[i])
		v_HImain.append(v_HImain1[i])
		v_HIclose.append(v_HIclose1[i])
		n.append(n1[i])
		Age.append(Age1[i])
		index.append(index1[i])
		ID.append(ID1[i])
		F160W.append(F160W1[i])
		F110W.append(F110W1[i])
		verr.append(verr1[i])
		v_s.append(v_s1[i])
print(len(v_s))
#writes data to new file
c1=fits.Column(name='RA', array=RA, format='11A')
c2=fits.Column(name='Dec', array=Dec, format='11A')
c3=fits.Column(name='v_s', array=v_s, format='D')
c4=fits.Column(name='Age', array=Age, format='K')
c5=fits.Column(name='HI_main', array=v_HImain, format='D')
c6=fits.Column(name='HI_close', array=v_HIclose, format='D')
c7=fits.Column(name='N', array=n, format='K')
c8=fits.Column(name='verr', array=verr, format='D')
c9=fits.Column(name='F160W', array=F160W, format='D')
c10=fits.Column(name='F110W', array=F110W, format='D')
c11=fits.Column(name='ID', array=ID, format='16A')
c12=fits.Column(name='Index in submasterSPLASH', array=index, format='K')
t= fits.BinTableHDU.from_columns([c1, c2, c3,c4,c5,c6,c7,c8,c9,c10,c11,c12])
t.writeto('/Users/amandaquirk/Desktop/submasterSPLASH_quirk_master.fits')
