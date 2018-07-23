import numpy as np
from astropy.io import fits
from astropy.table import Table
import matplotlib.pyplot as plt

data0=Table.read('/Users/amandaquirk/Documents/AsymmetricDrift/Data/submasterSPLASH.fits')
F160W=data0['F160W']
F110W=data0['F110W']

ind=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/RG_smoothed_chemin.txt', usecols=(9,), unpack=True)

RG_F110W=[]
RG_F160W=[]
for i in range(len(ind)):
	RG_F110W.append(F110W[int(ind[i])])
	RG_F160W.append(F160W[int(ind[i])])

np.savetxt('/Users/amandaquirk/Desktop/RG_color.txt', np.c_[RG_F110W, RG_F160W, ind], fmt='%1.16f', delimiter=' ', header='F110W, F160W, original ind')
# print(len(RA), len(Dec), len(v_s), len(HImain), len(HIclose), len(n_comp), len(Age), len(verr), len(ID))

# F160W_stars=[]
# F110W_stars=[]
# indices=[]
#for i in range(len(RA)):
#	for j in range(len(ID0)):
#		mag16=[]
#		mag10=[]
#		ind=[]	
#		if ID[i]==ID0[j]:
#			mag16.append(F160W[j])
#			mag10.append(F110W[j])
#			ind.append(j)
#	if len(mag16)>0:
#		F160W_stars.append(mag16[0])
#		F110W_stars.append(mag10[0])
#		indices.append(ind[0])
# for i in range(len(RA)):
# 	n=np.where(ID[i]==ID0)
# 	N=n[0]
# 	if len(N)==1:
# 		F160W_stars.append(F160W[N])
# 		F110W_stars.append(F110W[N])
# 		indices.append(N)
# 	if len(N)>1:
# 		a=N[0]
# 		F160W_stars.append(F160W[a])
# 		F110W_stars.append(F110W[a])
# 		indices.append(a)

# print(len(RA), len(F160W_stars))


# def color(color1, color2):
# 	return color1-color2

# RGB_color=[]
# RGB_mag=[]
# for i in range(len(Age)):
# 	if Age[i]==4:
# 		RGB_color.append(color(F110W_stars[i], F160W_stars[i]))
# 		RGB_mag.append(F160W_stars[i])

# plt.scatter(RGB_color, RGB_mag)
# plt.xlim(0,2)
# plt.ylim(17,24)
# plt.savefig('/Users/amandaquirk/Desktop/RGB_CMD.png')

# indices=np.loadtxt('/Users/amandaquirk/Desktop/indices.txt')
# F160W_stars=np.loadtxt('/Users/amandaquirk/Desktop/F160W.txt')
# F110W_stars=np.loadtxt('/Users/amandaquirk/Desktop/F110W.txt')

# #writes data to new file
# c1=fits.Column(name='RA', array=RA, format='11A')
# c2=fits.Column(name='Dec', array=Dec, format='11A')
# c3=fits.Column(name='v_s', array=v_s, format='D')
# c4=fits.Column(name='Age', array=Age, format='K')
# c5=fits.Column(name='HI_main', array=HImain, format='D')
# c6=fits.Column(name='HI_close', array=HIclose, format='D')
# c7=fits.Column(name='N', array=n_comp, format='K')
# c8=fits.Column(name='verr', array=verr, format='11A')
# c9=fits.Column(name='F160W', array=F160W_stars, format='11A')
# c10=fits.Column(name='F110W', array=F110W_stars, format='11A')
# c11=fits.Column(name='ID', array=ID, format='16A')
# c12=fits.Column(name='Index in submasterSPLASH', array=indices, format='K')
# t= fits.BinTableHDU.from_columns([c1, c2, c3,c4,c5,c6,c7,c8,c9,c10,c11,c12])
# t.writeto('/Users/amandaquirk/Desktop/submasterSPLASH_quirk_color.fits')

