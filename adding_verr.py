from astropy import units as u
from astropy.coordinates import SkyCoord
import numpy as np
from astropy.io import fits 
from astropy.table import Table 

# data0=Table.read('/Users/amandaquirk/Documents/AsymmetricDrift/Data/submasterSPLASH.fits')
#data from initial file
# RA0=data0['RA']
# Dec0=data0['DEC']
# name=data0['OBJNAME']
# z=data0['Z']

#my matched IDs from initial file
data_ID=Table.read('/Users/amandaquirk/Desktop/IDs.fits')
ID=data_ID['ID']

#read in the data and convert to a fits table
#stars i was using-- will be a subset at the end
data = Table.read('/Users/amandaquirk/Documents/AsymmetricDrift/Data/velocityhi-coords-amanda.dat', format='ascii')
#extract the data from the table
RA=data['RA'] #coordinate
Dec=data['DEC'] #coordinate
v_s=data['vopt'] #LOS velocity of star
v_HImain=data['vHIvfmain'] #velocity of HI main component 
v_HIclose=data['vHIclosestvopt'] #velocity of HI component closest to the star
n_comp=data['nbcomponentHI'] #number of peaks in HI velocity spectrum: 1-5
Age=data['age'] #age tag: 1-4, 1 for MS and 4 for RG

#below was used to get the ID of my stars
# ID=[]
# ind=[]
# for i in range(len(RA)):
# 	for j in range(len(RA0)):
# 		if RA0[j]==RA[i] and Dec0[j]==Dec[i] and abs(((3*10**5)*z[j])-v_s[i])<=.01:
# 			ID.append(name[j])
# 			ind.append(i)

# print(len(RA),len(ID), len(ind))
# c1=fits.Column(name='ID', array=ID, format='16A')
# c2=fits.Column(name='index', array=ind, format='K')
# t= fits.BinTableHDU.from_columns([c1, c2])
# t.writeto('/Users/amandaquirk/Desktop/IDs.fits')

#read in Claire's file
C_data= Table.read('/Users/amandaquirk/Desktop/ClaireIDs.fits')
#convert data to usable variables
err=C_data['err']
ID_C=C_data['ID']

#functions to compute the weighted stats
def normed_errors(error, sum_error):
	return error / sum_error

def weighted_mean(data,norm_error):
	return sum([a*b for a,b in zip(data, norm_error)])

#now I want to loop through all of the data and assign the appropriate velocity errors to the data i have
#v_err=[]
ID_unique=[]
RA_unique=[]
Dec_unique=[]
v_s_unique=[]
HImain_unique=[]
HIclose_unique=[]
n_unique=[]
agebin_unique=[]
check=[]
for i in range(len(ID)):
	if ID[i] in ID_C:
		n=np.where(ID[i]==ID_C)
		N=n[0]
		if len(N)==1: #good data
#			v_err.append(err[N])
			ID_unique.append(ID[i])
			RA_unique.append(RA[i])
			Dec_unique.append(Dec[i])
			v_s_unique.append(v_s[i])
			HImain_unique.append(v_HImain[i])
			HIclose_unique.append(v_HIclose[i])
			n_unique.append(n_comp[i])
			agebin_unique.append(Age[i])
			check.append(i) 
		if len(N)>2: #serendipitous stars
#			v_err.append(np.nan)
			ID_unique.append(ID[i])
			RA_unique.append(RA[i])
			Dec_unique.append(Dec[i])
			v_s_unique.append(v_s[i])
			HImain_unique.append(v_HImain[i])
			HIclose_unique.append(v_HIclose[i])
			n_unique.append(n_comp[i])
			agebin_unique.append(Age[i])
			check.append(i)
		#velocities=[]
		#errors_norm=[]
		#errors=[]
		if len(N)==2: #repeated stars-- want to combine into one
			velocities=[]
			errors_norm=[]
			errors=[]
			a=N[0]
			b=N[1]
			velocities.append(v_s[i])
			velocities.append(v_s[i])
			errors_norm.append(err[a]/(err[a]+err[b])) #will have normed errors
			errors_norm.append(err[b]/(err[a]+err[b]))
			errors.append(err[a])
			errors.append(err[b])
			ID_unique.append(ID[i])
			RA_unique.append(RA[i])
			Dec_unique.append(Dec[i])
			HImain_unique.append(v_HImain[i])
			HIclose_unique.append(v_HIclose[i])
			n_unique.append(n_comp[i])
			agebin_unique.append(Age[i])
			if len(velocities)==2:
				v_s_unique.append(weighted_mean(velocities,errors_norm))
#				v_err.append(np.median(errors))
				check.append(i)
	else: #the 151 stars i can't find for whatever reason!
#		v_err.append(np.nan)
		ID_unique.append(ID[i])
		RA_unique.append(RA[i])
		Dec_unique.append(Dec[i])
		v_s_unique.append(v_s[i])
		HImain_unique.append(v_HImain[i])
		HIclose_unique.append(v_HIclose[i])
		n_unique.append(n_comp[i])
		agebin_unique.append(Age[i])
		check.append(i)

	
#print(len(v_s_unique), len(v_err), len(RA_unique), len(check))

#file i formatted based on error assigments above
v_err=np.loadtxt('/Users/amandaquirk/Desktop/errors_quirk.txt')

v_err_final=[]
ID_final=[]
RA_final=[]
Dec_final=[]
v_s_final=[]
HImain_final=[]
HIclose_final=[]
age_final=[]
n_comp_final=[]
#remove nan data
for i in range(len(v_err)):
	if np.isnan(v_err[i])==False:
		v_err_final.append(v_err[i])
		ID_final.append(ID_unique[i])
		RA_final.append(RA_unique[i])
		Dec_final.append(Dec_unique[i])
		v_s_final.append(v_s_unique[i])
		HImain_final.append(HImain_unique[i])
		HIclose_final.append(HIclose_unique[i])
		age_final.append(agebin_unique[i])
		n_comp_final.append(n_comp[i])

#writes data to new file
c1=fits.Column(name='RA', array=RA_final, format='11A')
c2=fits.Column(name='Dec', array=Dec_final, format='11A')
c3=fits.Column(name='v_s', array=v_s_final, format='D')
c4=fits.Column(name='Age', array=age_final, format='K')
c5=fits.Column(name='HI_main', array=HImain_final, format='D')
c6=fits.Column(name='HI_close', array=HIclose_final, format='D')
c7=fits.Column(name='N', array=n_comp_final, format='K')
c8=fits.Column(name='v_err', array=v_err_final, format='11A')
c9=fits.Column(name='ID', array=ID_final, format='16A')
t= fits.BinTableHDU.from_columns([c1, c2, c3,c4,c5,c6,c7,c8,c9])
t.writeto('/Users/amandaquirk/Desktop/submasterSPLASH_quirk.fits')



