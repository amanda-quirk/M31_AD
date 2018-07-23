import matplotlib.pyplot as plt 
import numpy as np 
from functions import *

#importing all of the data
#xi (kpc), eta (kpc), average v(km/s), adjusted v err, var, n, HI main, HI close, ID, orginal index <-- header of data file to be read in 
MS_xi, MS_eta, MS_v, MS_HImain, MS_HIclose=np.loadtxt('/Users/amandaquirk/Desktop/MS_smoothed_chemin_no_outliers.txt', usecols=(0,1,2,6,7,), unpack=True)
AGy_xi, AGy_eta, AGy_v, AGy_HImain, AGy_HIclose=np.loadtxt('/Users/amandaquirk/Desktop/AGy_smoothed_chemin_no_outliers.txt', usecols=(0,1,2,6,7,), unpack=True)
AGo_xi, AGo_eta, AGo_v, AGo_HImain, AGo_HIclose=np.loadtxt('/Users/amandaquirk/Desktop/AGo_smoothed_chemin_no_outliers.txt', usecols=(0,1,2,6,7,), unpack=True)
RG_xi, RG_eta, RG_v, RG_HImain, RG_HIclose=np.loadtxt('/Users/amandaquirk/Desktop/RG_smoothed_chemin_no_outliers.txt', usecols=(0,1,2,6,7,), unpack=True)

#getting coordinates in shifted coordinate system (centered on M31)
MS_x=x_kpc(MS_xi, MS_eta) 
MS_y=y_kpc(MS_xi, MS_eta)

AGy_x=x_kpc(AGy_xi, AGy_eta) 
AGy_y=y_kpc(AGy_xi, AGy_eta)

AGo_x=x_kpc(AGo_xi, AGo_eta) 
AGo_y=y_kpc(AGo_xi, AGo_eta)

RG_x=x_kpc(RG_xi, RG_eta) 
RG_y=y_kpc(RG_xi, RG_eta)

#below will contain the deprojected distances of each star and gas cloud
MS_r=distance(MS_x,MS_y)
AGy_r=distance(AGy_x,AGy_y)
AGo_r=distance(AGo_x,AGo_y)
RG_r=distance(RG_x,RG_y)

#assign a PA to each star
MS_PA=[]
for i in range(len(MS_r)):
	MS_PA.append(PA(MS_x[i], MS_y[i]))
AGy_PA=[]
for i in range(len(AGy_r)):
	AGy_PA.append(PA(AGy_x[i], AGy_y[i]))
AGo_PA=[]
for i in range(len(AGo_r)):
	AGo_PA.append(PA(AGo_x[i], AGo_y[i]))
RG_PA=[]
for i in range(len(RG_r)):
	RG_PA.append(PA(RG_x[i], RG_y[i]))

#reads in the HI data
HI_r, HI_PA, HI_i, HI_v=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/HI_PA_i_vrot.txt', unpack=True)
def find_nearest_ring(radius):
	idx=[]
	for i in range(len(HI_r)):
		idx.append(np.abs(HI_r[i]-radius))
		n=np.argmin(idx)
	return n #index of the radius closest to the star's radius

#assigning a PA  and i to each star and HI component
MS_PA_ring=[]
MS_i=[]
for j in range(len(MS_r)):
	N=find_nearest_ring(MS_r[j])
	MS_PA_ring.append(HI_PA[N])
	MS_i.append(HI_i[N])

RG_PA_ring=[]
RG_i=[]
for j in range(len(RG_r)):
	N=find_nearest_ring(RG_r[j])
	RG_PA_ring.append(HI_PA[N])
	RG_i.append(HI_i[N])

AGy_PA_ring=[]
AGy_i=[]
for j in range(len(AGy_r)):
	N=find_nearest_ring(AGy_r[j])
	AGy_PA_ring.append(HI_PA[N])
	AGy_i.append(HI_i[N])

AGo_PA_ring=[]
AGo_i=[]
for j in range(len(AGo_r)):
	N=find_nearest_ring(AGo_r[j])
	AGo_PA_ring.append(HI_PA[N])
	AGo_i.append(HI_i[N])

#smoothed v and individual PA with tilted ring model
MS_si_vrot_tilt=Vrot_tilted_ring(MS_v,MS_PA_ring,MS_PA, MS_i)
RG_si_vrot_tilt=Vrot_tilted_ring(RG_v,RG_PA_ring,RG_PA, RG_i)
AGo_si_vrot_tilt=Vrot_tilted_ring(AGo_v,AGo_PA_ring,AGo_PA, AGo_i)
AGy_si_vrot_tilt=Vrot_tilted_ring(AGy_v,AGy_PA_ring,AGy_PA, AGy_i)

#HI with model
MS_HImain_vrot_tilt=Vrot_tilted_ring(MS_HImain,MS_PA_ring,MS_PA, MS_i)
MS_HIclose_vrot_tilt=Vrot_tilted_ring(MS_HIclose,MS_PA_ring,MS_PA, MS_i)
AGy_HImain_vrot_tilt=Vrot_tilted_ring(AGy_HImain,AGy_PA_ring,AGy_PA, AGy_i)
#AGy_HIclose_vrot_tilt=Vrot_tilted_ring(AGy_HIclose,AGy_PA_ring,AGy_PA, AGy_i)
AGo_HImain_vrot_tilt=Vrot_tilted_ring(AGo_HImain,AGo_PA_ring,AGo_PA, AGo_i)
#AGo_HIclose_vrot_tilt=Vrot_tilted_ring(AGo_HIclose,AGo_PA_ring,AGo_PA, AGo_i)
RG_HImain_vrot_tilt=Vrot_tilted_ring(RG_HImain,RG_PA_ring,RG_PA, RG_i)
#RG_HIclose_vrot_tilt=Vrot_tilted_ring(RG_HIclose,RG_PA_ring,RG_PA, RG_i)

#smooth v and individual PA with model
plt.scatter(MS_r, MS_si_vrot_tilt, s=1, c='b')
plt.scatter(MS_r, MS_HImain_vrot_tilt, s=1, c='k')
plt.xlim(4,20)
plt.ylim(100,300)
plt.xlabel('r [kpc]')
plt.ylabel('rotational v [km/s]')
plt.savefig('/Users/amandaquirk/Desktop/MS_no_outliers_vrot.png')
plt.close()

plt.scatter(RG_r, RG_si_vrot_tilt, s=1, c='r')
plt.scatter(RG_r, RG_HImain_vrot_tilt, s=1, c='k')
plt.xlim(4,20)
plt.ylim(100,300)
plt.xlabel('r [kpc]')
plt.ylabel('rotational v [km/s]')
plt.savefig('/Users/amandaquirk/Desktop/RG_no_outliers_vrot.png')
plt.close()

plt.scatter(AGo_r, AGo_si_vrot_tilt, s=1, c='m')
plt.scatter(AGo_r, AGo_HImain_vrot_tilt, s=1, c='k')
plt.xlim(4,20)
plt.ylim(100,300)
plt.xlabel('r [kpc]')
plt.ylabel('rotational v [km/s]')
plt.savefig('/Users/amandaquirk/Desktop/AGo_no_outliers_vrot.png')
plt.close()

plt.scatter(AGy_r, AGy_si_vrot_tilt, s=1, c='m')
plt.scatter(AGy_r, AGy_HImain_vrot_tilt, s=1, c='k')
plt.xlim(4,20)
plt.ylim(100,300)
plt.xlabel('r [kpc]')
plt.ylabel('rotational v [km/s]')
plt.savefig('/Users/amandaquirk/Desktop/AGy_no_outliers_vrot.png')
plt.close()

#writing data to files
file=open('/Users/amandaquirk/Desktop/MS_no_outliers_vrot.txt', 'w')
file.write('#r (kpc), v_rot tr model, HI main tr\n')
for i in range(len(MS_xi)):
	file.write('{} {} {}\n'.format(MS_r[i],MS_si_vrot_tilt[i], MS_HImain_vrot_tilt[i]))
file.close()

file=open('/Users/amandaquirk/Desktop/AGy_no_outliers_vrot.txt', 'w')
file.write('#r (kpc), v_rot tr model, HI main tr\n')
for i in range(len(AGy_xi)):
	file.write('{} {} {}\n'.format(AGy_r[i],AGy_si_vrot_tilt[i], AGy_HImain_vrot_tilt[i]))
file.close()

file=open('/Users/amandaquirk/Desktop/AGo_no_outliers_vrot.txt', 'w')
file.write('#r (kpc), v_rot tr model, HI main tr\n')
for i in range(len(AGo_xi)):
	file.write('{} {} {}\n'.format(AGo_r[i],AGo_si_vrot_tilt[i], AGo_HImain_vrot_tilt[i]))
file.close()

file=open('/Users/amandaquirk/Desktop/RG_no_outliers_vrot.txt', 'w')
file.write('#r (kpc), v_rot tr model, HI main tr\n')
for i in range(len(RG_xi)):
	file.write('{} {} {}\n'.format(RG_r[i],RG_si_vrot_tilt[i], RG_HImain_vrot_tilt[i]))
file.close()

