import matplotlib.pyplot as plt 
import numpy as np 
from matplotlib import rc,rcParams
from shapely.geometry import Point 
from shapely.geometry.polygon import Polygon
from functions import color, x_kpc, y_kpc, distance, PA, Vrot_tilted_ring

F336, F475, F814=np.loadtxt('/Users/amandaquirk/Desktop/MS_smooth_color.txt', unpack=True)
xi, eta, v, n, HImain, HIclose=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/MS_smoothed_chemin.txt', usecols=(0,1,2,5,6,7,), unpack=True)

#divide stars into 3 groups based on color

color=color(F336, F475)

MS_xi=[]
MS_eta=[]
MS_v=[]
MS_HI=[]
BL1_xi=[]
BL1_eta=[]
BL1_v=[]
BL1_HI=[]
BL2_xi=[]
BL2_eta=[]
BL2_v=[]
BL2_HI=[]
for i in range(len(color)):
	if color[i]<=-0.2 and F475[i]<=23.5:
		MS_xi.append(xi[i])
		MS_eta.append(eta[i])
		MS_v.append(v[i])
		MS_HI.append(HIclose[i])
	elif color[i]>=-0.2 and color[i]<0.2 and F475[i]<=23.5:
		BL1_xi.append(xi[i])
		BL1_eta.append(eta[i])
		BL1_v.append(v[i])
		BL1_HI.append(HIclose[i]) 
	elif color[i]<=1.4 and color[i]>0.2 and F475[i]<=23.5:
		BL2_xi.append(xi[i])
		BL2_eta.append(eta[i])
		BL2_v.append(v[i])
		BL2_HI.append(HIclose[i])

MS_x=x_kpc(MS_xi, MS_eta)
MS_y=y_kpc(MS_xi, MS_eta)
BL1_x=x_kpc(BL1_xi, BL1_eta)
BL1_y=y_kpc(BL1_xi, BL1_eta)
BL2_x=x_kpc(BL2_xi, BL2_eta)
BL2_y=y_kpc(BL2_xi, BL2_eta)

MS_r=distance(MS_x, MS_y)
BL1_r=distance(BL1_x, BL1_y)
BL2_r=distance(BL2_x, BL2_y)

MS_PA=[]
for i in range(len(MS_r)):
	MS_PA.append(PA(MS_x[i], MS_y[i]))
BL1_PA=[]
for i in range(len(BL1_r)):
	BL1_PA.append(PA(BL1_x[i], BL1_y[i]))
BL2_PA=[]
for i in range(len(BL2_r)):
	BL2_PA.append(PA(BL2_x[i], BL2_y[i]))

#reads in the HI data
HI_r, HI_PA, HI_i, HI_v=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/HI_PA_i_vrot.txt', unpack=True)

#below defines a function to determine which HI ring a star is cloest to
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

BL1_PA_ring=[]
BL1_i=[]
for j in range(len(BL1_r)):
	N=find_nearest_ring(BL1_r[j])
	BL1_PA_ring.append(HI_PA[N])
	BL1_i.append(HI_i[N])

BL2_PA_ring=[]
BL2_i=[]
for j in range(len(BL2_r)):
	N=find_nearest_ring(BL2_r[j])
	BL2_PA_ring.append(HI_PA[N])
	BL2_i.append(HI_i[N])

MS_vrot=Vrot_tilted_ring(MS_v,MS_PA_ring,MS_PA, MS_i)
BL1_vrot=Vrot_tilted_ring(BL1_v,BL1_PA_ring,BL1_PA, BL1_i)
BL2_vrot=Vrot_tilted_ring(BL2_v,BL2_PA_ring,BL2_PA, BL2_i)

MS_HI_vrot=Vrot_tilted_ring(MS_HI,MS_PA_ring,MS_PA, MS_i)
BL1_HI_vrot=Vrot_tilted_ring(BL1_HI,BL1_PA_ring,BL1_PA, BL1_i)
BL2_HI_vrot=Vrot_tilted_ring(BL2_HI,BL2_PA_ring,BL2_PA, BL2_i)


def rotation_curve(r, star_vrot, HI_vrot, description):
	plt.scatter(r, star_vrot, s=1, c='b')
	plt.scatter(r, HI_vrot, s=1, c='k')
	plt.xlim(4,20)
	plt.ylim(100,300)
	plt.xlabel('r [kpc]')
	plt.ylabel('rotational v [km/s]')
	plt.savefig('/Users/amandaquirk/Desktop/{}_vrot.png'.format(description))
	plt.close()

rotation_curve(MS_r, MS_vrot, MS_HI_vrot, 'MS' )
rotation_curve(BL1_r, BL1_vrot, BL1_HI_vrot, 'BL1')
rotation_curve(BL2_r, BL2_vrot, BL2_HI_vrot, 'BL2')


