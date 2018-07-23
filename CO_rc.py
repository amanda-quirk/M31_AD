import numpy as np 
import matplotlib.pyplot as plt 
from astropy import units as u
from astropy.coordinates import SkyCoord

#r (kpc), v_rot no model (km/s), v_rot tr model, smoothed PA, sPA v_rot no model, sPA v_rot tr model, n, HI main no model, HI main tr, HI close no model, HI close tr
RG_r, RG_vrot, HI_vrot = np.loadtxt('/Users/amandaquirk/Desktop/RG_master_vrot.txt', usecols=(0,2,8), unpack=True)
RG_xi, RG_eta=np.loadtxt('/Users/amandaquirk/Desktop/RG_smoothed_chemin.txt', usecols=(0,1,), unpack=True)
CO_xi, CO_eta, CO_v = np.loadtxt('/Users/amandaquirk/Desktop/CO/RG_CO.txt', unpack = True)

#need to do a cut on the location of the star data so that it doesn't go beyond the CO field
CO_field = RG_eta < 16.8 - 0.8 * RG_xi

RG_r = RG_r[CO_field]
RG_vrot = RG_vrot[CO_field]
HI_vrot = HI_vrot[CO_field]

#below defines a function that converts eta and xi and then uses the shifted coordinates to find the deprojected radius. to calculate the deprojected radius from the center, we must shift the xi and eta according to the PA of the disk. (This is a simple shift of coordinate axes.) We will assume the disk has a PA of 37 deg and an inclination of 77 deg 
def x(xi, eta):
	xi_deg=xi/13.67
	eta_deg=eta/13.67
	sine=np.sin(float(37*np.pi) / 180)
	cosine=np.cos(float(37*np.pi) / 180)
	x=(xi_deg*cosine)-(eta_deg*sine)
	return x 

def y(xi, eta):
	xi_deg=xi/13.67
	eta_deg=eta/13.67
	sine=np.sin(float(37*np.pi) / 180)
	cosine=np.cos(float(37*np.pi) / 180)
	y=(eta_deg*cosine)+(xi_deg*sine)
	return y

CO_x = x(CO_xi, CO_eta)
CO_y = y(CO_xi, CO_eta)

def distance(x, y):
	inclination_factor=np.cos(float(77*np.pi) / 180)**2
	ang_distance_sq=[(a**2)+(float(b**2)/inclination_factor) for a,b in zip(y,x)]
	ang_dist=[np.sqrt(a) for a in ang_distance_sq]
	dist=[a*13.67 for a in ang_dist]
	return dist

CO_r = distance(CO_x, CO_y)

def PA(x,y): #defined as west of north as positioned in maps (inverted x axis)
	if x>0:
		rad=np.arctan(float(y)/x)
		deg=90-(float(rad*180)/np.pi)
	else:
		rad=np.arctan(float(y)/x)
		deg=270-(float(rad*180)/np.pi)
	return deg + 37

CO_PA=[]
for i in range(len(CO_r)):
	CO_PA.append(PA(CO_x[i], CO_y[i]))

HI_r, HI_PA, HI_i, HI_v=np.loadtxt('/Users/amandaquirk/Desktop/HI_PA_i_vrot.txt', unpack=True)

#below defines a function to determine which HI ring a star is cloest to
def find_nearest_ring(radius):
	idx=[]
	for i in range(len(HI_r)):
		idx.append(np.abs(HI_r[i]-radius))
		n=np.argmin(idx)
	return n #index of the radius closest to the star's radius

#assigning a PA  and i to each star and HI component
CO_PA_ring=[]
CO_i=[]
for j in range(len(CO_r)):
	N=find_nearest_ring(CO_r[j])
	CO_PA_ring.append(HI_PA[N])
	CO_i.append(HI_i[N])

def Vrot_tilted_ring(v,PA_ring,PA_star, i_ring): 
	vsys= -300 #km/s, as defined in Claire's thesis
	A=[float(a-vsys)/np.sin(float(b*np.pi) / 180) for a, b in zip(v, i_ring)]
	B= [(np.tan(float((a-b)*np.pi) / 180))**2 for a,b in zip(PA_ring, PA_star)]
	C= [(np.cos(float(a*np.pi) / 180))**2 for a in i_ring]
	rotation_velocity= [a*np.sqrt(1+(float(b)/c)) for a,b,c in zip(A,B,C)]
	positive=[np.absolute(a) for a in rotation_velocity]
	return positive

CO_vrot = Vrot_tilted_ring(CO_v, CO_PA_ring, CO_PA, CO_i)

plt.scatter(RG_r, HI_vrot, s=5, c='darkgrey', alpha=0.4)
plt.scatter(CO_r, CO_vrot, s=5, c='darkcyan', marker='s', alpha=0.4)
plt.scatter(RG_r, RG_vrot, s=5, c='r', alpha=0.4)
plt.xlim(4,18)
plt.ylim(100,300)
plt.xlabel('r [kpc]')
plt.ylabel('rotational v [km/s]')
plt.savefig('RG_CO_rc.png')
plt.close()

def va(gas, stars):
	return gas - stars

RG_CO_ad = va(CO_vrot, RG_vrot)
RG_HI_ad = va(HI_vrot, RG_vrot)
CO_med = round(np.median(RG_CO_ad), 2)
HI_med = round(np.median(RG_HI_ad), 2)

plt.hist(RG_CO_ad, histtype='step', bins=range(-200, 200, 15), label='CO med= {}'.format(CO_med), linewidth=1.6,linestyle='--',stacked=True,fill=False, color='darkcyan')
plt.hist(RG_HI_ad, histtype='step', bins=range(-200, 200, 15),label='HI med= {}'.format(HI_med),linewidth=1.6,stacked=True,fill=False, color='darkgrey')
plt.legend(loc=2, frameon=False)
plt.xlabel('Asymmetric Drift (km/s)')
plt.xlim(-200, 200)
plt.savefig('RG_CO_ad.png')
