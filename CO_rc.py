import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib import rc 
from matplotlib.ticker import MaxNLocator
from functions import *

#r (kpc), v_rot no model (km/s), v_rot tr model, smoothed PA, sPA v_rot no model, sPA v_rot tr model, n, HI main no model, HI main tr, HI close no model, HI close tr
RG_r, RG_vrot, RG_HI_vrot = np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/RG_master_vrot.txt', usecols=(0,2,8), unpack=True)
RG_xi, RG_eta=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/RG_smoothed_chemin.txt', usecols=(0,1,), unpack=True)
RG_CO_xi, RG_CO_eta, RG_CO_v, RG_CO_RA, RG_CO_Dec = np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/RG_CO.txt', unpack = True)
MS_r, MS_vrot, MS_HI_vrot = np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/MS_master_vrot.txt', usecols=(0,2,8), unpack=True)
MS_xi, MS_eta=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/MS_smoothed_chemin.txt', usecols=(0,1,), unpack=True)
MS_CO_xi, MS_CO_eta, MS_CO_v, MS_CO_RA, MS_CO_Dec = np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/MS_CO.txt', unpack = True)
AGy_r, AGy_vrot, AGy_HI_vrot = np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/AGy_master_vrot.txt', usecols=(0,2,8), unpack=True)
AGy_xi, AGy_eta=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/AGy_smoothed_chemin.txt', usecols=(0,1,), unpack=True)
AGy_CO_xi, AGy_CO_eta, AGy_CO_v, AGy_CO_RA, AGy_CO_Dec = np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/AGy_CO.txt', unpack = True)
AGo_r, AGo_vrot, AGo_HI_vrot = np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/AGo_master_vrot.txt', usecols=(0,2,8), unpack=True)
AGo_xi, AGo_eta=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/AGo_smoothed_chemin.txt', usecols=(0,1,), unpack=True)
AGo_CO_xi, AGo_CO_eta, AGo_CO_v, AGo_CO_RA, AGo_CO_Dec = np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/AGo_CO.txt', unpack = True)

#need to do a cut on the location of the star data so that it doesn't go beyond the CO field
MS_CO_field = MS_eta < 16.8 - 0.8 * MS_xi
AGy_CO_field = AGy_eta < 16.8 - 0.8 * AGy_xi
AGo_CO_field = AGo_eta < 16.8 - 0.8 * AGo_xi
RG_CO_field = RG_eta < 16.8 - 0.8 * RG_xi

MS_r = MS_r[MS_CO_field]
MS_vrot = MS_vrot[MS_CO_field]
MS_HI_vrot = MS_HI_vrot[MS_CO_field]

AGy_r = AGy_r[AGy_CO_field]
AGy_vrot = AGy_vrot[AGy_CO_field]
AGy_HI_vrot = AGy_HI_vrot[AGy_CO_field]

AGo_r = AGo_r[AGo_CO_field]
AGo_vrot = AGo_vrot[AGo_CO_field]
AGo_HI_vrot = AGo_HI_vrot[AGo_CO_field]

RG_r = RG_r[RG_CO_field]
RG_vrot = RG_vrot[RG_CO_field]
HI_vrot = RG_HI_vrot[RG_CO_field]

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

MS_CO_x = x(MS_CO_xi, MS_CO_eta)
MS_CO_y = y(MS_CO_xi, MS_CO_eta)

AGy_CO_x = x(AGy_CO_xi, AGy_CO_eta)
AGy_CO_y = y(AGy_CO_xi, AGy_CO_eta)

AGo_CO_x = x(AGo_CO_xi, AGo_CO_eta)
AGo_CO_y = y(AGo_CO_xi, AGo_CO_eta)

RG_CO_x = x(RG_CO_xi, RG_CO_eta)
RG_CO_y = y(RG_CO_xi, RG_CO_eta)

def distance(x, y):
	inclination_factor=np.cos(float(77*np.pi) / 180)**2
	ang_distance_sq=[(a**2)+(float(b**2)/inclination_factor) for a,b in zip(y,x)]
	ang_dist=[np.sqrt(a) for a in ang_distance_sq]
	dist=[a*13.67 for a in ang_dist]
	return dist

MS_CO_r = distance(MS_CO_x, MS_CO_y)
AGy_CO_r = distance(AGy_CO_x, AGy_CO_y)
AGo_CO_r = distance(AGo_CO_x, AGo_CO_y)
RG_CO_r = distance(RG_CO_x, RG_CO_y)

def PA(x,y): #defined as west of north as positioned in maps (inverted x axis)
	if x>0:
		rad=np.arctan(float(y)/x)
		deg=90-(float(rad*180)/np.pi)
	else:
		rad=np.arctan(float(y)/x)
		deg=270-(float(rad*180)/np.pi)
	return deg + 37

MS_CO_PA=[]
for i in range(len(MS_CO_r)):
	MS_CO_PA.append(PA(MS_CO_x[i], MS_CO_y[i]))

AGy_CO_PA=[]
for i in range(len(AGy_CO_r)):
	AGy_CO_PA.append(PA(AGy_CO_x[i], AGy_CO_y[i]))

AGo_CO_PA=[]
for i in range(len(AGo_CO_r)):
	AGo_CO_PA.append(PA(AGo_CO_x[i], AGo_CO_y[i]))

RG_CO_PA=[]
for i in range(len(RG_CO_r)):
	RG_CO_PA.append(PA(RG_CO_x[i], RG_CO_y[i]))

HI_r, HI_PA, HI_i, HI_v=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/HI_PA_i_vrot.txt', unpack=True)

#below defines a function to determine which HI ring a star is cloest to
def find_nearest_ring(radius):
	idx=[]
	for i in range(len(HI_r)):
		idx.append(np.abs(HI_r[i]-radius))
		n=np.argmin(idx)
	return n #index of the radius closest to the star's radius

#assigning a PA  and i to each star and HI component
MS_CO_PA_ring=[]
MS_CO_i=[]
for j in range(len(MS_CO_r)):
	N=find_nearest_ring(MS_CO_r[j])
	MS_CO_PA_ring.append(HI_PA[N])
	MS_CO_i.append(HI_i[N])

AGy_CO_PA_ring=[]
AGy_CO_i=[]
for j in range(len(AGy_CO_r)):
	N=find_nearest_ring(AGy_CO_r[j])
	AGy_CO_PA_ring.append(HI_PA[N])
	AGy_CO_i.append(HI_i[N])

AGo_CO_PA_ring=[]
AGo_CO_i=[]
for j in range(len(AGo_CO_r)):
	N=find_nearest_ring(AGo_CO_r[j])
	AGo_CO_PA_ring.append(HI_PA[N])
	AGo_CO_i.append(HI_i[N])

RG_CO_PA_ring=[]
RG_CO_i=[]
for j in range(len(RG_CO_r)):
	N=find_nearest_ring(RG_CO_r[j])
	RG_CO_PA_ring.append(HI_PA[N])
	RG_CO_i.append(HI_i[N])

#convert the CO velocity from LSR to heliocentric
def LSR2Helio(v, RA, Dec): #km/s, deg, deg
	v_sun = 19.7 #km/s
	alpha_0 = 271 * np.pi / 180 #radians; B1950 coordinates
	delta_0 = 30 * np.pi / 180 #radians; B1950 coordinates

	#convert to this wack coordinate system
	def B1950_RA(RA, Dec):
		degrees = RA - 0.640265 + 0.278369 * np.sin(RA * np.pi / 180) * np.tan(Dec * np.pi /180)
		return degrees * np.pi / 180 #rad

	def B1950_Dec(RA, Dec):
		degrees = Dec - 0.278369 * np.cos(RA * np.pi / 180)
		return degrees * np.pi / 180 #rad

	alpha = B1950_RA(RA, Dec)
	delta = B1950_Dec(RA, Dec)

	trig = np.cos(alpha - alpha_0) * np.cos(delta) * np.cos(delta_0) + np.sin(delta) * np.sin(delta_0)

	return v - v_sun * trig 

MS_CO_v_helio = LSR2Helio(MS_CO_v, MS_CO_RA, MS_CO_Dec)
AGy_CO_v_helio = LSR2Helio(AGy_CO_v, AGy_CO_RA, AGy_CO_Dec)
AGo_CO_v_helio = LSR2Helio(AGo_CO_v, AGo_CO_RA, AGo_CO_Dec)
RG_CO_v_helio = LSR2Helio(RG_CO_v, RG_CO_RA, RG_CO_Dec)

def Vrot_tilted_ring(v,PA_ring,PA_star, i_ring): 
	vsys= -300 #km/s, as defined in Claire's thesis
	A=[float(a-vsys)/np.sin(float(b*np.pi) / 180) for a, b in zip(v, i_ring)]
	B= [(np.tan(float((a-b)*np.pi) / 180))**2 for a,b in zip(PA_ring, PA_star)]
	C= [(np.cos(float(a*np.pi) / 180))**2 for a in i_ring]
	rotation_velocity= [a*np.sqrt(1+(float(b)/c)) for a,b,c in zip(A,B,C)]
	positive=[np.absolute(a) for a in rotation_velocity]
	return positive

MS_CO_v = Vrot_tilted_ring(MS_CO_v_helio, MS_CO_PA_ring, MS_CO_PA, MS_CO_i)
AGy_CO_v = Vrot_tilted_ring(AGy_CO_v_helio, AGy_CO_PA_ring, AGy_CO_PA, AGy_CO_i)
AGo_CO_v = Vrot_tilted_ring(AGo_CO_v_helio, AGo_CO_PA_ring, AGo_CO_PA, AGo_CO_i)
RG_CO_v = Vrot_tilted_ring(RG_CO_v_helio, RG_CO_PA_ring, RG_CO_PA, RG_CO_i)

# plt.scatter(RG_r, HI_vrot, s=5, c='darkgrey', alpha=0.4)
# plt.scatter(CO_r, CO_vrot, s=5, c='darkcyan', marker='s', alpha=0.4)
# plt.scatter(RG_r, RG_vrot, s=5, c='r', alpha=0.4)
# plt.xlim(4,18)
# plt.ylim(100,300)
# plt.xlabel('r [kpc]')
# plt.ylabel('rotational v [km/s]')
# plt.savefig('/Users/amandaquirk/Desktop/RG_CO_rc.png')
# plt.close()

# from functions import *
# median_r = bins - delta_r / 2
# RG_vrot_med = median_line(RG_r, RG_vrot)
# MS_vrot_med = median_line(MS_r, MS_vrot)
# AGy_vrot_med = median_line(AGy_r, AGy_vrot)
# AGo_vrot_med = median_line(AGo_r, AGo_vrot)
# HI_RG_vrot_med = median_line(RG_r, RG_CO_v)
# HI_MS_vrot_med = median_line(MS_r, MS_CO_v)
# HI_AGy_vrot_med = median_line(AGy_r, AGy_CO_v)
# HI_AGo_vrot_med = median_line(AGo_r, AGo_CO_v)

# rc('font', family = 'serif')
# f, axes= plt.subplots(4,1, sharey=False, sharex=True, figsize=(4,9.8))
# axes[0].scatter(MS_r, MS_CO_v, s=2, c='darkcyan')
# axes[0].scatter(MS_r, MS_vrot, s=2, c='b', alpha=0.4)
# axes[1].scatter(AGy_r, AGy_CO_v, s=2, c='darkcyan')
# axes[1].scatter(AGy_r, AGy_vrot, s=2, c='m', alpha=0.4)
# axes[2].scatter(AGo_r, AGo_CO_v, s=2, c='darkcyan')
# axes[2].scatter(AGo_r, AGo_vrot, s=2, c='green', alpha=0.4)
# axes[3].scatter(RG_r, RG_CO_v, s=2, c='darkcyan')
# axes[3].scatter(RG_r, RG_vrot, s=2, c='r', alpha=0.4)
# axes[0].annotate('MS', xy=(9,115), horizontalalignment='left', fontsize=12)
# axes[1].annotate('young AGB', xy=(9,115), horizontalalignment='left', fontsize=12)
# axes[2].annotate('older AGB', xy=(9,115), horizontalalignment='left', fontsize=12)
# axes[3].annotate('RGB', xy=(9,115), horizontalalignment='left', fontsize=12)

# axes[0].plot(median_r, MS_vrot_med, linestyle='-', c='black', linewidth = 1.8, alpha=.85)
# axes[0].plot(median_r, HI_MS_vrot_med, linestyle='--', c='black', linewidth = 1.8, alpha=.85)
# axes[1].plot(median_r, AGy_vrot_med, linestyle='-', c='black', linewidth = 1.8)
# axes[1].plot(median_r, HI_AGy_vrot_med, linestyle='--', c='black', linewidth = 1.8)
# axes[2].plot(median_r, AGo_vrot_med, linestyle='-', c='black', linewidth = 1.8)
# axes[2].plot(median_r, HI_AGo_vrot_med, linestyle='--', c='black', linewidth = 1.8)
# axes[3].plot(median_r, RG_vrot_med, linestyle='-', c='black', linewidth = 1.8)
# axes[3].plot(median_r, HI_RG_vrot_med, linestyle='--', c='black', linewidth = 1.8)

# for ax in axes:
# 	ax.set_xlim(4, 20)
# 	#ax.set_ylabel(r'$\rm v_{\rm rot} \ (km\ s^{-1})$', fontsize=13)
# 	ax.set_ylim(100,300)
# 	ax.tick_params(axis='x',which='both',bottom='on',top='off', direction='out')
# 	ax.tick_params(axis='x',which='both',top='on', direction='in')
# 	ax.tick_params(axis='y',which='both',left='on',top='off', direction='out')
# 	ax.tick_params(axis='y',which='both',right='on', direction='in')
# 	ax.tick_params(which='both', width=2)
# 	ax.tick_params(which='major', length=7)
# 	ax.tick_params(which='minor', length=4)
# 	ax.tick_params(labelsize=12) 
# 	ax.minorticks_on()
# 	for axis in ['top','bottom','left','right']:
# 	        ax.spines[axis].set_linewidth(2)
# axes[3].set_xlabel(r'$\rm Radial\ Distance:\ \it r \ \rm(kpc)$', fontsize=13)
# nbins = len(axes[0].get_yticklabels())-1
# axes[3].yaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='upper'))
# axes[1].yaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='upper'))
# axes[2].yaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='upper'))
# f.subplots_adjust(left=0.17)
# f.text(0.008, 0.5, r'$\rm Rotation\ Velocity:\ \it v_{\rm rot}\ \rm(km\ s^{-1})$', va='center', rotation='vertical', fontsize=13)
# plt.subplots_adjust(wspace=0, hspace=0)
# #plt.savefig('/Users/amandaquirk/Desktop/CO_rotation_curves.pdf', bbox_inches='tight')

# def va(gas, stars):
# 	return gas - stars

# RG_CO_ad = va(CO_vrot, RG_vrot)
# RG_HI_ad = va(HI_vrot, RG_vrot)
# CO_med = round(np.median(RG_CO_ad), 2)
# HI_med = round(np.median(RG_HI_ad), 2)

# plt.hist(RG_CO_ad, histtype='step', bins=range(-200, 200, 15), label='CO med= {}'.format(CO_med), linewidth=1.6,linestyle='--',stacked=True,fill=False, color='darkcyan')
# plt.hist(RG_HI_ad, histtype='step', bins=range(-200, 200, 15),label='HI med= {}'.format(HI_med),linewidth=1.6,stacked=True,fill=False, color='darkgrey')
# plt.legend(loc=2, frameon=False)
# plt.xlabel('Asymmetric Drift (km/s)')
# plt.xlim(-200, 200)
# plt.savefig('/Users/amandaquirk/Desktop/RG_CO_ad.png')

#np.savetxt('/Users/amandaquirk/Desktop/RG_CO_ad.txt', np.c_[RG_r, CO_vrot, RG_CO_ad], fmt='%1.16f', delimiter=' ', header='r (kpc), v (km/s), ad (km/s)')

RG_AD = asymmetric_drift(RG_vrot, RG_CO_v)
MS_AD = asymmetric_drift(MS_vrot, MS_CO_v)
AGy_AD = asymmetric_drift(AGy_vrot, AGy_CO_v)
AGo_AD = asymmetric_drift(AGo_vrot, AGo_CO_v)

