import matplotlib.pyplot as plt 
import numpy as np 
from functions import *

#coordinate data
MS_xi, MS_eta=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/MS_smoothed_chemin.txt', usecols=(0,1), unpack=True)
AGy_xi, AGy_eta=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/AGy_smoothed_chemin.txt', usecols=(0,1), unpack=True)
AGo_xi, AGo_eta=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/AGo_smoothed_chemin.txt', usecols=(0,1), unpack=True)
RG_xi, RG_eta=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/RG_smoothed_chemin.txt', usecols=(0,1), unpack=True)

#import data
MS_r, MS_vrot, MS_HI=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/MS_master_vrot.txt', usecols=(0,2,8), unpack=True)
AGy_r, AGy_vrot, AGy_HI=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/AGy_master_vrot.txt', usecols=(0,2,8), unpack=True)
AGo_r, AGo_vrot, AGo_HI=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/AGo_master_vrot.txt', usecols=(0,2,8), unpack=True)
RG_r, RG_vrot, RG_HI=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/RG_master_vrot.txt', usecols=(0,2,8), unpack=True)

# median_r, MS_vrot_med, HI_MS_vrot_med = np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/MS_med.txt', unpack=True)
# AGy_vrot_med, HI_AGy_vrot_med = np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/AGy_med.txt', usecols=(1,2,), unpack=True)
# AGo_vrot_med, HI_AGo_vrot_med = np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/AGo_med.txt', usecols=(1,2,), unpack=True)
# RG_vrot_med, HI_RG_vrot_med = np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/RG_med.txt', usecols=(1,2,), unpack=True)


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

MS_x=x(np.array(MS_xi), np.array(MS_eta))
MS_y=y(np.array(MS_xi), np.array(MS_eta))

AGy_x=x(np.array(AGy_xi), np.array(AGy_eta))
AGy_y=y(np.array(AGy_xi), np.array(AGy_eta))

AGo_x=x(np.array(AGo_xi), np.array(AGo_eta))
AGo_y=y(np.array(AGo_xi), np.array(AGo_eta))

RG_x=x(np.array(RG_xi), np.array(RG_eta))
RG_y=y(np.array(RG_xi), np.array(RG_eta))

def PA(x,y): #defined as west of north as positioned in maps (inverted x axis)
	if x>0:
		rad=np.arctan(float(y)/x)
		deg=90-(float(rad*180)/np.pi)
	else:
		rad=np.arctan(float(y)/x)
		deg=270-(float(rad*180)/np.pi)
	return deg + 37

MS_PA=[]
for i in range(len(MS_x)):
	MS_PA.append(PA(MS_x[i], MS_y[i]))
AGy_PA=[]
for i in range(len(AGy_x)):
	AGy_PA.append(PA(AGy_x[i], AGy_y[i]))
AGo_PA=[]
for i in range(len(AGo_x)):
	AGo_PA.append(PA(AGo_x[i], AGo_y[i]))
RG_PA=[]
for i in range(len(RG_x)):
	RG_PA.append(PA(RG_x[i], RG_y[i]))

#deprojection factor from rotational velocity equation
def factor(PA_ring, PA_star, i_ring):
	B= [(np.tan(float((a-b)*np.pi) / 180))**2 for a,b in zip(PA_ring, PA_star)]
	C= [(np.cos(float(a*np.pi) / 180))**2 for a in i_ring]
	factor= [np.sqrt(1+(float(b)/c)) for b,c in zip(B,C)]
	return factor

#now we need to assign each star a PA_ring and i_ring using HI data. we define these values by seeing which HI ring a star is cloest to and assigning it the corresponding PA_ring and i_ring value (HI_PA and HI_i respectively). we also will pair a HI velocity to each star to later calculate the asymmetric drift
#reads in the HI data
HI_r, HI_PA, HI_i, HI_v=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/HI_PA_i_vrot.txt', unpack=True)

#below defines a function to determine which HI ring a star is cloest to
def find_nearest_ring(radius):
	idx=[]
	for i in range(len(HI_r)):
		idx.append(np.abs(HI_r[i]-radius))
		n=np.argmin(idx)
	return n #index of the radius closest to the star's radius

MS_PA_ring=[]
MS_i=[]
for j in range(len(MS_r)):
	N=find_nearest_ring(MS_r[j])
	MS_PA_ring.append(HI_PA[N])
	MS_i.append(HI_i[N])

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

RG_PA_ring=[]
RG_i=[]
for j in range(len(RG_r)):
	N=find_nearest_ring(RG_r[j])
	RG_PA_ring.append(HI_PA[N])
	RG_i.append(HI_i[N])

#want the deprojection factor 
MS_factor=factor(MS_PA_ring, MS_PA, MS_i)
AGy_factor=factor(AGy_PA_ring, AGy_PA, AGy_i)
AGo_factor=factor(AGo_PA_ring, AGo_PA, AGo_i)
RG_factor=factor(RG_PA_ring, RG_PA, RG_i)

#median of the deprojection factors
MS_median=np.median(MS_factor)
AGy_median=np.median(AGy_factor)
AGo_median=np.median(AGo_factor)
RG_median=np.median(RG_factor)

#splits data into two groups: ones below and ones above the median of the factor
MS_xi_lowfactor=[a for a, b in zip(MS_xi, MS_factor) if b<MS_median]
AGy_xi_lowfactor=[a for a, b in zip(AGy_xi, AGy_factor) if b<AGy_median]
AGo_xi_lowfactor=[a for a, b in zip(AGo_xi, AGo_factor) if b<AGo_median]
RG_xi_lowfactor=[a for a, b in zip(RG_xi, RG_factor) if b<RG_median]

MS_eta_lowfactor=[a for a, b in zip(MS_eta, MS_factor) if b<MS_median]
AGy_eta_lowfactor=[a for a, b in zip(AGy_eta, AGy_factor) if b<AGy_median]
AGo_eta_lowfactor=[a for a, b in zip(AGo_eta, AGo_factor) if b<AGo_median]
RG_eta_lowfactor=[a for a, b in zip(RG_eta, RG_factor) if b<RG_median]

MS_xi_highfactor=[a for a, b in zip(MS_xi, MS_factor) if b>=MS_median]
AGy_xi_highfactor=[a for a, b in zip(AGy_xi, AGy_factor) if b>=AGy_median]
AGo_xi_highfactor=[a for a, b in zip(AGo_xi, AGo_factor) if b>=AGo_median]
RG_xi_highfactor=[a for a, b in zip(RG_xi, RG_factor) if b>=RG_median]

MS_eta_highfactor=[a for a, b in zip(MS_eta, MS_factor) if b>=MS_median]
AGy_eta_highfactor=[a for a, b in zip(AGy_eta, AGy_factor) if b>=AGy_median]
AGo_eta_highfactor=[a for a, b in zip(AGo_eta, AGo_factor) if b>=AGo_median]
RG_eta_highfactor=[a for a, b in zip(RG_eta, RG_factor) if b>=RG_median]

MS_r_lowfactor=[a for a, b in zip(MS_r, MS_factor) if b<MS_median]
AGy_r_lowfactor=[a for a, b in zip(AGy_r, AGy_factor) if b<AGy_median]
AGo_r_lowfactor=[a for a, b in zip(AGo_r, AGo_factor) if b<AGo_median]
RG_r_lowfactor=[a for a, b in zip(RG_r, RG_factor) if b<RG_median]

MS_vrot_lowfactor=[a for a, b in zip(MS_vrot, MS_factor) if b<MS_median]
AGy_vrot_lowfactor=[a for a, b in zip(AGy_vrot, AGy_factor) if b<AGy_median]
AGo_vrot_lowfactor=[a for a, b in zip(AGo_vrot, AGo_factor) if b<AGo_median]
RG_vrot_lowfactor=[a for a, b in zip(RG_vrot, RG_factor) if b<RG_median]

MS_HI_lowfactor=[a for a, b in zip(MS_HI, MS_factor) if b<MS_median]
AGy_HI_lowfactor=[a for a, b in zip(AGy_HI, AGy_factor) if b<AGy_median]
AGo_HI_lowfactor=[a for a, b in zip(AGo_HI, AGo_factor) if b<AGo_median]
RG_HI_lowfactor=[a for a, b in zip(RG_HI, RG_factor) if b<RG_median]

MS_r_highfactor=[a for a, b in zip(MS_r, MS_factor) if b>=MS_median]
AGy_r_highfactor=[a for a, b in zip(AGy_r, AGy_factor) if b>=AGy_median]
AGo_r_highfactor=[a for a, b in zip(AGo_r, AGo_factor) if b>=AGo_median]
RG_r_highfactor=[a for a, b in zip(RG_r, RG_factor) if b>=RG_median]

MS_vrot_highfactor=[a for a, b in zip(MS_vrot, MS_factor) if b>=MS_median]
AGy_vrot_highfactor=[a for a, b in zip(AGy_vrot, AGy_factor) if b>=AGy_median]
AGo_vrot_highfactor=[a for a, b in zip(AGo_vrot, AGo_factor) if b>=AGo_median]
RG_vrot_highfactor=[a for a, b in zip(RG_vrot, RG_factor) if b>=RG_median]

MS_HI_highfactor=[a for a, b in zip(MS_HI, MS_factor) if b>=MS_median]
AGy_HI_highfactor=[a for a, b in zip(AGy_HI, AGy_factor) if b>=AGy_median]
AGo_HI_highfactor=[a for a, b in zip(AGo_HI, AGo_factor) if b>=AGo_median]
RG_HI_highfactor=[a for a, b in zip(RG_HI, RG_factor) if b>=RG_median]

#rotation curve plot
def rotation_curve(r, star_vrot, HI_vrot, age, description):
	if age=='MS':
		color='b'
	if age=='RG':
		color='r'
	if age!='MS' and age !='RG':
		color='m'
	plt.scatter(r, star_vrot, s=1, c='{}'.format(color)) #should put quotes around the bracket?
	plt.scatter(r, HI_vrot, s=1, c='k')
	plt.xlim(4,20)
	plt.ylim(100,300)
	plt.xlabel('r [kpc]')
	plt.ylabel('rotational v [km/s]')
	plt.savefig('/Users/amandaquirk/Desktop/{}_vrot_{}.png'.format(age, description))
	plt.close()

# rotation_curve(MS_r_lowfactor, MS_vrot_lowfactor, MS_HI_lowfactor, 'MS', 'lowfactor')
# rotation_curve(AGy_r_lowfactor, AGy_vrot_lowfactor, AGy_HI_lowfactor, 'AGy', 'lowfactor')
# rotation_curve(AGo_r_lowfactor, AGo_vrot_lowfactor, AGo_HI_lowfactor, 'AGo', 'lowfactor')
# rotation_curve(RG_r_lowfactor, RG_vrot_lowfactor, RG_HI_lowfactor, 'RG', 'lowfactor')

# rotation_curve(MS_r_highfactor, MS_vrot_highfactor, MS_HI_highfactor, 'MS', 'highfactor')
# rotation_curve(AGy_r_highfactor, AGy_vrot_highfactor, AGy_HI_highfactor, 'AGy', 'highfactor')
# rotation_curve(AGo_r_highfactor, AGo_vrot_highfactor, AGo_HI_highfactor, 'AGo', 'highfactor')
# rotation_curve(RG_r_highfactor, RG_vrot_highfactor, RG_HI_highfactor, 'RG', 'highfactor')
# fig, ax=plt.subplots(1)
# for axis in ['top','bottom','left','right']:
#         ax.spines[axis].set_linewidth(2)
# ax.tick_params(axis='x',which='minor',bottom='on')
# plt.tick_params(which='both', width=2)
# plt.tick_params(which='major', length=7)
# plt.tick_params(which='minor', length=4)
# plt.tick_params(labelsize=10) 
# plt.minorticks_on()
# #plt.figure(figsize=(4, 6))
# plt.scatter(MS_xi_lowfactor, MS_eta_lowfactor, s=1, c='b')
# plt.scatter(AGy_xi_lowfactor, AGy_eta_lowfactor, s=1, c='b')
# plt.scatter(AGo_xi_lowfactor, AGo_eta_lowfactor, s=1, c='b')
# plt.scatter(RG_xi_lowfactor, RG_eta_lowfactor, s=1, c='b')
# plt.scatter(MS_xi_highfactor, MS_eta_highfactor, s=1, c='r')
# plt.scatter(AGy_xi_highfactor, AGy_eta_highfactor, s=1, c='r')
# plt.scatter(AGo_xi_highfactor, AGo_eta_highfactor, s=1, c='r')
# plt.scatter(RG_xi_highfactor, RG_eta_highfactor, s=1, c='r')
# plt.scatter(0,0, s=2, c='k')
# plt.gca().invert_xaxis()
# plt.xlim(14,-8)
# plt.ylim(-12,16)
# plt.xlabel('xi (kpc)')
# plt.ylabel('eta (kpc)')
# plt.show()

from matplotlib import rc 
from matplotlib.ticker import MaxNLocator

RG_r_med, RG_vrot_med_high = median_line(RG_r_highfactor, RG_vrot_highfactor)
MS_r_med, MS_vrot_med_high = median_line(MS_r_highfactor, MS_vrot_highfactor)
AGy_r_med, AGy_vrot_med_high = median_line(AGy_r_highfactor, AGy_vrot_highfactor)
AGo_r_med, AGo_vrot_med_high = median_line(AGo_r_highfactor, AGo_vrot_highfactor)
HI_RG_r_med, HI_RG_vrot_med_high = median_line(RG_r_highfactor, RG_HI_highfactor)
HI_MS_r_med, HI_MS_vrot_med_high = median_line(MS_r_highfactor, MS_HI_highfactor)
HI_AGy_r_med, HI_AGy_vrot_med_high = median_line(AGy_r_highfactor, AGy_HI_highfactor)
HI_AGo_r_med, HI_AGo_vrot_med_high = median_line(AGo_r_highfactor, AGo_HI_highfactor)

#compute median AD 
RG_AD = asymmetric_drift(RG_vrot_highfactor, RG_HI_highfactor)
AGy_AD = asymmetric_drift(AGy_vrot_highfactor, AGy_HI_highfactor)
AGo_AD = asymmetric_drift(AGo_vrot_highfactor, AGo_HI_highfactor)
MS_AD = asymmetric_drift(MS_vrot_highfactor, MS_HI_highfactor)

rc('font', family = 'serif')
f, axes= plt.subplots(4,1, sharey=False, sharex=True, figsize=(4,9.8))
axes[0].scatter(MS_r_highfactor, MS_HI_highfactor, s=2, c='darkgray')
axes[0].scatter(MS_r_highfactor, MS_vrot_highfactor, s=2, c='b', alpha=0.4)
axes[1].scatter(AGy_r_highfactor, AGy_HI_highfactor, s=2, c='darkgray')
axes[1].scatter(AGy_r_highfactor, AGy_vrot_highfactor, s=2, c='m', alpha=0.4)
axes[2].scatter(AGo_r_highfactor, AGo_HI_highfactor, s=2, c='darkgray')
axes[2].scatter(AGo_r_highfactor, AGo_vrot_highfactor, s=2, c='green', alpha=0.4)
axes[3].scatter(RG_r_highfactor, RG_HI_highfactor, s=2, c='darkgray')
axes[3].scatter(RG_r_highfactor, RG_vrot_highfactor, s=2, c='r', alpha=0.4)
axes[0].annotate('MS', xy=(19,128), horizontalalignment='right', fontsize=10)
axes[0].annotate(r'$\overline{v}_{\rm a}=$'+ '${}$'.format(round(np.median(MS_AD),2)) + r'$\rm km\ s^{-1}$', xy=(19,112), horizontalalignment='right', fontsize=10)
axes[0].plot(MS_r_med, MS_vrot_med_high, linestyle='-', c='black', linewidth = 1.8, alpha=.85)
axes[0].plot(HI_MS_r_med, HI_MS_vrot_med_high, linestyle='--', c='black', linewidth = 1.8, alpha=.85)
axes[1].plot(AGy_r_med, AGy_vrot_med_high, linestyle='-', c='black', linewidth = 1.8)
axes[1].plot(HI_AGy_r_med, HI_AGy_vrot_med_high, linestyle='--', c='black', linewidth = 1.8)
axes[2].plot(AGo_r_med, AGo_vrot_med_high, linestyle='-', c='black', linewidth = 1.8)
axes[2].plot(HI_AGo_r_med, HI_AGo_vrot_med_high, linestyle='--', c='black', linewidth = 1.8)
axes[3].plot(RG_r_med, RG_vrot_med_high, linestyle='-', c='black', linewidth = 1.8)
axes[3].plot(HI_RG_r_med, HI_RG_vrot_med_high, linestyle='--', c='black', linewidth = 1.8)
axes[1].annotate('young AGB', xy=(19,128), horizontalalignment='right', fontsize=10)
axes[1].annotate(r'$\overline{v}_{\rm a}=$'+ '${}$'.format(round(np.median(AGy_AD),2)) + r'$\rm km\ s^{-1}$', xy=(19,112), horizontalalignment='right', fontsize=10)
axes[2].annotate('older AGB', xy=(19,128), horizontalalignment='right', fontsize=10)
axes[2].annotate(r'$\overline{v}_{\rm a}=$'+ '${}$'.format(round(np.median(AGo_AD),2)) + r'$\rm km\ s^{-1}$', xy=(19,112), horizontalalignment='right', fontsize=10)
axes[3].annotate('RGB', xy=(19,128), horizontalalignment='right', fontsize=10)
axes[3].annotate(r'$\overline{v}_{\rm a}=$'+ '${}$'.format(round(np.median(RG_AD),2)) + r'$\rm km\ s^{-1}$', xy=(19,112), horizontalalignment='right', fontsize=10)

for ax in axes:
	ax.set_xlim(4, 20)
	#ax.set_ylabel(r'$\rm v_{\rm rot} \ (km\ s^{-1})$', fontsize=13)
	ax.set_ylim(100,300)
	ax.tick_params(axis='x',which='both',bottom='on',top='off', direction='out')
	ax.tick_params(axis='x',which='both',top='on', direction='in')
	ax.tick_params(axis='y',which='both',left='on',top='off', direction='out')
	ax.tick_params(axis='y',which='both',right='on', direction='in')
	ax.tick_params(which='both', width=2)
	ax.tick_params(which='major', length=7)
	ax.tick_params(which='minor', length=4)
	ax.tick_params(labelsize=12) 
	ax.minorticks_on()
	for axis in ['top','bottom','left','right']:
	        ax.spines[axis].set_linewidth(2)
axes[3].set_xlabel(r'$\rm Radial\ Distance:\ \it r \ \rm(kpc)$', fontsize=13)
nbins = len(axes[0].get_yticklabels())-1
axes[3].yaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='upper'))
axes[1].yaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='upper'))
axes[2].yaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='upper'))
f.subplots_adjust(left=0.17)
f.text(0.008, 0.5, r'$\rm Rotation\ Velocity:\ \itv_{\rm rot}\ \rm(km\ s^{-1})$', va='center', rotation='vertical', fontsize=13)
plt.subplots_adjust(wspace=0, hspace=0)
plt.savefig('/Users/amandaquirk/Desktop/high_rotation_curves.pdf', bbox_inches='tight')
