import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib import rc
from matplotlib.ticker import MaxNLocator
from functions import *

#import data
#importing all of the data
#xi (kpc), eta (kpc), average v(km/s), v err,var, n, HI main, HI close <-- header of data file to be read in 
MS_xi, MS_eta, MS_v, MS_dispersion, MS_n, MS_HImain, MS_HIclose=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/MS_smoothed_chemin.txt', usecols=(0,1,2,4,5,6,7,), unpack=True)
AGy_xi, AGy_eta, AGy_v, AGy_dispersion,AGy_n, AGy_HImain, AGy_HIclose=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/AGy_smoothed_chemin.txt', usecols=(0,1,2,4,5,6,7,), unpack=True)
AGo_xi, AGo_eta, AGo_v, AGo_dispersion, AGo_n, AGo_HImain, AGo_HIclose=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/AGo_smoothed_chemin.txt', usecols=(0,1,2,4,5,6,7,), unpack=True)
RG_xi, RG_eta, RG_v, RG_dispersion, RG_n, RG_HImain, RG_HIclose=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/RG_smoothed_chemin.txt', usecols=(0,1,2,4,5,6,7,), unpack=True)

#import rotation curve data
#r (kpc), v_rot no model (km/s), v_rot tr model, smoothed PA, sPA v_rot no model, sPA v_rot tr model, n, HI main no model, HI main tr, HI close no model, HI close tr
MS_r, MS_vrot_tr, MS_HImain_vrot_tr, MS_HIclose_vrot_tr=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/MS_master_vrot.txt', usecols=(0,2,8,10), unpack=True)
AGy_r, AGy_vrot_tr, AGy_HImain_vrot_tr, AGy_HIclose_vrot_tr=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/AGy_master_vrot.txt', usecols=(0,2,8,10), unpack=True)
AGo_r, AGo_vrot_tr, AGo_HImain_vrot_tr, AGo_HIclose_vrot_tr=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/AGo_master_vrot.txt', usecols=(0,2,8,10), unpack=True)
RG_r, RG_vrot_tr, RG_HImain_vrot_tr, RG_HIclose_vrot_tr=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/RG_master_vrot.txt', usecols=(0,2,8,10), unpack=True)

MS_Av=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/MS_dust.txt', usecols=(2,), unpack=True)
AGy_Av=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/AGy_dust.txt', usecols=(2,), unpack=True)
AGo_Av=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/AGo_dust.txt', usecols=(2,), unpack=True)
RG_Av=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/RG_dust.txt', usecols=(2,), unpack=True)

def dispersion_n_plt(n, dispersion, age):
	plt.scatter(n, dispersion, alpha=.3)
	plt.ylim(10,150)
	plt.ylabel('Dispersion (km/s)')
	plt.xlabel('Number of HI Peaks')
	plt.savefig('/Users/amandaquirk/Desktop/{}_var_n.png'.format(age))
	plt.close()

# dispersion_n_plt(MS_n, MS_dispersion, 'MS')
# dispersion_n_plt(AGy_n, AGy_dispersion, 'AGy')
# dispersion_n_plt(AGo_n, AGo_dispersion, 'AGo')
# dispersion_n_plt(RG_n, RG_dispersion, 'RG')

def HI_diff(HImain, HIclose):
	return (HImain-HIclose)

# MS_HI_diff=HI_diff(MS_HImain,MS_HIclose)
# AGy_HI_diff=HI_diff(AGy_HImain,AGy_HIclose)
# AGo_HI_diff=HI_diff(AGo_HImain,AGo_HIclose)
# RG_HI_diff=HI_diff(RG_HImain,RG_HIclose)

MS_HI_diff_d=HI_diff(MS_HImain_vrot_tr,MS_HIclose_vrot_tr)
AGy_HI_diff_d=HI_diff(AGy_HImain_vrot_tr,AGy_HIclose_vrot_tr)
AGo_HI_diff_d=HI_diff(AGo_HImain_vrot_tr,AGo_HIclose_vrot_tr)
RG_HI_diff_d=HI_diff(RG_HImain_vrot_tr,RG_HIclose_vrot_tr)

def dispersion_HI_diff(diff, dispersion, age, deprojected=True):
	if deprojected==True:
		label='deprojected_vrot'
		plt.xlim(-500,500)
	else:
		label='projected_LOS'
	plt.scatter(diff, dispersion, alpha=.3)
	plt.ylim(10,150)
	plt.ylabel('Dispersion (km/s)')
	plt.xlabel('{} HI_main-HI_close (km/s)'.format(label))
	plt.savefig('/Users/amandaquirk/Desktop/{}_var_HIdiff_{}.png'.format(age,label))
	plt.close()

# dispersion_HI_diff(MS_HI_diff, MS_dispersion, 'MS', deprojected=False)
# dispersion_HI_diff( AGy_HI_diff,AGy_dispersion, 'AGy', deprojected=False)
# dispersion_HI_diff( AGo_HI_diff,AGo_dispersion, 'AGo', deprojected=False)
# dispersion_HI_diff(RG_HI_diff, RG_dispersion, 'RG', deprojected=False)

# dispersion_HI_diff(MS_HI_diff_d, MS_dispersion, 'MS')
# dispersion_HI_diff(AGy_HI_diff_d, AGy_dispersion, 'AGy')
# dispersion_HI_diff(AGo_HI_diff_d, AGo_dispersion, 'AGo')
# dispersion_HI_diff(RG_HI_diff_d, RG_dispersion, 'RG')


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

def dust_rotation_curve(Av, r, vrot, HI_vrot, age, dusty=True):
	if age=='MS':
		color='b'
	if age=='RG':
		color='r'
	if age!='MS' and age !='RG':
		color='m'
	r_dust=[]
	r_clear=[]
	vrot_dust=[]
	vrot_clear=[]
	HI_vrot_clear=[]
	HI_vrot_dust=[]
	for i in range(len(r)):
		if Av[i]<1:
			r_clear.append(r[i])
			vrot_clear.append(vrot[i])
			HI_vrot_clear.append(HI_vrot[i])
		else:
			r_dust.append(r[i])
			vrot_dust.append(vrot[i])
			HI_vrot_dust.append(HI_vrot[i])
	if dusty==True:
		description='dusty'
		x=r_dust
		star_y=vrot_dust
		HI_y=HI_vrot_dust
	if dusty==False:
		description='translucent'
		x=r_clear
		star_y=vrot_clear
		HI_y=HI_vrot_clear 

	plt.scatter(x, star_y, s=1, c='{}'.format(color))
	plt.scatter(x, HI_y, s=1, c='k')
	plt.xlim(4,20)
	plt.ylim(100,300)
	plt.xlabel('r [kpc]')
	plt.ylabel('rotational v [km/s]')
	plt.savefig('/Users/amandaquirk/Desktop/{}_rotationcurve_{}.png'.format(description,age))
	plt.close()

MS_r_dust=[a for a,b in zip(MS_r, MS_Av) if b>=1]
AGy_r_dust=[a for a,b in zip(AGy_r, AGy_Av) if b>=1]
AGo_r_dust=[a for a,b in zip(AGo_r, AGo_Av) if b>=1]
RG_r_dust=[a for a,b in zip(RG_r, RG_Av) if b>=1]
MS_r_clear=[a for a,b in zip(MS_r, MS_Av) if b<1]
AGy_r_clear=[a for a,b in zip(AGy_r, AGy_Av) if b<1]
AGo_r_clear=[a for a,b in zip(AGo_r, AGo_Av) if b<1]
RG_r_clear=[a for a,b in zip(RG_r, RG_Av) if b<1]

MS_vrot_dust=[a for a,b in zip(MS_vrot_tr, MS_Av) if b>=1]
AGy_vrot_dust=[a for a,b in zip(AGy_vrot_tr, AGy_Av) if b>=1]
AGo_vrot_dust=[a for a,b in zip(AGo_vrot_tr, AGo_Av) if b>=1]
RG_vrot_dust=[a for a,b in zip(RG_vrot_tr, RG_Av) if b>=1]
MS_vrot_clear=[a for a,b in zip(MS_vrot_tr, MS_Av) if b<1]
AGy_vrot_clear=[a for a,b in zip(AGy_vrot_tr, AGy_Av) if b<1]
AGo_vrot_clear=[a for a,b in zip(AGo_vrot_tr, AGo_Av) if b<1]
RG_vrot_clear=[a for a,b in zip(RG_vrot_tr, RG_Av) if b<1]

MS_HIvrot_dust=[a for a,b in zip(MS_HImain_vrot_tr, MS_Av) if b>=1]
AGy_HIvrot_dust=[a for a,b in zip(AGy_HImain_vrot_tr, AGy_Av) if b>=1]
AGo_HIvrot_dust=[a for a,b in zip(AGo_HImain_vrot_tr, AGo_Av) if b>=1]
RG_HIvrot_dust=[a for a,b in zip(RG_HImain_vrot_tr, RG_Av) if b>=1]
MS_HIvrot_clear=[a for a,b in zip(MS_HImain_vrot_tr, MS_Av) if b<1]
AGy_HIvrot_clear=[a for a,b in zip(AGy_HImain_vrot_tr, AGy_Av) if b<1]
AGo_HIvrot_clear=[a for a,b in zip(AGo_HImain_vrot_tr, AGo_Av) if b<1]
RG_HIvrot_clear=[a for a,b in zip(RG_HImain_vrot_tr, RG_Av) if b<1]

AGo_r_med, AGo_vrot_med = median_line(AGo_r_clear, AGo_vrot_clear)
MS_r_med, MS_vrot_med = median_line(MS_r_clear, MS_vrot_clear)
AGy_r_med, AGy_vrot_med = median_line(AGy_r_clear, AGy_vrot_clear)
RG_r_med, RG_vrot_med = median_line(RG_r_clear, RG_vrot_clear)
AGo_r_med, HI_AGo_vrot_med = median_line(AGo_r_clear, AGo_HIvrot_clear)
MS_r_med, HI_MS_vrot_med = median_line(MS_r_clear, MS_HIvrot_clear)
AGy_r_med, HI_AGy_vrot_med = median_line(AGy_r_clear, AGy_HIvrot_clear)
RG_r_med, HI_RG_vrot_med = median_line(RG_r_clear, RG_HIvrot_clear)

AGo_AD = asymmetric_drift(AGo_vrot_clear, AGo_HIvrot_clear)
MS_AD = asymmetric_drift(MS_vrot_clear, MS_HIvrot_clear)
AGy_AD = asymmetric_drift(AGy_vrot_clear, AGy_HIvrot_clear)
RG_AD = asymmetric_drift(RG_vrot_clear, RG_HIvrot_clear)

rc('font', family = 'serif')
f, axes= plt.subplots(4,1, sharey=False, sharex=True, figsize=(4,9.8))
axes[0].scatter(MS_r_clear, MS_HIvrot_clear, s=2, c='darkgray')
axes[0].scatter(MS_r_clear, MS_vrot_clear, s=2, c='b', alpha=0.4)
axes[1].scatter(AGy_r_clear, AGy_HIvrot_clear, s=2, c='darkgray')
axes[1].scatter(AGy_r_clear, AGy_vrot_clear, s=2, c='m', alpha=0.4)
axes[2].scatter(AGo_r_clear, AGo_HIvrot_clear, s=2, c='darkgray')
axes[2].scatter(AGo_r_clear, AGo_vrot_clear, s=2, c='green', alpha=0.4)
axes[3].scatter(RG_r_clear, RG_HIvrot_clear, s=2, c='darkgray')
axes[3].scatter(RG_r_clear, RG_vrot_clear, s=2, c='r', alpha=0.4)
axes[0].annotate('MS', xy=(19,128), horizontalalignment='right', fontsize=10)
axes[0].annotate(r'$\overline{v}_{\rm a}=$'+ '${}$'.format(round(np.median(MS_AD),2)) + r'$\rm km\ s^{-1}$', xy=(19,112), horizontalalignment='right', fontsize=10)
axes[0].plot(MS_r_med, MS_vrot_med, linestyle='-', c='black', linewidth = 1.8, alpha=.85)
axes[0].plot(MS_r_med, HI_MS_vrot_med, linestyle='--', c='black', linewidth = 1.8, alpha=.85)
axes[1].plot(AGy_r_med, AGy_vrot_med, linestyle='-', c='black', linewidth = 1.8)
axes[1].plot(AGy_r_med, HI_AGy_vrot_med, linestyle='--', c='black', linewidth = 1.8)
axes[2].plot(AGo_r_med, AGo_vrot_med, linestyle='-', c='black', linewidth = 1.8)
axes[2].plot(AGo_r_med, HI_AGo_vrot_med, linestyle='--', c='black', linewidth = 1.8)
axes[3].plot(RG_r_med, RG_vrot_med, linestyle='-', c='black', linewidth = 1.8)
axes[3].plot(RG_r_med, HI_RG_vrot_med, linestyle='--', c='black', linewidth = 1.8)
axes[1].annotate('young AGB', xy=(19,128), horizontalalignment='right', fontsize=10)
axes[1].annotate(r'$\overline{v}_{\rm a}=$'+ '${}$'.format(round(np.median(AGy_AD),2)) + r'$\rm km\ s^{-1}$', xy=(19,112), horizontalalignment='right', fontsize=10)
axes[2].annotate('older AGB', xy=(19,128), horizontalalignment='right', fontsize=10)
axes[2].annotate(r'$\overline{v}_{\rm a}=$'+ '${}$'.format(round(np.median(AGo_AD),2)) + r'$\rm km\ s^{-1}$', xy=(19,112), horizontalalignment='right', fontsize=10)
axes[3].annotate('RGB', xy=(19,128), horizontalalignment='right', fontsize=10)
axes[3].annotate(r'$\overline{v}_{\rm a}=$'+ '${}$'.format(round(np.median(RG_AD),2)) + r'$\rm km\ s^{-1}$', xy=(19,112), horizontalalignment='right', fontsize=10)

axes[0].annotate(r'$\rm f=0.57$', xy=(5,112), horizontalalignment='left', fontsize=10)
axes[1].annotate(r'$\rm f=0.78$', xy=(5,112), horizontalalignment='left', fontsize=10)
axes[2].annotate(r'$\rm f=0.74$', xy=(5,112), horizontalalignment='left', fontsize=10)
axes[3].annotate(r'$\rm f=0.77$', xy=(5,112), horizontalalignment='left', fontsize=10)


for ax in axes:
	ax.set_xlim(4, 20)
	#ax.set_ylabel(r'$\rm v_{\rm rot} \ (km^{-1})$', fontsize=13)
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
axes[3].set_xlabel(r'$\rm Radial\ Distance:\ \itr \ \rm(kpc)$', fontsize=13)
nbins = len(axes[0].get_yticklabels())-1
axes[3].yaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='upper'))
axes[1].yaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='upper'))
axes[2].yaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='upper'))
f.subplots_adjust(left=0.17)
f.text(0.008, 0.5, r'$\rm Rotation\ Velocity:\ \itv_{\rm rot}\ \rm(km\ s^{-1})$', va='center', rotation='vertical', fontsize=13)
plt.subplots_adjust(wspace=0, hspace=0)
plt.savefig('/Users/amandaquirk/Desktop/clear_rcs.pdf', bbox_inches='tight')

# dust_rotation_curve(MS_Av, MS_r, MS_vrot_tr, MS_HImain_vrot_tr, 'MS', dusty=True)
# dust_rotation_curve(MS_Av, MS_r, MS_vrot_tr, MS_HImain_vrot_tr, 'MS', dusty=False)
# dust_rotation_curve(AGy_Av, AGy_r, AGy_vrot_tr, AGy_HImain_vrot_tr, 'AGy', dusty=True)
# dust_rotation_curve(AGy_Av, AGy_r, AGy_vrot_tr, AGy_HImain_vrot_tr, 'AGy', dusty=False)
# dust_rotation_curve(AGo_Av, AGo_r, AGo_vrot_tr, AGo_HImain_vrot_tr, 'AGo', dusty=True)
# dust_rotation_curve(AGo_Av, AGo_r, AGo_vrot_tr, AGo_HImain_vrot_tr, 'AGo', dusty=False)
# dust_rotation_curve(RG_Av, RG_r, RG_vrot_tr, RG_HImain_vrot_tr, 'RG', dusty=True)
# dust_rotation_curve(RG_Av, RG_r, RG_vrot_tr, RG_HImain_vrot_tr, 'RG', dusty=False)

MS_r_0=[] #points where HImain=HIclose
MS_vrot_0=[]
MS_HImain_0=[]
MS_r_not0=[] #points where HImain=/=HIclose
MS_vrot_not0=[]
MS_HImain_not0=[]
for i in range(len(MS_HIclose_vrot_tr)):
	if MS_HI_diff_d[i]==0:
		MS_r_0.append(MS_r[i])
		MS_vrot_0.append(MS_vrot_tr[i])
		MS_HImain_0.append(MS_HImain_vrot_tr[i])
	else:
		MS_r_not0.append(MS_r[i])
		MS_vrot_not0.append(MS_vrot_tr[i])
		MS_HImain_not0.append(MS_HImain_vrot_tr[i])

AGy_r_0=[] #points where HImain=HIclose
AGy_vrot_0=[]
AGy_HImain_0=[]
AGy_r_not0=[] #points where HImain=/=HIclose
AGy_vrot_not0=[]
AGy_HImain_not0=[]
for i in range(len(AGy_HIclose_vrot_tr)):
	if AGy_HI_diff_d[i]==0:
		AGy_r_0.append(AGy_r[i])
		AGy_vrot_0.append(AGy_vrot_tr[i])
		AGy_HImain_0.append(AGy_HImain_vrot_tr[i])
	else:
		AGy_r_not0.append(AGy_r[i])
		AGy_vrot_not0.append(AGy_vrot_tr[i])
		AGy_HImain_not0.append(AGy_HImain_vrot_tr[i])

AGo_r_0=[] #points where HImain=HIclose
AGo_vrot_0=[]
AGo_HImain_0=[]
AGo_r_not0=[] #points where HImain=/=HIclose
AGo_vrot_not0=[]
AGo_HImain_not0=[]
for i in range(len(AGo_HIclose_vrot_tr)):
	if AGo_HI_diff_d[i]==0:
		AGo_r_0.append(AGo_r[i])
		AGo_vrot_0.append(AGo_vrot_tr[i])
		AGo_HImain_0.append(AGo_HImain_vrot_tr[i])
	else:
		AGo_r_not0.append(AGo_r[i])
		AGo_vrot_not0.append(AGo_vrot_tr[i])
		AGo_HImain_not0.append(AGo_HImain_vrot_tr[i])

RG_r_0=[] #points where HImain=HIclose
RG_vrot_0=[]
RG_HImain_0=[]
RG_r_not0=[] #points where HImain=/=HIclose
RG_vrot_not0=[]
RG_HImain_not0=[]
for i in range(len(RG_HIclose_vrot_tr)):
	if RG_HI_diff_d[i]==0:
		RG_r_0.append(RG_r[i])
		RG_vrot_0.append(RG_vrot_tr[i])
		RG_HImain_0.append(RG_HImain_vrot_tr[i])
	else:
		RG_r_not0.append(RG_r[i])
		RG_vrot_not0.append(RG_vrot_tr[i])
		RG_HImain_not0.append(RG_HImain_vrot_tr[i])

# AGo_r_med, AGo_vrot_med = median_line(AGo_r_0, AGo_vrot_0)
# MS_r_med, MS_vrot_med = median_line(MS_r_0, MS_vrot_0)
# AGy_r_med, AGy_vrot_med = median_line(AGy_r_0, AGy_vrot_0)
# RG_r_med, RG_vrot_med = median_line(RG_r_0, RG_vrot_0)
# AGo_r_med, HI_AGo_vrot_med = median_line(AGo_r_0, AGo_HImain_0)
# MS_r_med, HI_MS_vrot_med = median_line(MS_r_0, MS_HImain_0)
# AGy_r_med, HI_AGy_vrot_med = median_line(AGy_r_0, AGy_HImain_0)
# RG_r_med, HI_RG_vrot_med = median_line(RG_r_0, RG_HImain_0)

# AGo_AD = asymmetric_drift(AGo_vrot_0, AGo_HImain_0)
# MS_AD = asymmetric_drift(MS_vrot_0, MS_HImain_0)
# AGy_AD = asymmetric_drift(AGy_vrot_0, AGy_HImain_0)
# RG_AD = asymmetric_drift(RG_vrot_0, RG_HImain_0)

# rc('font', family = 'serif')
# f, axes= plt.subplots(4,1, sharey=False, sharex=True, figsize=(4,9.8))
# axes[0].scatter(MS_r_0, MS_HImain_0, s=2, c='darkgray')
# axes[0].scatter(MS_r_0, MS_vrot_0, s=2, c='b', alpha=0.4)
# axes[1].scatter(AGy_r_0, AGy_HImain_0, s=2, c='darkgray')
# axes[1].scatter(AGy_r_0, AGy_vrot_0, s=2, c='m', alpha=0.4)
# axes[2].scatter(AGo_r_0, AGo_HImain_0, s=2, c='darkgray')
# axes[2].scatter(AGo_r_0, AGo_vrot_0, s=2, c='green', alpha=0.4)
# axes[3].scatter(RG_r_0, RG_HImain_0, s=2, c='darkgray')
# axes[3].scatter(RG_r_0, RG_vrot_0, s=2, c='r', alpha=0.4)
# axes[0].annotate('MS', xy=(19,128), horizontalalignment='right', fontsize=10)
# axes[0].annotate(r'$\overline{v}_{\rm a}=$'+ '${}$'.format(round(np.median(MS_AD),2)) + r'$\rm km\ s^{-1}$', xy=(19,112), horizontalalignment='right', fontsize=10)
# axes[0].plot(MS_r_med, MS_vrot_med, linestyle='-', c='black', linewidth = 1.8, alpha=.85)
# axes[0].plot(MS_r_med, HI_MS_vrot_med, linestyle='--', c='black', linewidth = 1.8, alpha=.85)
# axes[1].plot(AGy_r_med, AGy_vrot_med, linestyle='-', c='black', linewidth = 1.8)
# axes[1].plot(AGy_r_med, HI_AGy_vrot_med, linestyle='--', c='black', linewidth = 1.8)
# axes[2].plot(AGo_r_med, AGo_vrot_med, linestyle='-', c='black', linewidth = 1.8)
# axes[2].plot(AGo_r_med, HI_AGo_vrot_med, linestyle='--', c='black', linewidth = 1.8)
# axes[3].plot(RG_r_med, RG_vrot_med, linestyle='-', c='black', linewidth = 1.8)
# axes[3].plot(RG_r_med, HI_RG_vrot_med, linestyle='--', c='black', linewidth = 1.8)
# axes[1].annotate('young AGB', xy=(19,128), horizontalalignment='right', fontsize=10)
# axes[1].annotate(r'$\overline{v}_{\rm a}=$'+ '${}$'.format(round(np.median(AGy_AD),2)) + r'$\rm km\ s^{-1}$', xy=(19,112), horizontalalignment='right', fontsize=10)
# axes[2].annotate('older AGB', xy=(19,128), horizontalalignment='right', fontsize=10)
# axes[2].annotate(r'$\overline{v}_{\rm a}=$'+ '${}$'.format(round(np.median(AGo_AD),2)) + r'$\rm km\ s^{-1}$', xy=(19,112), horizontalalignment='right', fontsize=10)
# axes[3].annotate('RGB', xy=(19,128), horizontalalignment='right', fontsize=10)
# axes[3].annotate(r'$\overline{v}_{\rm a}=$'+ '${}$'.format(round(np.median(RG_AD),2)) + r'$\rm km\ s^{-1}$', xy=(19,112), horizontalalignment='right', fontsize=10)

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
# f.text(0.008, 0.5, r'$\rm Rotation\ Velocity:\ \itv_{\rm rot}\ \rm(km\ s^{-1})$', va='center', rotation='vertical', fontsize=13)
# plt.subplots_adjust(wspace=0, hspace=0)
# plt.savefig('/Users/amandaquirk/Desktop/HIsame_curves.pdf', bbox_inches='tight')

# rotation_curve(MS_r_0, MS_vrot_0, MS_HImain_0, 'MS', 0)
# rotation_curve(MS_r_not0, MS_vrot_not0, MS_HImain_not0, 'MS', 'not0')
# rotation_curve(AGy_r_0, AGy_vrot_0, AGy_HImain_0, 'AGy', 0)
# rotation_curve(AGy_r_not0, AGy_vrot_not0, AGy_HImain_not0, 'AGy', 'not0')
# rotation_curve(AGo_r_0, AGo_vrot_0, AGo_HImain_0, 'AGo', 0)
# rotation_curve(AGo_r_not0, AGo_vrot_not0, AGo_HImain_not0, 'AGo', 'not0')
# rotation_curve(RG_r_0, RG_vrot_0, RG_HImain_0, 'RG', 0)
# rotation_curve(RG_r_not0, RG_vrot_not0, RG_HImain_not0, 'RG', 'not0')

#doing the rotation curves for different n
# MS_r_1=[] #points where n=1
# MS_vrot_1=[]
# MS_HImain_1=[]
# MS_r_2=[] #points where n=2
# MS_vrot_2=[]
# MS_HImain_2=[]
# MS_r_3=[] #points where n=3
# MS_vrot_3=[]
# MS_HImain_3=[]
# MS_r_4=[] #points where n=4
# MS_vrot_4=[]
# MS_HImain_4=[]
# MS_r_5=[] #points where n=5
# MS_vrot_5=[]
# MS_HImain_5=[]
# for i in range(len(MS_r)):
# 	if MS_n[i]==1:
# 		MS_r_1.append(MS_r[i])
# 		MS_vrot_1.append(MS_vrot_tr[i])
# 		MS_HImain_1.append(MS_HImain_vrot_tr[i])
# 	if MS_n[i]==2:
# 		MS_r_2.append(MS_r[i])
# 		MS_vrot_2.append(MS_vrot_tr[i])
# 		MS_HImain_2.append(MS_HImain_vrot_tr[i])
# 	if MS_n[i]==3:
# 		MS_r_3.append(MS_r[i])
# 		MS_vrot_3.append(MS_vrot_tr[i])
# 		MS_HImain_3.append(MS_HImain_vrot_tr[i])
# 	if MS_n[i]==4:
# 		MS_r_4.append(MS_r[i])
# 		MS_vrot_4.append(MS_vrot_tr[i])
# 		MS_HImain_4.append(MS_HImain_vrot_tr[i])
# 	if MS_n[i]==5:
# 		MS_r_5.append(MS_r[i])
# 		MS_vrot_5.append(MS_vrot_tr[i])
# 		MS_HImain_5.append(MS_HImain_vrot_tr[i])

# AGy_r_1=[] #points where n=1
# AGy_vrot_1=[]
# AGy_HImain_1=[]
# AGy_r_2=[] #points where n=2
# AGy_vrot_2=[]
# AGy_HImain_2=[]
# AGy_r_3=[] #points where n=3
# AGy_vrot_3=[]
# AGy_HImain_3=[]
# AGy_r_4=[] #points where n=4
# AGy_vrot_4=[]
# AGy_HImain_4=[]
# AGy_r_5=[] #points where n=5
# AGy_vrot_5=[]
# AGy_HImain_5=[]
# for i in range(len(AGy_r)):
# 	if AGy_n[i]==1:
# 		AGy_r_1.append(AGy_r[i])
# 		AGy_vrot_1.append(AGy_vrot_tr[i])
# 		AGy_HImain_1.append(AGy_HImain_vrot_tr[i])
# 	if AGy_n[i]==2:
# 		AGy_r_2.append(AGy_r[i])
# 		AGy_vrot_2.append(AGy_vrot_tr[i])
# 		AGy_HImain_2.append(AGy_HImain_vrot_tr[i])
# 	if AGy_n[i]==3:
# 		AGy_r_3.append(AGy_r[i])
# 		AGy_vrot_3.append(AGy_vrot_tr[i])
# 		AGy_HImain_3.append(AGy_HImain_vrot_tr[i])
# 	if AGy_n[i]==4:
# 		AGy_r_4.append(AGy_r[i])
# 		AGy_vrot_4.append(AGy_vrot_tr[i])
# 		AGy_HImain_4.append(AGy_HImain_vrot_tr[i])
# 	if AGy_n[i]==5:
# 		AGy_r_5.append(AGy_r[i])
# 		AGy_vrot_5.append(AGy_vrot_tr[i])
# 		AGy_HImain_5.append(AGy_HImain_vrot_tr[i])

# AGo_r_1=[] #points where n=1
# AGo_vrot_1=[]
# AGo_HImain_1=[]
# AGo_r_2=[] #points where n=2
# AGo_vrot_2=[]
# AGo_HImain_2=[]
# AGo_r_3=[] #points where n=3
# AGo_vrot_3=[]
# AGo_HImain_3=[]
# AGo_r_4=[] #points where n=4
# AGo_vrot_4=[]
# AGo_HImain_4=[]
# AGo_r_5=[] #points where n=5
# AGo_vrot_5=[]
# AGo_HImain_5=[]
# for i in range(len(AGo_r)):
# 	if AGo_n[i]==1:
# 		AGo_r_1.append(AGo_r[i])
# 		AGo_vrot_1.append(AGo_vrot_tr[i])
# 		AGo_HImain_1.append(AGo_HImain_vrot_tr[i])
# 	if AGo_n[i]==2:
# 		AGo_r_2.append(AGo_r[i])
# 		AGo_vrot_2.append(AGo_vrot_tr[i])
# 		AGo_HImain_2.append(AGo_HImain_vrot_tr[i])
# 	if AGo_n[i]==3:
# 		AGo_r_3.append(AGo_r[i])
# 		AGo_vrot_3.append(AGo_vrot_tr[i])
# 		AGo_HImain_3.append(AGo_HImain_vrot_tr[i])
# 	if AGo_n[i]==4:
# 		AGo_r_4.append(AGo_r[i])
# 		AGo_vrot_4.append(AGo_vrot_tr[i])
# 		AGo_HImain_4.append(AGo_HImain_vrot_tr[i])
# 	if AGo_n[i]==5:
# 		AGo_r_5.append(AGo_r[i])
# 		AGo_vrot_5.append(AGo_vrot_tr[i])
# 		AGo_HImain_5.append(AGo_HImain_vrot_tr[i])

# RG_r_1=[] #points where n=1
# RG_vrot_1=[]
# RG_HImain_1=[]
# RG_r_2=[] #points where n=2
# RG_vrot_2=[]
# RG_HImain_2=[]
# RG_r_3=[] #points where n=3
# RG_vrot_3=[]
# RG_HImain_3=[]
# RG_r_4=[] #points where n=4
# RG_vrot_4=[]
# RG_HImain_4=[]
# RG_r_5=[] #points where n=5
# RG_vrot_5=[]
# RG_HImain_5=[]
# for i in range(len(RG_r)):
# 	if RG_n[i]==1:
# 		RG_r_1.append(RG_r[i])
# 		RG_vrot_1.append(RG_vrot_tr[i])
# 		RG_HImain_1.append(RG_HImain_vrot_tr[i])
# 	if RG_n[i]==2:
# 		RG_r_2.append(RG_r[i])
# 		RG_vrot_2.append(RG_vrot_tr[i])
# 		RG_HImain_2.append(RG_HImain_vrot_tr[i])
# 	if RG_n[i]==3:
# 		RG_r_3.append(RG_r[i])
# 		RG_vrot_3.append(RG_vrot_tr[i])
# 		RG_HImain_3.append(RG_HImain_vrot_tr[i])
# 	if RG_n[i]==4:
# 		RG_r_4.append(RG_r[i])
# 		RG_vrot_4.append(RG_vrot_tr[i])
# 		RG_HImain_4.append(RG_HImain_vrot_tr[i])
# 	if RG_n[i]==5:
# 		RG_r_5.append(RG_r[i])
# 		RG_vrot_5.append(RG_vrot_tr[i])
# 		RG_HImain_5.append(RG_HImain_vrot_tr[i])

# from functions import *
# AGo_r_med1, AGo_vrot_med1 = median_line(AGo_r_1, AGo_vrot_1)
# AGo_r_med2, AGo_vrot_med2 = median_line(AGo_r_2, AGo_vrot_2)
# AGo_r_med3, AGo_vrot_med3 = median_line(AGo_r_3, AGo_vrot_3)
# AGo_r_med4, AGo_vrot_med4 = median_line(AGo_r_4, AGo_vrot_4)
# AGo_r_med5, AGo_vrot_med5 = median_line(AGo_r_5, AGo_vrot_5)
# AGo_r_med1, AGo_HIvrot_med1 = median_line(AGo_r_1, AGo_HImain_1)
# AGo_r_med2, AGo_HIvrot_med2 = median_line(AGo_r_2, AGo_HImain_2)
# AGo_r_med3, AGo_HIvrot_med3 = median_line(AGo_r_3, AGo_HImain_3)
# AGo_r_med4, AGo_HIvrot_med4 = median_line(AGo_r_4, AGo_HImain_4)
# AGo_r_med5, AGo_HIvrot_med5 = median_line(AGo_r_5, AGo_HImain_5)

# AGo_AD1 = asymmetric_drift(AGo_vrot_1, AGo_HImain_1)
# AGo_AD2 = asymmetric_drift(AGo_vrot_2, AGo_HImain_2)
# AGo_AD3 = asymmetric_drift(AGo_vrot_3, AGo_HImain_3)
# AGo_AD4 = asymmetric_drift(AGo_vrot_4, AGo_HImain_4)
# AGo_AD5 = asymmetric_drift(AGo_vrot_5, AGo_HImain_5)

# rc('font', family = 'serif')
# f, axes= plt.subplots(2,3, sharey=True, sharex=False, figsize=(10,5))
# axes[0,0].scatter(AGo_r_1, AGo_HImain_1, s=2, c='darkgray')
# axes[0,0].scatter(AGo_r_1, AGo_vrot_1, s=2, c='green', alpha=0.4)
# axes[0,1].scatter(AGo_r_2, AGo_HImain_2, s=2, c='darkgray')
# axes[0,1].scatter(AGo_r_2, AGo_vrot_2, s=2, c='green', alpha=0.4)
# axes[0,2].scatter(AGo_r_3, AGo_HImain_3, s=2, c='darkgray')
# axes[0,2].scatter(AGo_r_3, AGo_vrot_3, s=2, c='green', alpha=0.4)
# axes[1,0].scatter(AGo_r_4, AGo_HImain_4, s=2, c='darkgray')
# axes[1,0].scatter(AGo_r_4, AGo_vrot_4, s=2, c='green', alpha=0.4)
# axes[1,1].scatter(AGo_r_5, AGo_HImain_5, s=2, c='darkgray')
# axes[1,1].scatter(AGo_r_5, AGo_vrot_5, s=2, c='green', alpha=0.4)
# axes[0,0].annotate(r'$N_{\rm HI}=1$', xy=(5,115), horizontalalignment='left', fontsize=12)
# axes[0,1].annotate(r'$N_{\rm HI}=2$', xy=(5,115), horizontalalignment='left', fontsize=12)
# axes[0,2].annotate(r'$N_{\rm HI}=3$', xy=(5,115), horizontalalignment='left', fontsize=12)
# axes[1,0].annotate(r'$N_{\rm HI}=4$', xy=(5,115), horizontalalignment='left', fontsize=12)
# axes[1,1].annotate(r'$N_{\rm HI}=5$', xy=(5,115), horizontalalignment='left', fontsize=12)

# axes[0,0].plot(AGo_r_med1, AGo_vrot_med1, linestyle='-', c='black', linewidth = 1.8)
# axes[0,0].plot(AGo_r_med1, AGo_HIvrot_med1, linestyle='--', c='black', linewidth = 1.8)
# axes[0,1].plot(AGo_r_med2, AGo_vrot_med2, linestyle='-', c='black', linewidth = 1.8)
# axes[0,1].plot(AGo_r_med2, AGo_HIvrot_med2, linestyle='--', c='black', linewidth = 1.8)
# axes[0,2].plot(AGo_r_med3, AGo_vrot_med3, linestyle='-', c='black', linewidth = 1.8)
# axes[0,2].plot(AGo_r_med3, AGo_HIvrot_med3, linestyle='--', c='black', linewidth = 1.8)
# axes[1,0].plot(AGo_r_med4, AGo_vrot_med4, linestyle='-', c='black', linewidth = 1.8)
# axes[1,0].plot(AGo_r_med4, AGo_HIvrot_med4, linestyle='--', c='black', linewidth = 1.8)
# axes[1,1].plot(AGo_r_med5, AGo_vrot_med5, linestyle='-', c='black', linewidth = 1.8)
# axes[1,1].plot(AGo_r_med5, AGo_HIvrot_med5, linestyle='--', c='black', linewidth = 1.8)

# axes[0,0].annotate(r'$\overline{v}_{\rm a}=$'+ '${}$'.format(round(np.median(AGo_AD1),2)) + r'$\rm km\ s^{-1}$', xy=(10,115), horizontalalignment='left', fontsize=12)
# axes[0,1].annotate(r'$\overline{v}_{\rm a}=$'+ '${}$'.format(round(np.median(AGo_AD2),2)) + r'$\rm km\ s^{-1}$', xy=(10,115), horizontalalignment='left', fontsize=12)
# axes[0,2].annotate(r'$\overline{v}_{\rm a}=$'+ '${}$'.format(round(np.median(AGo_AD3),2)) + r'$\rm km\ s^{-1}$', xy=(10,115), horizontalalignment='left', fontsize=12)
# axes[1,0].annotate(r'$\overline{v}_{\rm a}=$'+ '${}$'.format(round(np.median(AGo_AD4),2)) + r'$\rm km\ s^{-1}$', xy=(10,115), horizontalalignment='left', fontsize=12)
# axes[1,1].annotate(r'$\overline{v}_{\rm a}=$'+ '${}$'.format(round(np.median(AGo_AD5),2)) + r'$\rm km\ s^{-1}$', xy=(10,115), horizontalalignment='left', fontsize=12)


# for ax in axes[0,:]:
# 	ax.set_xlim(4, 20)
# 	ax.tick_params(axis='x',which='both',bottom='on',top='off', direction='out')
# 	ax.tick_params(axis='x',which='both',top='on', direction='in')
# 	ax.tick_params(which='both', width=2)
# 	ax.tick_params(which='major', length=7)
# 	ax.tick_params(which='minor', length=4)
# 	ax.tick_params(labelsize=12) 
# 	ax.minorticks_on()
# 	for axis in ['top','bottom','left','right']:
# 	        ax.spines[axis].set_linewidth(2)

# for ax in axes[1,:]:
# 	ax.set_xlim(4, 20)
# 	ax.set_xlabel(r'$r \ \rm(kpc)$', fontsize=13)
# 	ax.tick_params(axis='x',which='both',bottom='on',top='off', direction='out')
# 	ax.tick_params(axis='x',which='both',top='on', direction='in')
# 	ax.tick_params(which='both', width=2)
# 	ax.tick_params(which='major', length=7)
# 	ax.tick_params(which='minor', length=4)
# 	ax.tick_params(labelsize=12) 
# 	ax.minorticks_on()
# 	for axis in ['top','bottom','left','right']:
# 	        ax.spines[axis].set_linewidth(2)

# for ax in axes[:,0]:
# 	#ax.set_ylabel(r'$ v_{\rm rot} \ \rm(km\ s^{-1})$', fontsize=13)
# 	ax.set_ylim(100,300)
# 	ax.tick_params(axis='y',which='both',left='on',top='off', direction='out')
# 	ax.tick_params(axis='y',which='both',right='on', direction='in')
# 	ax.tick_params(which='both', width=2)
# 	ax.tick_params(which='major', length=7)
# 	ax.tick_params(which='minor', length=4)
# 	ax.tick_params(labelsize=12) 
# 	ax.minorticks_on()
# 	for axis in ['top','bottom','left','right']:
# 	        ax.spines[axis].set_linewidth(2)

# for ax in axes[:,1]:
# 	#ax.set_ylabel(r'$\rm v_{\rm rot} \ (km\ s^{-1})$', fontsize=13)
# 	ax.set_ylim(100,300)
# 	ax.tick_params(axis='y',which='both',left='on',top='off', direction='out')
# 	ax.tick_params(axis='y',which='both',right='on', direction='in')
# 	ax.tick_params(which='both', width=2)
# 	ax.tick_params(which='major', length=7)
# 	ax.tick_params(which='minor', length=4)
# 	ax.tick_params(labelsize=12) 
# 	ax.minorticks_on()
# 	for axis in ['top','bottom','left','right']:
# 	        ax.spines[axis].set_linewidth(2)

# for ax in axes[:,2]:
# 	#ax.set_ylabel(r'$\rm v_{\rm rot} \ (km\ s^{-1})$', fontsize=13)
# 	ax.set_ylim(100,300)
# 	ax.tick_params(axis='y',which='both',left='on',top='off', direction='out')
# 	ax.tick_params(axis='y',which='both',right='on', direction='in')
# 	ax.tick_params(which='both', width=2)
# 	ax.tick_params(which='major', length=7)
# 	ax.tick_params(which='minor', length=4)
# 	ax.tick_params(labelsize=12) 
# 	ax.minorticks_on()
# 	for axis in ['top','bottom','left','right']:
# 	        ax.spines[axis].set_linewidth(2)

# axes[-1,-1].axis('off')
# axes[0,1].tick_params(labelsize=0)
# axes[0,2].tick_params(labelsize=0)
# nbins = 4
# axes[0,0].yaxis.set_major_locator(MaxNLocator(nbins=nbins, prune=None))
# axes[0, 0].set_yticks([100, 150, 200, 250, 300])
# axes[1,0].yaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='upper'))

# nbins= 6
# axes[0,0].xaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='upper'))
# axes[1,0].xaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='upper'))
# axes[0,1].xaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='upper'))
# axes[1,1].xaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='upper'))
# axes[0,2].xaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='upper'))

# f.subplots_adjust(left=0.17)
# f.text(0.10, 0.5, r'$\rm Rotation\ Velocity:\ \it v_{\rm rot}\ \rm(km\ s^{-1})$', va='center', rotation='vertical', fontsize=13)

# plt.subplots_adjust(wspace=0, hspace=0)
#plt.savefig('/Users/amandaquirk/Desktop/AGo_rotation_curves.pdf', bbox_inches='tight')

# rotation_curve(MS_r_1, MS_vrot_1, MS_HImain_1, 'MS', 1)
# rotation_curve(MS_r_2, MS_vrot_2, MS_HImain_2, 'MS', 2)
# rotation_curve(MS_r_3, MS_vrot_3, MS_HImain_3, 'MS', 3)
# rotation_curve(MS_r_4, MS_vrot_4, MS_HImain_4, 'MS', 4)
# rotation_curve(MS_r_5, MS_vrot_5, MS_HImain_5, 'MS', 5)
# rotation_curve(AGy_r_1, AGy_vrot_1, AGy_HImain_1, 'AGy', 1)
# rotation_curve(AGy_r_2, AGy_vrot_2, AGy_HImain_2, 'AGy', 2)
# rotation_curve(AGy_r_3, AGy_vrot_3, AGy_HImain_3, 'AGy', 3)
# rotation_curve(AGy_r_4, AGy_vrot_4, AGy_HImain_4, 'AGy', 4)
# rotation_curve(AGy_r_5, AGy_vrot_5, AGy_HImain_5, 'AGy', 5)
# rotation_curve(AGo_r_1, AGo_vrot_1, AGo_HImain_1, 'AGo', 1)
# rotation_curve(AGo_r_2, AGo_vrot_2, AGo_HImain_2, 'AGo', 2)
# rotation_curve(AGo_r_3, AGo_vrot_3, AGo_HImain_3, 'AGo', 3)
# rotation_curve(AGo_r_4, AGo_vrot_4, AGo_HImain_4, 'AGo', 4)
# rotation_curve(AGo_r_5, AGo_vrot_5, AGo_HImain_5, 'AGo', 5)
# rotation_curve(RG_r_1, RG_vrot_1, RG_HImain_1, 'RG', 1)
# rotation_curve(RG_r_2, RG_vrot_2, RG_HImain_2, 'RG', 2)
# rotation_curve(RG_r_3, RG_vrot_3, RG_HImain_3, 'RG', 3)
# rotation_curve(RG_r_4, RG_vrot_4, RG_HImain_4, 'RG', 4)
# rotation_curve(RG_r_5, RG_vrot_5, RG_HImain_5, 'RG', 5)
