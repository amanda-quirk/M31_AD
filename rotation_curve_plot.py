import matplotlib.pyplot as plt
import numpy as np 
from matplotlib import rc 
from matplotlib.ticker import MaxNLocator


MS_r, MS_vrot_tr, MS_HImain_vrot_tr, MS_HIclose_vrot_tr=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/MS_master_vrot.txt', usecols=(0,2,8,10), unpack=True)
AGy_r, AGy_vrot_tr, AGy_HImain_vrot_tr, AGy_HIclose_vrot_tr=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/AGy_master_vrot.txt', usecols=(0,2,8,10), unpack=True)
AGo_r, AGo_vrot_tr, AGo_HImain_vrot_tr, AGo_HIclose_vrot_tr=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/AGo_master_vrot.txt', usecols=(0,2,8,10), unpack=True)
RG_r, RG_vrot_tr, RG_HImain_vrot_tr, RG_HIclose_vrot_tr=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/RG_master_vrot.txt', usecols=(0,2,8,10), unpack=True)

rc('font', family = 'serif')
f, axes= plt.subplots(4,1, sharey=False, sharex=True, figsize=(4,9.8))
axes[0].scatter(MS_r, MS_HImain_vrot_tr, s=2, c='darkgray')
axes[0].scatter(MS_r, MS_vrot_tr, s=2, c='b', alpha=0.4)
axes[1].scatter(AGy_r, AGy_HImain_vrot_tr, s=2, c='darkgray')
axes[1].scatter(AGy_r, AGy_vrot_tr, s=2, c='m', alpha=0.4)
axes[2].scatter(AGo_r, AGo_HImain_vrot_tr, s=2, c='darkgray')
axes[2].scatter(AGo_r, AGo_vrot_tr, s=2, c='k', alpha=0.4)
axes[3].scatter(RG_r, RG_HImain_vrot_tr, s=2, c='darkgray')
axes[3].scatter(RG_r, RG_vrot_tr, s=2, c='r', alpha=0.4)
axes[0].annotate('MS', xy=(19,115), horizontalalignment='right', fontsize=12)
axes[1].annotate('young AGB', xy=(19,115), horizontalalignment='right', fontsize=12)
axes[2].annotate('older AGB', xy=(19,115), horizontalalignment='right', fontsize=12)
axes[3].annotate('RGB', xy=(19,115), horizontalalignment='right', fontsize=12)

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
f.text(0.008, 0.5, r'$\rm Rotation\ Velocity:\ \it v_{\rm rot}\ \rm(km\ s^{-1})$', va='center', rotation='vertical', fontsize=13)
plt.subplots_adjust(wspace=0, hspace=0)
plt.savefig('/Users/amandaquirk/Desktop/rotation_curves.pdf', bbox_inches='tight')