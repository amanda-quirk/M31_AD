import numpy as np 
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import rc,rcParams
from matplotlib.ticker import MaxNLocator
from matplotlib import patches

#import rotation curve data
#r (kpc), v_rot no model (km/s), v_rot tr model, smoothed PA, sPA v_rot no model, sPA v_rot tr model, n, HI main no model, HI main tr, HI close no model, HI close tr
MS_r, MS_vrot_tr, MS_HImain_vrot_tr, MS_HIclose_vrot_tr=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/MS_master_vrot.txt', usecols=(0,2,8,10), unpack=True)
AGy_r, AGy_vrot_tr, AGy_HImain_vrot_tr, AGy_HIclose_vrot_tr=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/AGy_master_vrot.txt', usecols=(0,2,8,10), unpack=True)
AGo_r, AGo_vrot_tr, AGo_HImain_vrot_tr, AGo_HIclose_vrot_tr=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/AGo_master_vrot.txt', usecols=(0,2,8,10), unpack=True)
RG_r, RG_vrot_tr, RG_HImain_vrot_tr, RG_HIclose_vrot_tr=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/RG_master_vrot.txt', usecols=(0,2,8,10), unpack=True)

MS_xi, MS_eta, MS_v, MS_dispersion, MS_n, MS_HImain, MS_HIclose=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/MS_smoothed_chemin.txt', usecols=(0,1,2,4,5,6,7,), unpack=True)
AGy_xi, AGy_eta, AGy_v, AGy_dispersion,AGy_n, AGy_HImain, AGy_HIclose=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/AGy_smoothed_chemin.txt', usecols=(0,1,2,4,5,6,7,), unpack=True)
AGo_xi, AGo_eta, AGo_v, AGo_dispersion, AGo_n, AGo_HImain, AGo_HIclose=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/AGo_smoothed_chemin.txt', usecols=(0,1,2,4,5,6,7,), unpack=True)
RG_xi, RG_eta, RG_v, RG_dispersion, RG_n, RG_HImain, RG_HIclose=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/RG_smoothed_chemin.txt', usecols=(0,1,2,4,5,6,7,), unpack=True)

MS_Av=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/MS_dust.txt', usecols=(2,), unpack=True)
AGy_Av=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/AGy_dust.txt', usecols=(2,), unpack=True)
AGo_Av=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/AGo_dust.txt', usecols=(2,), unpack=True)
RG_Av=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/RG_dust.txt', usecols=(2,), unpack=True)

MS_ad_CO= np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/MS_CO_ad.txt', usecols=(2,), unpack= True)
AGy_ad_CO= np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/AGy_CO_ad.txt', usecols=(2,), unpack= True)
AGo_ad_CO= np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/AGo_CO_ad.txt', usecols=(2,), unpack= True)
RG_ad_CO= np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/RG_CO_ad.txt', usecols=(2,), unpack= True)

def asymmetric_drift(v_star, v_gas):
	return v_gas-v_star

MS_ad=asymmetric_drift(MS_vrot_tr, MS_HImain_vrot_tr)
AGy_ad=asymmetric_drift(AGy_vrot_tr, AGy_HImain_vrot_tr)
AGo_ad=asymmetric_drift(AGo_vrot_tr, AGo_HImain_vrot_tr)
RG_ad=asymmetric_drift(RG_vrot_tr, RG_HImain_vrot_tr)

MS_range=max(MS_ad)-min(MS_ad)
AGy_range=max(AGy_ad)-min(AGy_ad)
AGo_range=max(AGo_ad)-min(AGo_ad)
RG_range=max(RG_ad)-min(RG_ad)

rc('font', family = 'serif')
fig, ax=plt.subplots(1)
for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(2)
ax.tick_params(axis='x',which='both',bottom='on',top='off', direction='out')
ax.tick_params(axis='x',which='both',top='on', direction='in')
ax.tick_params(axis='y',which='both',left='on',top='off', direction='out')
ax.tick_params(axis='y',which='both',right='on', direction='in')
plt.tick_params(which='both', width=2)
plt.tick_params(which='major', length=7)
plt.tick_params(which='minor', length=4)
plt.tick_params(labelsize=12) 
plt.minorticks_on()
plt.hist(MS_ad_CO, bins=range(-200, 300, 20), label='MS={} stars'.format(len(MS_ad_CO)),normed=1, histtype='step', linewidth=1.6,linestyle='--',stacked=True,fill=False, color='b')
plt.hist(AGy_ad_CO,  bins=range(-200, 300, 20),label='young AGB={} stars'.format(len(AGy_ad_CO)),normed=1, histtype='step', linewidth=1.6,stacked=True,fill=False, color='m')
plt.hist(AGo_ad_CO,  bins=range(-200, 300, 20),label='old AGB={} stars'.format(len(AGo_ad_CO)),normed=1,histtype='step', linewidth=1.6,stacked=True,fill=False, color='k')
plt.hist(RG_ad_CO, bins=range(-200, 300, 20), label='RGB={} stars'.format(len(RG_ad_CO)),normed=1,histtype='step', linewidth=1.6,linestyle='--',stacked=True,fill=False, color='r')
plt.legend(loc=1, frameon=False)
plt.xlim(-200,300)
plt.xlabel(r'$ \rm Asymmetric\ Drift:\ \itv_{a}\ \rm(km\ s^{-1})$', fontsize=13)
plt.savefig('/Users/amandaquirk/Desktop/CO_asymmetric_drift_zoom.pdf')

def lag_n_plt(n, lag, age):
	plt.scatter(n, lag, alpha=.3)
	plt.ylim(-200,200)
	plt.ylabel('Asymmetric Drift (km/s)')
	plt.xlabel('Number of HI Peaks')
	plt.savefig('/Users/amandaquirk/Desktop/{}_lag_n.png'.format(age))
	plt.close()

def dispersion_lag_plt(lag, dispersion, age):
	plt.scatter(dispersion, lag, alpha=.3)
	plt.xlim(10,150)
	plt.ylim(-200,200)
	plt.xlabel('Dispersion (km/s)')
	plt.ylabel('Asymmetric Drift (km/s)')
	plt.savefig('/Users/amandaquirk/Desktop/{}_var_lag.png'.format(age))
	plt.close()

def lag_r_n_plt(n, lag, radius, age):
	plt.ylim(-200,200)
	plt.xlim(0,20)
	plt.ylabel('Asymmetric Drift (km/s)')
	plt.xlabel('Radius [kpc]')
	cmap=plt.cm.rainbow
	norm = matplotlib.colors.BoundaryNorm([0.5,1.5,2.5,3.5,4.5,5.5], cmap.N)
	plt.scatter(radius, lag, alpha=.7, c=n,cmap=cmap,norm=norm,edgecolor='none')
	clb=plt.colorbar(ticks=np.linspace(1,5,5))
	clb.set_label('Number of HI Peaks')
	plt.savefig('/Users/amandaquirk/Desktop/{}_lag_r_n.png'.format(age))
	plt.close()


# rc('font', family = 'serif')
# cmap=plt.cm.rainbow
# norm = matplotlib.colors.BoundaryNorm([0.5,1.5,2.5,3.5,4.5,5.5], cmap.N)
# f, axes= plt.subplots(1,4, sharex=False, sharey=True, figsize=(15,3))
# im=axes[0].scatter(MS_r, MS_ad, c=MS_n, cmap=cmap, norm=norm, edgecolor='none')
# axes[1].scatter(AGy_r, AGy_ad, c=AGy_n, cmap=cmap, norm=norm, edgecolor='none')
# axes[2].scatter(AGo_r, AGo_ad, c=AGo_n, cmap=cmap, norm=norm, edgecolor='none')
# axes[3].scatter(RG_r, RG_ad, c=RG_n, cmap=cmap, norm=norm, edgecolor='none')
# axes[0].annotate('MS', xy=(1,-180), horizontalalignment='left', fontsize=12)
# axes[1].annotate('young AGB', xy=(1,-180), horizontalalignment='left', fontsize=12)
# axes[2].annotate('older AGB', xy=(1,-180), horizontalalignment='left', fontsize=12)
# axes[3].annotate('RGB', xy=(1,-180), horizontalalignment='left', fontsize=12)

# for ax in axes:
# 	ax.set_xlim(0, 20)
# 	#ax.set_xlabel(r'$r \ \rm(kpc)$', fontsize=13)
# 	ax.set_ylim(-200,200)
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
# axes[0].set_ylabel(r'$\rm Asymmetric\ Drift:\ \it v_{a} \ \rm(km\ s^{-1})$', fontsize=13)
# nbins = len(axes[0].get_xticklabels())-1
# axes[0].xaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='upper'))
# axes[1].xaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='upper'))
# axes[2].xaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='upper'))
# f.subplots_adjust(right=0.86)
# cax= f.add_axes([0.87, 0.13, 0.009, 0.72])
# clb=f.colorbar(im, cax=cax, ticks=np.linspace(1,5,5))
# clb.set_label(r'$\rm Multiplicity\ of\ HI:\ \itN_{\rm HI}$', fontsize=13)

# f.subplots_adjust(bottom=0.13)
# f.text(0.5, 0.008, r'$\rm Radial\ Distance:\ \it r \ \rm(kpc)$', ha='center', fontsize=13)

# plt.subplots_adjust(wspace=0, hspace=0)
# plt.savefig('/Users/amandaquirk/Desktop/lag_r_n.pdf', bbox_inches='tight')

# def dust_lag_plt(dust, lag, age):
# 	plt.scatter(dust,lag, alpha=.3)
# 	plt.ylim(-200,200)
# 	plt.ylabel('Asymmetric Drift (km/s)')
# 	plt.xlabel('Extinction')
# 	plt.savefig('/Users/amandaquirk/Desktop/dust_lag_{}'.format(age))
# 	plt.close()

# def dust_dispersion_plt(dust, dispersion, age):
# 	plt.scatter(dust,dispersion, alpha=.3)
# 	plt.ylim(10,150)
# 	plt.ylabel('Dispersion (km/s)')
# 	plt.xlabel('Extinction')
# 	plt.savefig('/Users/amandaquirk/Desktop/dust_dispersion_{}'.format(age))
# 	plt.close()

#rc('font', family = 'serif')
#f, axes= plt.subplots(1,4, sharex=False, sharey=True, figsize=(15,3))
#im=axes[0].scatter(MS_Av, MS_dispersion, c='b',alpha=0.3)
#axes[1].scatter(AGy_Av, AGy_dispersion, c='m',alpha=0.3)
#axes[2].scatter(AGo_Av, AGo_dispersion, c='k',alpha=0.3)
#axes[3].scatter(RG_Av, RG_dispersion, c='r',alpha=0.3)

#for ax in axes:
#	ax.set_xlim(0,5)
#	ax.set_ylim(10,150)
#	ax.tick_params(axis='x',which='both',bottom='on',top='off', direction='out')
#	ax.tick_params(axis='x',which='both',top='on', direction='in')
#	ax.tick_params(axis='y',which='both',left='on',top='off', direction='out')
#	ax.tick_params(axis='y',which='both',right='on', direction='in')
#	ax.tick_params(which='both', width=2)
#	ax.tick_params(which='major', length=7)
#	ax.tick_params(which='minor', length=4)
#	ax.tick_params(labelsize=12) 
#	ax.minorticks_on()
#	for axis in ['top','bottom','left','right']:
#	        ax.spines[axis].set_linewidth(2)
#axes[0].annotate('MS', xy=(4.8,20), horizontalalignment='right', fontsize=12) #-175 for va, 20 for sigma
#axes[1].annotate('young AGB', xy=(4.8,20), horizontalalignment='right', fontsize=12)
#axes[2].annotate('older AGB', xy=(4.8,20), horizontalalignment='right', fontsize=12)
#axes[3].annotate('RGB', xy=(4.8,20), horizontalalignment='right', fontsize=12)
#nbins = len(axes[0].get_xticklabels())
#axes[0].xaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='upper'))
#axes[1].xaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='upper'))
#axes[2].xaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='upper'))
#axes[3].xaxis.set_major_locator(MaxNLocator(nbins=nbins))	        
#axes[0].set_ylabel(r'$\rm Velocity\ Dispersion:\ \sigma\ (km\ s^{-1})$', fontsize=13, labelpad=15)
#axes[0].set_ylabel(r'$\rm Asymmetric\ Drift:\ \itv_{a}\ \rm(km\ s^{-1})$', fontsize=13)
#f.subplots_adjust(bottom=0.13)
#f.text(0.5, 0.008, r'$\rm Dust\ Extinction:\ \itA_{v}$', ha='center', fontsize=13)
#plt.subplots_adjust(wspace=0, hspace=0)
#plt.savefig('/Users/amandaquirk/Desktop/disp_dust.pdf', bbox_inches='tight')

# lag_n_plt(MS_n, MS_ad, 'MS')
# lag_n_plt(AGy_n, AGy_ad, 'AGy')
# lag_n_plt(AGo_n, AGo_ad, 'AGo')
# lag_n_plt(RG_n, RG_ad, 'RG')
# dispersion_lag_plt(MS_ad, MS_dispersion, 'MS')
# dispersion_lag_plt(AGy_ad, AGy_dispersion, 'AGy')
# dispersion_lag_plt(AGo_ad, AGo_dispersion, 'AGo')
# dispersion_lag_plt(RG_ad, RG_dispersion, 'RG')
# lag_r_n_plt(MS_n, MS_ad, MS_r, 'MS')
# lag_r_n_plt(AGy_n, AGy_ad, AGy_r, 'AGy')
# lag_r_n_plt(AGo_n, AGo_ad, AGo_r, 'AGo')
# lag_r_n_plt(RG_n, RG_ad, RG_r, 'RG')
# dust_lag_plt(MS_Av,MS_ad ,'MS')
# dust_lag_plt(AGy_Av,AGy_ad, 'AGy')
# dust_lag_plt(AGo_Av,AGo_ad, 'AGo')
# dust_lag_plt(RG_Av,RG_ad, 'RG')
# dust_dispersion_plt(MS_Av,MS_dispersion ,'MS')
# dust_dispersion_plt(AGy_Av,AGy_dispersion, 'AGy')
# dust_dispersion_plt(AGo_Av,AGo_dispersion, 'AGo')
# dust_dispersion_plt(RG_Av,RG_dispersion, 'RG')

# file=open('/Users/amandaquirk/Desktop/MS_ad.txt', 'w')
# file.write('#ad (km/s)\n')
# for i in range(len(MS_ad)):
# 	file.write('{}\n'.format(MS_ad[i]))
# file.close()

# file=open('/Users/amandaquirk/Desktop/AGy_ad.txt', 'w')
# file.write('#ad (km/s)\n')
# for i in range(len(AGy_ad)):
# 	file.write('{}\n'.format(AGy_ad[i]))
# file.close()

# file=open('/Users/amandaquirk/Desktop/AGo_ad.txt', 'w')
# file.write('#ad (km/s)\n')
# for i in range(len(AGo_ad)):
# 	file.write('{}\n'.format(AGo_ad[i]))
# file.close()

# file=open('/Users/amandaquirk/Desktop/RG_ad.txt', 'w')
# file.write('#ad (km/s)\n')
# for i in range(len(RG_ad)):
# 	file.write('{}\n'.format(RG_ad[i]))
# file.close()

rc('font', family = 'serif')
f, axes= plt.subplots(1,4, sharey=True, sharex=False, figsize=(11.1,3.42))
ylength=37.5
xlength=ylength*np.cos(float(73.7*np.pi) / 180)
e =  patches.Ellipse((0, 0), xlength, ylength, angle=180-37.1, linewidth=.5, fill=False, zorder=0)
e1 = patches.Ellipse((0, 0), xlength, ylength, angle=180-37.1, linewidth=.5, fill=False, zorder=0)
e2 = patches.Ellipse((0, 0), xlength, ylength, angle=180-37.1, linewidth=.5, fill=False, zorder=0)
e3 = patches.Ellipse((0, 0), xlength, ylength, angle=180-37.1, linewidth=.5, fill=False, zorder=0)
axes[3].add_patch(e3)
axes[2].add_patch(e2)
axes[1].add_patch(e1)
axes[0].add_patch(e)
axes[0].scatter(MS_xi,MS_eta, c=MS_ad, cmap='plasma', s=7, vmin=-100, vmax=150) 
axes[1].scatter(AGy_xi,AGy_eta, c=AGy_ad, cmap='plasma', s=7, vmin=-100, vmax=150) 
axes[2].scatter(AGo_xi,AGo_eta, c=AGo_ad, cmap='plasma', s=7, vmin=-100, vmax=150) 
im=axes[3].scatter(RG_xi,RG_eta, c=RG_ad, cmap='plasma', s=7, vmin=-100, vmax=150)
for ax in axes:
	ax.set_xlim(13, -1.5)
	ax.set_xlabel(r'$\xi\ (kpc)$', fontsize=13)
	ax.set_ylim(-2.5,15.5)
	ax.tick_params(axis='x',which='both',bottom='on',top='off', direction='out')
	ax.tick_params(axis='x',which='both',top='on', direction='in')
	ax.tick_params(axis='y',which='both',left='on',top='off', direction='out')
	ax.tick_params(axis='y',which='both',right='on', direction='in')
	ax.tick_params(which='both', width=1)
	ax.tick_params(which='major', length=7)
	ax.tick_params(which='minor', length=4)
	ax.tick_params(labelsize=12) 
	ax.minorticks_on()
	for axis in ['top','bottom','left','right']:
		ax.spines[axis].set_linewidth(1)
axes[0].set_ylabel(r'$\eta\ (kpc)$', fontsize=13)
f.subplots_adjust(right=0.885)
cbar_ax = f.add_axes([0.89, 0.13, 0.015, 0.72])
clb=f.colorbar(im, cax=cbar_ax)
clb.set_label(r'$\rm Asymmetric\ Drift:\ v_{a}$', fontsize=13)
axes[0].scatter(0,0,marker='+', c='b', linewidth=2)
axes[1].scatter(0,0,marker='+', c='b', linewidth=2)
axes[2].scatter(0,0,marker='+', c='b', linewidth=2)
axes[3].scatter(0,0,marker='+', c='b', linewidth=2)
axes[0].annotate('MS', xy=(-.55,13.5), horizontalalignment='right', fontsize=12)
axes[1].annotate('young AGB', xy=(-.55,13.5), horizontalalignment='right', fontsize=12)
axes[2].annotate('older AGB', xy=(-.55,13.5), horizontalalignment='right', fontsize=12)
axes[3].annotate('RGB', xy=(-.55,13.5), horizontalalignment='right', fontsize=12)
plt.subplots_adjust(wspace=0, hspace=0)
plt.savefig('/Users/amandaquirk/Desktop/ad_maps.pdf', bbox_inches='tight')
