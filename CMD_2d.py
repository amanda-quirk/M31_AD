import matplotlib.pyplot as plt 
from matplotlib.colors import LogNorm
import numpy as np
from astropy.io import fits
from astropy.table import Table
from matplotlib import rc 

#data
F110W, F160W=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/RG_color.txt', usecols=(0,1), unpack=True)
RG_r, RG_vrot_tr, RG_HImain_vrot_tr, RG_HIclose_vrot_tr=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/RG_master_vrot.txt', usecols=(0,2,8,10), unpack=True)
RG_Av=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/RG_dust.txt', usecols=(2,), unpack=True)


#get data for RGB stars
color=[]
mag=[]
r=[]
vrot=[]
HI_vrot=[]
Av=[]
for i in range(len(F110W)):
	if np.isnan(F110W[i])==False and np.isnan(F160W[i])==False:# and F160W[i]>18.1 and F160W[i]<21:
		color.append(float(F110W[i])-float(F160W[i]))
		mag.append(float(F160W[i]))
		r.append(RG_r[i])
		vrot.append(RG_vrot_tr[i])
		HI_vrot.append(RG_HImain_vrot_tr[i])
		Av.append(RG_Av[i])

#plotting
# rc('font', family = 'serif')
# fig, ax=plt.subplots(1, figsize=(4,6))
# for axis in ['top','bottom','left','right']:
#         ax.spines[axis].set_linewidth(2)
# ax.tick_params(axis='x',which='both',bottom='on',top='off', direction='out')
# ax.tick_params(axis='x',which='both',top='on', direction='in')
# ax.tick_params(axis='y',which='both',left='on',top='off', direction='out')
# ax.tick_params(axis='y',which='both',right='on', direction='in')
# plt.tick_params(which='both', width=2)
# plt.tick_params(which='major', length=7)
# plt.tick_params(which='minor', length=4)
# plt.tick_params(labelsize=12) 
# plt.minorticks_on()
# im=ax.hist2d(color, mag, bins=900, cmap=plt.cm.jet,  norm=LogNorm())
# x= np.linspace(0,1.5,10)
# y= -20*x + 38.5
# y1=-20*x + 36.5
# y2=-20*x + 40.5
# plt.colorbar(im[3], ax=ax)
# ax.plot(x,y, color='w', linewidth=3.5)
# ax.plot(x,y1, color='w', linewidth=3.5)
# ax.plot(x,y2, color='w', linewidth=3.5)
# ax.plot(x,y, color='k', linewidth=2)
# ax.plot(x,y1, color='k', linewidth=2, linestyle='--')
# ax.plot(x,y2, color='k', linewidth=2, linestyle='--')
# plt.xlim(0.45,1.5)
# plt.xlabel(r'$\rm F110W-F160W$', fontsize=13)
# plt.ylim(17,22)
# plt.gca().invert_yaxis()
# plt.ylabel(r'$\rm F160W$', fontsize=13)
# ax.annotate(r'$\bf I$', xy=(.6, 21.8), horizontalalignment='right', fontsize=14, color='hotpink') 
# ax.annotate(r'$\bf II$', xy=(.88, 21.9), horizontalalignment='right', fontsize=26, color='w')
# ax.annotate(r'$\bf II$', xy=(.88, 21.8), horizontalalignment='right', fontsize=14, color='hotpink')
# ax.annotate(r'$\bf III$', xy=(1.3,21.8), horizontalalignment='right', fontsize=14, color='hotpink')
# plt.savefig('/Users/amandaquirk/Desktop/RGB_CMD_scatter.pdf')

#diving the pop into 3 different regions
color1=[]
mag1=[]
r1=[]
vrot1=[]
HI_vrot1=[]
Av1=[]
color2=[]
mag2=[]
r2=[]
vrot2=[]
HI_vrot2=[]
Av2=[]
color3=[]
mag3=[]
r3=[]
vrot3=[]
HI_vrot3=[]
Av3=[]
for i in range(len(r)):
	if mag[i]<-20*color[i]+36.5:
		color1.append(color[i])
		mag1.append(mag[i])
		r1.append(r[i])
		vrot1.append(vrot[i])
		HI_vrot1.append(HI_vrot[i])
		Av1.append(Av[i])
	elif mag[i]>-20*color[i]+40.5:
		color3.append(color[i])
		mag3.append(mag[i])
		r3.append(r[i])
		vrot3.append(vrot[i])
		HI_vrot3.append(HI_vrot[i])
		Av3.append(Av[i])
	else:
		color2.append(color[i])
		mag2.append(mag[i])
		r2.append(r[i])
		vrot2.append(vrot[i])
		HI_vrot2.append(HI_vrot[i])
		Av2.append(Av[i])

# plt.figure(figsize=(4, 6))
# plt.scatter(color1,mag1, c='r')
# plt.scatter(color2,mag2, c='b')
# plt.scatter(color3,mag3, c='k')
# plt.plot(x,y, color='k', linewidth=2)
# plt.plot(x,y1, color='k', linewidth=2, linestyle='--')
# plt.plot(x,y2, color='k', linewidth=2, linestyle='--')
# plt.xlim(0,1.5)
# plt.ylim(17,22)
# plt.gca().invert_yaxis()
# plt.show()

#exploring the different regions
# rc('font', family = 'serif')
# fig, ax=plt.subplots(1)
# for axis in ['top','bottom','left','right']:
#         ax.spines[axis].set_linewidth(2)
# ax.tick_params(axis='x',which='both',bottom='on',top='off', direction='out')
# ax.tick_params(axis='x',which='both',top='on', direction='in')
# ax.tick_params(axis='y',which='both',left='on',top='off', direction='out')
# ax.tick_params(axis='y',which='both',right='on', direction='in')
# plt.tick_params(which='both', width=2)
# plt.tick_params(which='major', length=7)
# plt.tick_params(which='minor', length=4)
# plt.tick_params(labelsize=12) 
# plt.minorticks_on()
# weights1 = np.ones_like(Av1)/float(len(Av1))
# weights2 = np.ones_like(Av2)/float(len(Av2))
# weights3 = np.ones_like(Av3)/float(len(Av3))
# plt.hist(Av1, weights=weights1, bins=20, range=(0, 5),label=r'$\rm Region\ I$', normed=False, histtype='step', linewidth=1.6,stacked=True,fill=False, color='b')
# plt.hist(Av2, weights=weights2, bins=20, range=(0, 5),label=r'$\rm Region\ II$', normed=False, histtype='step', linewidth=1.6,linestyle='--',stacked=True,fill=False, color='k')
# plt.hist(Av3, weights=weights3, bins=20, range=(0, 5), label=r'$\rm Region\ III$', normed=False, histtype='step', linewidth=1.6,stacked=True,fill=False, color='r')
# plt.legend(loc=1, frameon=False)
# plt.xlim(-0.5,5)
# plt.xlabel(r'$ \rm Dust\ Extinction:\ \itA_{v}$', fontsize=13, labelpad=5)
# plt.savefig('/Users/amandaquirk/Desktop/RGB_av.pdf')
# plt.close()

# def rotation_curve(r, star_vrot, HI_vrot, age, region):
# 	if age=='MS':
# 		color='b'
# 	if age=='RG':
# 		color='r'
# 	if age!='MS' and age !='RG':
# 		color='m'

# 	plt.scatter(r, star_vrot, s=1, c='{}'.format(color))
# 	plt.scatter(r, HI_vrot, s=1, c='k')
# 	plt.xlim(4,20)
# 	plt.ylim(100,300)
# 	plt.xlabel('r [kpc]')
# 	plt.ylabel('rotational v [km/s]')
# 	plt.savefig('/Users/amandaquirk/Desktop/{}_vrot_{}.png'.format(age, region))
# 	plt.close()

# rotation_curve(r1, vrot1, HI_vrot1, 'RG', 1)
# rotation_curve(r2, vrot2, HI_vrot2, 'RG', 2)
# rotation_curve(r3, vrot3, HI_vrot3, 'RG', 3)

#divde region 1 into two sub regions, according to the double peaks in the Av histogram
vrot1_dust=[]
vrot1_clear=[]
r1_dust=[]
r1_clear=[]
HI_vrot1_dust=[]
HI_vrot1_clear=[]
for i in range(len(Av1)):
	if Av1[i]<1:
		vrot1_clear.append(vrot1[i])
		r1_clear.append(r1[i])
		HI_vrot1_clear.append(HI_vrot1[i])
	else:
		vrot1_dust.append(vrot1[i])
		r1_dust.append(r1[i])
		HI_vrot1_dust.append(HI_vrot1[i])

# plt.scatter(r1_clear, vrot1_clear, s=1, c='r', label='Av<1')
# plt.scatter(r1_clear, HI_vrot1_clear, s=1, c='k')
# plt.scatter(r1_dust, vrot1_dust, s=12, c='r', marker='*',label='Av>1')
# plt.scatter(r1_dust, HI_vrot1_dust, s=12, c='k', marker='*')
# plt.xlim(4,20)
# plt.ylim(100,300)
# plt.legend(loc=3)
# plt.xlabel('r [kpc]')
# plt.ylabel('rotational v [km/s]')
# plt.savefig('/Users/amandaquirk/Desktop/RG_vrot_1.png')
# plt.close()

from matplotlib.ticker import MaxNLocator

rc('font', family = 'serif')
f, axes= plt.subplots(3,1, sharey=False, sharex=True, figsize=(4,9.8))
axes[0].scatter(r1_clear, vrot1_clear, s=2, c='r', alpha=0.4, label=r'$A_{v} \rm < 1$')
axes[0].scatter(r1_dust, vrot1_dust, s=8, c='r', marker='s', label=r'$A_{v} \rm \geq 1$')
axes[0].scatter(r1_clear, HI_vrot1_clear, s=2, c='darkgray')
axes[0].scatter(r1_dust, HI_vrot1_dust, s=8, marker='s', c='darkgray')
axes[1].scatter(r2, HI_vrot2, s=2, c='darkgray')
axes[1].scatter(r2, vrot2, s=2, c='r', alpha=0.4)
axes[2].scatter(r3, HI_vrot3, s=2, c='darkgray')
axes[2].scatter(r3, vrot3, s=2, c='r', alpha=0.4)
axes[0].annotate(r'$ \rm Region\ I$', xy=(19,115), horizontalalignment='right', fontsize=12)
axes[1].annotate(r'$\rm Region\ II$', xy=(19,115), horizontalalignment='right', fontsize=12)
axes[2].annotate(r'$\rm Region\ III$', xy=(19,115), horizontalalignment='right', fontsize=12)

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
axes[2].set_xlabel(r'$\rm Radial\ Distance:\ \itr \ \rm(kpc)$', fontsize=13)
nbins = len(axes[0].get_yticklabels())-1
axes[1].yaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='upper'))
axes[2].yaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='upper'))
f.subplots_adjust(left=0.17)
f.text(0.008, 0.5, r'$\rm Rotation\ Velocity:\ \itv_{\rm rot}\ \rm(km\ s^{-1})$', va='center', rotation='vertical', fontsize=13)
plt.subplots_adjust(wspace=0, hspace=0)
axes[0].legend(loc=3, frameon=False)
plt.savefig('/Users/amandaquirk/Desktop/RG_regions_rotation_curves.pdf', bbox_inches='tight')


def asymmetric_drift(v_star, v_gas):
	return [a-b for a,b in zip(v_gas,v_star)]

ad1=asymmetric_drift(vrot1, HI_vrot1)
ad2=asymmetric_drift(vrot2, HI_vrot2)
ad3=asymmetric_drift(vrot3, HI_vrot3)

# rc('font', family = 'serif')
# fig, ax=plt.subplots(1)
# for axis in ['top','bottom','left','right']:
#         ax.spines[axis].set_linewidth(2)
# ax.tick_params(axis='x',which='both',bottom='on',top='off', direction='out')
# ax.tick_params(axis='x',which='both',top='on', direction='in')
# ax.tick_params(axis='y',which='both',left='on',top='off', direction='out')
# ax.tick_params(axis='y',which='both',right='on', direction='in')
# plt.tick_params(which='both', width=2)
# plt.tick_params(which='major', length=7)
# plt.tick_params(which='minor', length=4)
# plt.tick_params(labelsize=12) 
# plt.minorticks_on()

# weights1 = np.ones_like(Av1)/float(len(Av1))
# weights2 = np.ones_like(Av2)/float(len(Av2))
# weights3 = np.ones_like(Av3)/float(len(Av3))
# plt.hist(ad1, weights=weights1, bins=range(-100, 200, 15),label=r'$\rm Region\ I$',normed=False, histtype='step', linewidth=1.6,stacked=True,fill=False, color='b')
# plt.hist(ad2, weights=weights2, bins=range(-100, 200, 15),label=r'$\rm Region\ II$',normed=False, histtype='step', linewidth=1.6,linestyle='--',stacked=True,fill=False, color='k')
# plt.hist(ad3, weights=weights3, bins=range(-100, 200, 15), label=r'$\rm Region\ III$',normed=False, histtype='step', linewidth=1.6,stacked=True,fill=False, color='r')
# plt.legend(loc=1, frameon=False)
# plt.xlim(-100,200)
# plt.xlabel(r'$\rm Asymmetric\ Drift:\ \it v_{a}\ \rm(km\ s^{-1})$', fontsize=13)
# plt.savefig('/Users/amandaquirk/Desktop/RGB__regions_ad.pdf')
# plt.close()

