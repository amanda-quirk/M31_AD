import matplotlib.pyplot as plt 
import numpy as np 
from matplotlib import rc
from matplotlib import patches

MS_xi, MS_eta, MS_Av=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/MS_dust.txt', usecols=(0,1,2,), unpack=True)
AGy_xi, AGy_eta, AGy_Av=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/AGy_dust.txt', usecols=(0,1,2,), unpack=True)
AGo_xi, AGo_eta, AGo_Av=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/AGo_dust.txt', usecols=(0,1,2,), unpack=True)
RG_xi, RG_eta, RG_Av=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/RG_dust.txt', usecols=(0,1,2,), unpack=True)


#setting bounds for the histogram
Av_max=max(max(MS_Av), max(AGy_Av), max(AGo_Av), max(RG_Av))
Av_min=min(min(MS_Av), min(AGy_Av), min(AGo_Av), min(RG_Av))

MS_bins=np.linspace(min(MS_Av), max(MS_Av)+0.5, len(MS_xi))
AGy_bins=np.linspace(min(AGy_Av), max(AGy_Av)+0.5, len(AGy_xi))
AGo_bins=np.linspace(min(AGo_Av), max(AGo_Av)+0.5, len(AGo_xi))
RG_bins=np.linspace(min(RG_Av), max(RG_Av)+0.5, len(RG_xi))

#sorting data into bins
MS_y=np.zeros(len(MS_bins))
AGy_y=np.zeros(len(AGy_bins))
AGo_y=np.zeros(len(AGo_bins))
RG_y=np.zeros(len(RG_bins))
for i in range(len(MS_bins)):
	MS_data=[a for a in MS_Av if a<MS_bins[i]]
	MS_y[i]=len(MS_data)/len(MS_Av)
	if abs(1-MS_bins[i])<0.001:
		print(MS_y[i])
print("one")	
for i in range(len(AGy_bins)):
	AGy_data=[a for a in AGy_Av if a<AGy_bins[i]]
	AGy_y[i]=len(AGy_data)/len(AGy_Av)
	if abs(1-AGy_bins[i])<0.007:
		print(AGy_y[i])
print("two")
for i in range(len(AGo_bins)):
	AGo_data=[a for a in AGo_Av if a<AGo_bins[i]]
	AGo_y[i]=len(AGo_data)/len(AGo_Av)
	if abs(1-AGo_bins[i])<0.005:
		print(AGo_y[i])
print("three")
for i in range(len(RG_bins)):
	RG_data=[a for a in RG_Av if a<RG_bins[i]]
	RG_y[i]=len(RG_data)/len(RG_Av)
	if abs(1-RG_bins[i])<0.001:
		print(RG_y[i])

#plotting
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

plt.plot(MS_bins, MS_y, c='b', label='MS')
plt.plot(AGy_bins, AGy_y, c='m', label='young AGB')
plt.plot(AGo_bins, AGo_y, c='k', label='old AGB')
plt.plot(RG_bins, RG_y, c='r', label='RGB')
plt.plot([1,1], [-1,2], c='darkgrey', linestyle=':')
plt.xlim(-0.05, 5)
plt.ylim(-0.05, 1.05)
plt.xlabel(r'$\rm Dust\ Extinction:\ \itA_{v}$', fontsize=14)
plt.ylabel(r'$\rm N(< \itA_{v}\rm)$', fontsize=14)
plt.legend(frameon=False, fontsize=10)
plt.savefig('/Users/amandaquirk/Desktop/dust_hists.pdf', bbox='tight')

# rc('font', family = 'serif')
# f, axes= plt.subplots(1,4, sharey=True, sharex=False, figsize=(11.1,3.42))
# ylength=37.5
# xlength=ylength*np.cos(float(73.7*np.pi) / 180)
# e =  patches.Ellipse((0, 0), xlength, ylength, angle=180-37.1, linewidth=.5, fill=False, zorder=0)
# e1 = patches.Ellipse((0, 0), xlength, ylength, angle=180-37.1, linewidth=.5, fill=False, zorder=0)
# e2 = patches.Ellipse((0, 0), xlength, ylength, angle=180-37.1, linewidth=.5, fill=False, zorder=0)
# e3 = patches.Ellipse((0, 0), xlength, ylength, angle=180-37.1, linewidth=.5, fill=False, zorder=0)
# axes[3].add_patch(e3)
# axes[2].add_patch(e2)
# axes[1].add_patch(e1)
# axes[0].add_patch(e)
# axes[0].scatter(MS_xi,MS_eta, c=MS_Av, cmap='gist_heat_r', s=7, vmin=Av_min, vmax=Av_max+0.2) 
# axes[1].scatter(AGy_xi,AGy_eta, c=AGy_Av, cmap='gist_heat_r', s=7, vmin=Av_min, vmax=Av_max+0.2) 
# axes[2].scatter(AGo_xi,AGo_eta, c=AGo_Av, cmap='gist_heat_r', s=7, vmin=Av_min, vmax=Av_max+0.2) 
# im=axes[3].scatter(RG_xi,RG_eta, c=RG_Av, cmap='gist_heat_r', s=7, vmin=Av_min, vmax=Av_max+0.2)
# for ax in axes:
# 	ax.set_xlim(13, -1.5)
# 	ax.set_xlabel(r'$\xi\ (kpc)$', fontsize=13)
# 	ax.set_ylim(-2.5,15.5)
# 	ax.tick_params(axis='x',which='both',bottom='on',top='off', direction='out')
# 	ax.tick_params(axis='x',which='both',top='on', direction='in')
# 	ax.tick_params(axis='y',which='both',left='on',top='off', direction='out')
# 	ax.tick_params(axis='y',which='both',right='on', direction='in')
# 	ax.tick_params(which='both', width=1)
# 	ax.tick_params(which='major', length=7)
# 	ax.tick_params(which='minor', length=4)
# 	ax.tick_params(labelsize=12) 
# 	ax.minorticks_on()
# 	for axis in ['top','bottom','left','right']:
# 	        ax.spines[axis].set_linewidth(1)
# axes[0].set_ylabel(r'$\eta\ (kpc)$', fontsize=13)
# f.subplots_adjust(right=0.885)
# cbar_ax = f.add_axes([0.89, 0.13, 0.015, 0.72])
# clb=f.colorbar(im, cax=cbar_ax)
# clb.set_label(r'$\rm Dust\ Extinction:\ A_{v}$', fontsize=13)
# axes[0].scatter(0,0,marker='+', c='b', linewidth=2)
# axes[1].scatter(0,0,marker='+', c='b', linewidth=2)
# axes[2].scatter(0,0,marker='+', c='b', linewidth=2)
# axes[3].scatter(0,0,marker='+', c='b', linewidth=2)
# axes[0].annotate('MS', xy=(-.55,13.5), horizontalalignment='right', fontsize=12)
# axes[1].annotate('young AGB', xy=(-.55,13.5), horizontalalignment='right', fontsize=12)
# axes[2].annotate('older AGB', xy=(-.55,13.5), horizontalalignment='right', fontsize=12)
# axes[3].annotate('RGB', xy=(-.55,13.5), horizontalalignment='right', fontsize=12)
# plt.subplots_adjust(wspace=0, hspace=0)
# plt.savefig('/Users/amandaquirk/Desktop/dust_maps.pdf', bbox_inches='tight')



