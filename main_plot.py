import matplotlib.pyplot as plt 
from matplotlib import rc 
import numpy as np 
from matplotlib.ticker import MaxNLocator
from matplotlib import patches

MS_avg_xi, MS_avg_eta,  MS_avg_v, MS_var=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/MS_smoothed_chemin.txt', usecols=(0,1,2,4,), unpack=True)
AGy_avg_xi, AGy_avg_eta, AGy_avg_v, AGy_var=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/AGy_smoothed_chemin.txt', usecols=(0,1,2,4,), unpack=True)
AGo_avg_xi, AGo_avg_eta, AGo_avg_v, AGo_var=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/AGo_smoothed_chemin.txt', usecols=(0,1,2,4,), unpack=True)
RG_avg_xi, RG_avg_eta, RG_avg_v, RG_var=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/RG_smoothed_chemin.txt', usecols=(0,1,2,4,), unpack=True)

MS_xi, MS_eta,  MS_v=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/MS_individual_chemin.txt', usecols=(0,1,2,), unpack=True)
AGy_xi, AGy_eta, AGy_v=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/AGy_individual_chemin.txt', usecols=(0,1,2,), unpack=True)
AGo_xi, AGo_eta, AGo_v=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/AGo_individual_chemin.txt', usecols=(0,1,2,), unpack=True)
RG_xi, RG_eta, RG_v=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/RG_individual_chemin.txt', usecols=(0,1,2,), unpack=True)

MS_xi=[13.67*a for a in MS_xi]
MS_eta=[13.67*a for a in MS_eta]
AGy_xi=[13.67*a for a in AGy_xi]
AGy_eta=[13.67*a for a in AGy_eta]
AGo_xi=[13.67*a for a in AGo_xi]
AGo_eta=[13.67*a for a in AGo_eta]
RG_xi=[13.67*a for a in RG_xi]
RG_eta=[13.67*a for a in RG_eta]

ylength=37.5
xlength=ylength*np.cos(float(73.7*np.pi) / 180)
e =  patches.Ellipse((0, 0), xlength, ylength, angle=180-37.1, linewidth=.5, fill=False, zorder=0)
e1 = patches.Ellipse((0, 0), xlength, ylength, angle=180-37.1, linewidth=.5, fill=False, zorder=0)
e2 = patches.Ellipse((0, 0), xlength, ylength, angle=180-37.1, linewidth=.5, fill=False, zorder=0)
e3 = patches.Ellipse((0, 0), xlength, ylength, angle=180-37.1, linewidth=.5, fill=False, zorder=0)
e4 = patches.Ellipse((0, 0), xlength, ylength, angle=180-37.1, linewidth=.5, fill=False, zorder=0)
e5 =  patches.Ellipse((0, 0), xlength, ylength, angle=180-37.1, linewidth=.5, fill=False, zorder=0)
e6 = patches.Ellipse((0, 0), xlength, ylength, angle=180-37.1, linewidth=.5, fill=False, zorder=0)
e7 = patches.Ellipse((0, 0), xlength, ylength, angle=180-37.1, linewidth=.5, fill=False, zorder=0)
e8 = patches.Ellipse((0, 0), xlength, ylength, angle=180-37.1, linewidth=.5, fill=False, zorder=0)
e9 = patches.Ellipse((0, 0), xlength, ylength, angle=180-37.1, linewidth=.5, fill=False, zorder=0)
e10 =  patches.Ellipse((0, 0), xlength, ylength, angle=180-37.1, linewidth=.5, fill=False, zorder=0)
e11 = patches.Ellipse((0, 0), xlength, ylength, angle=180-37.1, linewidth=.5, fill=False, zorder=0)

rc('font', family = 'serif')
f, axes= plt.subplots(3,4, sharey=True, sharex=False, figsize=(11.1,10.26))
axes[0,0].add_patch(e)
axes[0,0].scatter(MS_xi,MS_eta, c=MS_v, cmap='plasma', s=4, vmin=-300,vmax=0) 
axes[0,1].add_patch(e1)
axes[0,1].scatter(AGy_xi,AGy_eta, c=AGy_v, cmap='plasma', s=4, vmin=-300,vmax=0) 
axes[0,2].add_patch(e2)
axes[0,2].scatter(AGo_xi,AGo_eta, c=AGo_v, cmap='plasma', s=4, vmin=-300,vmax=0) 
axes[0,3].add_patch(e3)
im0=axes[0,3].scatter(RG_xi,RG_eta, c=RG_v, cmap='plasma', s=4, vmin=-300,vmax=0)

axes[0,0].annotate('MS', xy=(-.55,13.5), horizontalalignment='right', fontsize=12)
axes[0,1].annotate('young AGB', xy=(-.550,13.5), horizontalalignment='right', fontsize=12)
axes[0,2].annotate('older AGB', xy=(-.550,13.5), horizontalalignment='right', fontsize=12)
axes[0,3].annotate('RGB', xy=(-.550,13.5), horizontalalignment='right', fontsize=12)

axes[1,0].add_patch(e4)
axes[1,0].scatter(MS_avg_xi,MS_avg_eta, c=MS_avg_v, cmap='plasma', s=4, vmin=-300,vmax=0) 
axes[1,1].add_patch(e5)
axes[1,1].scatter(AGy_avg_xi,AGy_avg_eta, c=AGy_avg_v, cmap='plasma', s=4, vmin=-300,vmax=0)
axes[1,2].add_patch(e6) 
axes[1,2].scatter(AGo_avg_xi,AGo_avg_eta, c=AGo_avg_v, cmap='plasma', s=4, vmin=-300,vmax=0) 
axes[1,3].add_patch(e7)
im1=axes[1,3].scatter(RG_avg_xi,RG_avg_eta, c=RG_avg_v, cmap='plasma', s=4, vmin=-300,vmax=0)

axes[1,0].annotate('MS', xy=(-.550,13.5), horizontalalignment='right', fontsize=12)
axes[1,1].annotate('young AGB', xy=(-.550,13.5), horizontalalignment='right', fontsize=12)
axes[1,2].annotate('older AGB', xy=(-.550,13.5), horizontalalignment='right', fontsize=12)
axes[1,3].annotate('RGB', xy=(-.550,13.5), horizontalalignment='right', fontsize=12)

axes[2,0].add_patch(e8)
axes[2,0].scatter(MS_avg_xi,MS_avg_eta, c=MS_var, cmap='copper', s=4, vmin=10,vmax=150) 
axes[2,1].add_patch(e9)
axes[2,1].scatter(AGy_avg_xi,AGy_avg_eta, c=AGy_var, cmap='copper', s=4, vmin=10,vmax=150)
axes[2,2].add_patch(e10)
axes[2,2].scatter(AGo_avg_xi,AGo_avg_eta, c=AGo_var, cmap='copper', s=4,vmin=10,vmax=150) 
axes[2,3].add_patch(e11)
im2=axes[2,3].scatter(RG_avg_xi,RG_avg_eta, c=RG_var, cmap='copper', s=4,vmin=10,vmax=150)

axes[2,0].annotate('MS', xy=(-.550,13.5), horizontalalignment='right', fontsize=12)
axes[2,1].annotate('young AGB', xy=(-.550,13.5), horizontalalignment='right', fontsize=12)
axes[2,2].annotate('older AGB', xy=(-.550,13.5), horizontalalignment='right', fontsize=12)
axes[2,3].annotate('RGB', xy=(-.550,13.5), horizontalalignment='right', fontsize=12)

for ax in axes[0,:]:
	ax.set_xlim(13, -1.5)
	ax.tick_params(axis='x',which='both',bottom='on',top='off', direction='out')
	ax.tick_params(axis='x',which='both',top='on', direction='in')
	ax.tick_params(which='both', width=1)
	ax.tick_params(which='major', length=7)
	ax.tick_params(which='minor', length=4)
	ax.tick_params(labelsize=12) 
	ax.minorticks_on()
	for axis in ['top','bottom','left','right']:
	        ax.spines[axis].set_linewidth(1)

for ax in axes[1,:]:
	ax.set_xlim(13, -1.5)
	ax.tick_params(axis='x',which='both',bottom='on',top='off', direction='out')
	ax.tick_params(axis='x',which='both',top='on', direction='in')
	ax.tick_params(which='both', width=1)
	ax.tick_params(which='major', length=7)
	ax.tick_params(which='minor', length=4)
	ax.tick_params(labelsize=12) 
	ax.minorticks_on()
	for axis in ['top','bottom','left','right']:
	        ax.spines[axis].set_linewidth(1)

for ax in axes[2,:]:
	ax.set_xlim(13, -1.5)
	ax.set_xlabel(r'$\xi\ (kpc)$', fontsize=13)
	ax.tick_params(axis='x',which='both',bottom='on',top='off', direction='out')
	ax.tick_params(axis='x',which='both',top='on', direction='in')
	ax.tick_params(which='both', width=1)
	ax.tick_params(which='major', length=7)
	ax.tick_params(which='minor', length=4)
	ax.tick_params(labelsize=12) 
	ax.minorticks_on()
	for axis in ['top','bottom','left','right']:
	        ax.spines[axis].set_linewidth(1)

for ax in axes[:,0]:
	ax.set_ylabel(r'$\eta\ (kpc)$', fontsize=13)
	ax.set_ylim(-2.5,15.5)
	ax.tick_params(axis='y',which='both',left='on',top='off', direction='out')
	ax.tick_params(axis='y',which='both',right='on', direction='in')
	ax.tick_params(which='both', width=1)
	ax.tick_params(which='major', length=7)
	ax.tick_params(which='minor', length=4)
	ax.tick_params(labelsize=12) 
	nbins = 7
	ax.yaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='upper'))
	ax.minorticks_on()
	for axis in ['top','bottom','left','right']:
	        ax.spines[axis].set_linewidth(1)

for ax in axes[:,1]:
	ax.set_ylim(-2.5,15.5)
	ax.tick_params(axis='y',which='both',left='on',top='off', direction='out')
	ax.tick_params(axis='y',which='both',right='on', direction='in')
	ax.tick_params(which='both', width=1)
	ax.tick_params(which='major', length=7)
	ax.tick_params(which='minor', length=4)
	ax.tick_params(labelsize=12) 
	ax.minorticks_on()
	for axis in ['top','bottom','left','right']:
	        ax.spines[axis].set_linewidth(1)

for ax in axes[:,2]:
	ax.set_ylim(-2.5,15.5)
	ax.tick_params(axis='y',which='both',left='on',top='off', direction='out')
	ax.tick_params(axis='y',which='both',right='on', direction='in')
	ax.tick_params(which='both', width=1)
	ax.tick_params(which='major', length=7)
	ax.tick_params(which='minor', length=4)
	ax.tick_params(labelsize=12) 
	ax.minorticks_on()
	for axis in ['top','bottom','left','right']:
	        ax.spines[axis].set_linewidth(1)

for ax in axes[:,3]:
	ax.set_ylim(-2.5,15.5)
	ax.tick_params(axis='y',which='both',left='on',top='off', direction='out')
	ax.tick_params(axis='y',which='both',right='on', direction='in')
	ax.tick_params(which='both', width=1)
	ax.tick_params(which='major', length=7)
	ax.tick_params(which='minor', length=4)
	ax.tick_params(labelsize=12) 
	ax.minorticks_on()
	for axis in ['top','bottom','left','right']:
	        ax.spines[axis].set_linewidth(1)

f.subplots_adjust(right=0.885)
cbar_ax1 = f.add_axes([0.89,0.375,0.015,0.512])
cbar_ax2 = f.add_axes([0.89,0.1,0.015,0.27])
clb1=f.colorbar(im1, cax=cbar_ax1)
clb2=f.colorbar(im2, cax=cbar_ax2)
clb1.set_label(r'$\rm Individual,\ Mean\ LOS\ velocity:\ v, \ \overline{v}\ (km\ s^{-1})$', fontsize=13)
clb2.set_label(r'$\rm Velocity\ Dispersion: \sigma\ (km\ s^{-1})$', fontsize=13, labelpad=12)
axes[0,0].scatter(0,0,marker='+', c='dodgerblue', linewidth=2)
axes[0,1].scatter(0,0,marker='+', c='dodgerblue', linewidth=2)
axes[0,2].scatter(0,0,marker='+', c='dodgerblue', linewidth=2)
axes[0,3].scatter(0,0,marker='+', c='dodgerblue', linewidth=2)
axes[1,0].scatter(0,0,marker='+', c='dodgerblue', linewidth=2)
axes[1,1].scatter(0,0,marker='+', c='dodgerblue', linewidth=2)
axes[1,2].scatter(0,0,marker='+', c='dodgerblue', linewidth=2)
axes[1,3].scatter(0,0,marker='+', c='dodgerblue', linewidth=2)
axes[2,0].scatter(0,0,marker='+', c='dodgerblue', linewidth=2)
axes[2,1].scatter(0,0,marker='+', c='dodgerblue', linewidth=2)
axes[2,2].scatter(0,0,marker='+', c='dodgerblue', linewidth=2)
axes[2,3].scatter(0,0,marker='+', c='dodgerblue', linewidth=2)
plt.subplots_adjust(wspace=0, hspace=0)
plt.savefig('/Users/amandaquirk/Desktop/maps.pdf', bbox_inches='tight')
