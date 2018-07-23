import matplotlib.pyplot as plt
import numpy as np 
from matplotlib import rc
from matplotlib.ticker import MaxNLocator

#import data
MS_dispersion, MS_n=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/MS_smoothed_chemin.txt', usecols=(4,5,), unpack=True)
MS_lag=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/MS_ad.txt')
AGy_dispersion, AGy_n=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/AGy_smoothed_chemin.txt', usecols=(4,5,), unpack=True)
AGy_lag=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/AGy_ad.txt')
AGo_dispersion, AGo_n=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/AGo_smoothed_chemin.txt', usecols=(4,5,), unpack=True)
AGo_lag=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/AGo_ad.txt')
RG_dispersion, RG_n=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/RG_smoothed_chemin.txt', usecols=(4,5,), unpack=True)
RG_lag=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/RG_ad.txt')

#sorting the data into bins based on number of componenst of HI
MS_dispersion_1=[a for a,b in zip(MS_dispersion, MS_n) if b==1]
MS_dispersion_2=[a for a,b in zip(MS_dispersion, MS_n) if b==2]
MS_dispersion_3=[a for a,b in zip(MS_dispersion, MS_n) if b==3]
MS_dispersion_4=[a for a,b in zip(MS_dispersion, MS_n) if b==4]
MS_dispersion_5=[a for a,b in zip(MS_dispersion, MS_n) if b==5]
AGy_dispersion_1=[a for a,b in zip(AGy_dispersion, AGy_n) if b==1]
AGy_dispersion_2=[a for a,b in zip(AGy_dispersion, AGy_n) if b==2]
AGy_dispersion_3=[a for a,b in zip(AGy_dispersion, AGy_n) if b==3]
AGy_dispersion_4=[a for a,b in zip(AGy_dispersion, AGy_n) if b==4]
AGy_dispersion_5=[a for a,b in zip(AGy_dispersion, AGy_n) if b==5]
AGo_dispersion_1=[a for a,b in zip(AGo_dispersion, AGo_n) if b==1]
AGo_dispersion_2=[a for a,b in zip(AGo_dispersion, AGo_n) if b==2]
AGo_dispersion_3=[a for a,b in zip(AGo_dispersion, AGo_n) if b==3]
AGo_dispersion_4=[a for a,b in zip(AGo_dispersion, AGo_n) if b==4]
AGo_dispersion_5=[a for a,b in zip(AGo_dispersion, AGo_n) if b==5]
RG_dispersion_1=[a for a,b in zip(RG_dispersion, RG_n) if b==1]
RG_dispersion_2=[a for a,b in zip(RG_dispersion, RG_n) if b==2]
RG_dispersion_3=[a for a,b in zip(RG_dispersion, RG_n) if b==3]
RG_dispersion_4=[a for a,b in zip(RG_dispersion, RG_n) if b==4]
RG_dispersion_5=[a for a,b in zip(RG_dispersion, RG_n) if b==5]

MS_lag_1=[a for a,b in zip(MS_lag, MS_n) if b==1]
MS_lag_2=[a for a,b in zip(MS_lag, MS_n) if b==2]
MS_lag_3=[a for a,b in zip(MS_lag, MS_n) if b==3]
MS_lag_4=[a for a,b in zip(MS_lag, MS_n) if b==4]
MS_lag_5=[a for a,b in zip(MS_lag, MS_n) if b==5]

AGy_lag_1=[a for a,b in zip(AGy_lag, AGy_n) if b==1]
AGy_lag_2=[a for a,b in zip(AGy_lag, AGy_n) if b==2]
AGy_lag_3=[a for a,b in zip(AGy_lag, AGy_n) if b==3]
AGy_lag_4=[a for a,b in zip(AGy_lag, AGy_n) if b==4]
AGy_lag_5=[a for a,b in zip(AGy_lag, AGy_n) if b==5]

AGo_lag_1=[a for a,b in zip(AGo_lag, AGo_n) if b==1]
AGo_lag_2=[a for a,b in zip(AGo_lag, AGo_n) if b==2]
AGo_lag_3=[a for a,b in zip(AGo_lag, AGo_n) if b==3]
AGo_lag_4=[a for a,b in zip(AGo_lag, AGo_n) if b==4]
AGo_lag_5=[a for a,b in zip(AGo_lag, AGo_n) if b==5]

RG_lag_1=[a for a,b in zip(RG_lag, RG_n) if b==1]
RG_lag_2=[a for a,b in zip(RG_lag, RG_n) if b==2]
RG_lag_3=[a for a,b in zip(RG_lag, RG_n) if b==3]
RG_lag_4=[a for a,b in zip(RG_lag, RG_n) if b==4]
RG_lag_5=[a for a,b in zip(RG_lag, RG_n) if b==5]

dispersion_max=max(max(MS_dispersion), max(AGy_dispersion), max(AGo_dispersion), max(RG_dispersion))
dispersion_min=min(min(MS_dispersion), min(AGy_dispersion), min(AGo_dispersion), min(RG_dispersion))
lag_max=max(max(MS_lag), max(AGy_lag), max(AGo_lag), max(RG_lag))
lag_min=min(min(MS_lag), min(AGy_lag), min(AGo_lag), min(RG_lag))

dis_MS_bins=np.linspace(min(MS_dispersion), max(MS_dispersion)+0.5, len(MS_n))
dis_AGy_bins=np.linspace(min(AGy_dispersion), max(AGy_dispersion)+0.5, len(AGy_n))
dis_AGo_bins=np.linspace(min(AGo_dispersion), max(AGo_dispersion)+0.5, len(AGo_n))
dis_RG_bins=np.linspace(min(RG_dispersion), max(RG_dispersion)+0.5, len(RG_n))

lag_MS_bins=np.linspace(min(MS_lag), max(MS_lag)+0.5, len(MS_n))
lag_AGy_bins=np.linspace(min(AGy_lag), max(AGy_lag)+0.5, len(AGy_n))
lag_AGo_bins=np.linspace(min(AGo_lag), max(AGo_lag)+0.5, len(AGo_n))
lag_RG_bins=np.linspace(min(RG_lag), max(RG_lag)+0.5, len(RG_n))

MS_dispersion_y_1=np.zeros(len(dis_MS_bins))
MS_dispersion_y_2=np.zeros(len(dis_MS_bins))
MS_dispersion_y_3=np.zeros(len(dis_MS_bins))
MS_dispersion_y_4=np.zeros(len(dis_MS_bins))
MS_dispersion_y_5=np.zeros(len(dis_MS_bins))

AGy_dispersion_y_1=np.zeros(len(dis_AGy_bins))
AGy_dispersion_y_2=np.zeros(len(dis_AGy_bins))
AGy_dispersion_y_3=np.zeros(len(dis_AGy_bins))
AGy_dispersion_y_4=np.zeros(len(dis_AGy_bins))
AGy_dispersion_y_5=np.zeros(len(dis_AGy_bins))

AGo_dispersion_y_1=np.zeros(len(dis_AGo_bins))
AGo_dispersion_y_2=np.zeros(len(dis_AGo_bins))
AGo_dispersion_y_3=np.zeros(len(dis_AGo_bins))
AGo_dispersion_y_4=np.zeros(len(dis_AGo_bins))
AGo_dispersion_y_5=np.zeros(len(dis_AGo_bins))

RG_dispersion_y_1=np.zeros(len(dis_RG_bins))
RG_dispersion_y_2=np.zeros(len(dis_RG_bins))
RG_dispersion_y_3=np.zeros(len(dis_RG_bins))
RG_dispersion_y_4=np.zeros(len(dis_RG_bins))
RG_dispersion_y_5=np.zeros(len(dis_RG_bins))
for i in range(len(dis_MS_bins)):
	MS_dispersion_1_data=[a for a in MS_dispersion_1 if a<dis_MS_bins[i]]
	MS_dispersion_2_data=[a for a in MS_dispersion_2 if a<dis_MS_bins[i]]
	MS_dispersion_3_data=[a for a in MS_dispersion_3 if a<dis_MS_bins[i]]
	MS_dispersion_4_data=[a for a in MS_dispersion_4 if a<dis_MS_bins[i]]
	MS_dispersion_5_data=[a for a in MS_dispersion_5 if a<dis_MS_bins[i]]
	MS_dispersion_y_1[i]=len(MS_dispersion_1_data)/len(MS_dispersion_1)
	MS_dispersion_y_2[i]=len(MS_dispersion_2_data)/len(MS_dispersion_2)
	MS_dispersion_y_3[i]=len(MS_dispersion_3_data)/len(MS_dispersion_3)
	MS_dispersion_y_4[i]=len(MS_dispersion_4_data)/len(MS_dispersion_4)
	MS_dispersion_y_5[i]=len(MS_dispersion_5_data)/len(MS_dispersion_5)

for i in range(len(dis_AGy_bins)):
	AGy_dispersion_1_data=[a for a in AGy_dispersion_1 if a<dis_AGy_bins[i]]
	AGy_dispersion_2_data=[a for a in AGy_dispersion_2 if a<dis_AGy_bins[i]]
	AGy_dispersion_3_data=[a for a in AGy_dispersion_3 if a<dis_AGy_bins[i]]
	AGy_dispersion_4_data=[a for a in AGy_dispersion_4 if a<dis_AGy_bins[i]]
	AGy_dispersion_5_data=[a for a in AGy_dispersion_5 if a<dis_AGy_bins[i]]
	AGy_dispersion_y_1[i]=len(AGy_dispersion_1_data)/len(AGy_dispersion_1)
	AGy_dispersion_y_2[i]=len(AGy_dispersion_2_data)/len(AGy_dispersion_2)
	AGy_dispersion_y_3[i]=len(AGy_dispersion_3_data)/len(AGy_dispersion_3)
	AGy_dispersion_y_4[i]=len(AGy_dispersion_4_data)/len(AGy_dispersion_4)
	AGy_dispersion_y_5[i]=len(AGy_dispersion_5_data)/len(AGy_dispersion_5)

for i in range(len(dis_AGo_bins)):
	AGo_dispersion_1_data=[a for a in AGo_dispersion_1 if a<dis_AGo_bins[i]]
	AGo_dispersion_2_data=[a for a in AGo_dispersion_2 if a<dis_AGo_bins[i]]
	AGo_dispersion_3_data=[a for a in AGo_dispersion_3 if a<dis_AGo_bins[i]]
	AGo_dispersion_4_data=[a for a in AGo_dispersion_4 if a<dis_AGo_bins[i]]
	AGo_dispersion_5_data=[a for a in AGo_dispersion_5 if a<dis_AGo_bins[i]]
	AGo_dispersion_y_1[i]=len(AGo_dispersion_1_data)/len(AGo_dispersion_1)
	AGo_dispersion_y_2[i]=len(AGo_dispersion_2_data)/len(AGo_dispersion_2)
	AGo_dispersion_y_3[i]=len(AGo_dispersion_3_data)/len(AGo_dispersion_3)
	AGo_dispersion_y_4[i]=len(AGo_dispersion_4_data)/len(AGo_dispersion_4)
	AGo_dispersion_y_5[i]=len(AGo_dispersion_5_data)/len(AGo_dispersion_5)

for i in range(len(dis_RG_bins)):
	RG_dispersion_1_data=[a for a in RG_dispersion_1 if a<dis_RG_bins[i]]
	RG_dispersion_2_data=[a for a in RG_dispersion_2 if a<dis_RG_bins[i]]
	RG_dispersion_3_data=[a for a in RG_dispersion_3 if a<dis_RG_bins[i]]
	RG_dispersion_4_data=[a for a in RG_dispersion_4 if a<dis_RG_bins[i]]
	RG_dispersion_5_data=[a for a in RG_dispersion_5 if a<dis_RG_bins[i]]
	RG_dispersion_y_1[i]=len(RG_dispersion_1_data)/len(RG_dispersion_1)
	RG_dispersion_y_2[i]=len(RG_dispersion_2_data)/len(RG_dispersion_2)
	RG_dispersion_y_3[i]=len(RG_dispersion_3_data)/len(RG_dispersion_3)
	RG_dispersion_y_4[i]=len(RG_dispersion_4_data)/len(RG_dispersion_4)
	RG_dispersion_y_5[i]=len(RG_dispersion_5_data)/len(RG_dispersion_5)

MS_lag_y_1=np.zeros(len(lag_MS_bins))
MS_lag_y_2=np.zeros(len(lag_MS_bins))
MS_lag_y_3=np.zeros(len(lag_MS_bins))
MS_lag_y_4=np.zeros(len(lag_MS_bins))
MS_lag_y_5=np.zeros(len(lag_MS_bins))

AGy_lag_y_1=np.zeros(len(lag_AGy_bins))
AGy_lag_y_2=np.zeros(len(lag_AGy_bins))
AGy_lag_y_3=np.zeros(len(lag_AGy_bins))
AGy_lag_y_4=np.zeros(len(lag_AGy_bins))
AGy_lag_y_5=np.zeros(len(lag_AGy_bins))

AGo_lag_y_1=np.zeros(len(lag_AGo_bins))
AGo_lag_y_2=np.zeros(len(lag_AGo_bins))
AGo_lag_y_3=np.zeros(len(lag_AGo_bins))
AGo_lag_y_4=np.zeros(len(lag_AGo_bins))
AGo_lag_y_5=np.zeros(len(lag_AGo_bins))

RG_lag_y_1=np.zeros(len(lag_RG_bins))
RG_lag_y_2=np.zeros(len(lag_RG_bins))
RG_lag_y_3=np.zeros(len(lag_RG_bins))
RG_lag_y_4=np.zeros(len(lag_RG_bins))
RG_lag_y_5=np.zeros(len(lag_RG_bins))
for i in range(len(lag_MS_bins)):
	MS_lag_1_data=[a for a in MS_lag_1 if a<lag_MS_bins[i]]
	MS_lag_2_data=[a for a in MS_lag_2 if a<lag_MS_bins[i]]
	MS_lag_3_data=[a for a in MS_lag_3 if a<lag_MS_bins[i]]
	MS_lag_4_data=[a for a in MS_lag_4 if a<lag_MS_bins[i]]
	MS_lag_5_data=[a for a in MS_lag_5 if a<lag_MS_bins[i]]
	MS_lag_y_1[i]=len(MS_lag_1_data)/len(MS_lag_1)
	MS_lag_y_2[i]=len(MS_lag_2_data)/len(MS_lag_2)
	MS_lag_y_3[i]=len(MS_lag_3_data)/len(MS_lag_3)
	MS_lag_y_4[i]=len(MS_lag_4_data)/len(MS_lag_4)
	MS_lag_y_5[i]=len(MS_lag_5_data)/len(MS_lag_5)

for i in range(len(lag_AGy_bins)):
	AGy_lag_1_data=[a for a in AGy_lag_1 if a<lag_AGy_bins[i]]
	AGy_lag_2_data=[a for a in AGy_lag_2 if a<lag_AGy_bins[i]]
	AGy_lag_3_data=[a for a in AGy_lag_3 if a<lag_AGy_bins[i]]
	AGy_lag_4_data=[a for a in AGy_lag_4 if a<lag_AGy_bins[i]]
	AGy_lag_5_data=[a for a in AGy_lag_5 if a<lag_AGy_bins[i]]
	AGy_lag_y_1[i]=len(AGy_lag_1_data)/len(AGy_lag_1)
	AGy_lag_y_2[i]=len(AGy_lag_2_data)/len(AGy_lag_2)
	AGy_lag_y_3[i]=len(AGy_lag_3_data)/len(AGy_lag_3)
	AGy_lag_y_4[i]=len(AGy_lag_4_data)/len(AGy_lag_4)
	AGy_lag_y_5[i]=len(AGy_lag_5_data)/len(AGy_lag_5)

for i in range(len(lag_AGo_bins)):
	AGo_lag_1_data=[a for a in AGo_lag_1 if a<lag_AGo_bins[i]]
	AGo_lag_2_data=[a for a in AGo_lag_2 if a<lag_AGo_bins[i]]
	AGo_lag_3_data=[a for a in AGo_lag_3 if a<lag_AGo_bins[i]]
	AGo_lag_4_data=[a for a in AGo_lag_4 if a<lag_AGo_bins[i]]
	AGo_lag_5_data=[a for a in AGo_lag_5 if a<lag_AGo_bins[i]]
	AGo_lag_y_1[i]=len(AGo_lag_1_data)/len(AGo_lag_1)
	AGo_lag_y_2[i]=len(AGo_lag_2_data)/len(AGo_lag_2)
	AGo_lag_y_3[i]=len(AGo_lag_3_data)/len(AGo_lag_3)
	AGo_lag_y_4[i]=len(AGo_lag_4_data)/len(AGo_lag_4)
	AGo_lag_y_5[i]=len(AGo_lag_5_data)/len(AGo_lag_5)

for i in range(len(lag_RG_bins)):
	RG_lag_1_data=[a for a in RG_lag_1 if a<lag_RG_bins[i]]
	RG_lag_2_data=[a for a in RG_lag_2 if a<lag_RG_bins[i]]
	RG_lag_3_data=[a for a in RG_lag_3 if a<lag_RG_bins[i]]
	RG_lag_4_data=[a for a in RG_lag_4 if a<lag_RG_bins[i]]
	RG_lag_5_data=[a for a in RG_lag_5 if a<lag_RG_bins[i]]
	RG_lag_y_1[i]=len(RG_lag_1_data)/len(RG_lag_1)
	RG_lag_y_2[i]=len(RG_lag_2_data)/len(RG_lag_2)
	RG_lag_y_3[i]=len(RG_lag_3_data)/len(RG_lag_3)
	RG_lag_y_4[i]=len(RG_lag_4_data)/len(RG_lag_4)
	RG_lag_y_5[i]=len(RG_lag_5_data)/len(RG_lag_5)

rc('font', family = 'serif')
f, axes= plt.subplots(1,4, sharey=True, sharex=False, figsize=(10,3.5))

axes[0].plot([0,0], [-0.2,1.2], c='darkgrey', linestyle='--')
#axes[0].plot([75,75], [-0.2,1.2], c='darkgrey', linestyle='--')
axes[0].plot(lag_MS_bins, MS_lag_y_1, c='mediumpurple', label='n=1')
axes[0].plot(lag_MS_bins, MS_lag_y_2, c='deepskyblue', label='n=2')
axes[0].plot(lag_MS_bins, MS_lag_y_3, c='lightgreen',label='n=3')
axes[0].plot(lag_MS_bins, MS_lag_y_4, c='orange', label='n=4')
axes[0].plot(lag_MS_bins, MS_lag_y_5, c='r', label='n=5')

axes[1].plot([0,0], [-0.2,1.2], c='darkgrey', linestyle='--')
#axes[1].plot([75,75], [-0.2,1.2], c='darkgrey', linestyle='--')
axes[1].plot(lag_AGy_bins, AGy_lag_y_1, c='mediumpurple', label='n=1')
axes[1].plot(lag_AGy_bins, AGy_lag_y_2, c='deepskyblue', label='n=2')
axes[1].plot(lag_AGy_bins, AGy_lag_y_3, c='lightgreen',label='n=3')
axes[1].plot(lag_AGy_bins, AGy_lag_y_4, c='orange', label='n=4')
axes[1].plot(lag_AGy_bins, AGy_lag_y_5, c='r', label='n=5')

axes[2].plot([0,0], [-0.2,1.2], c='darkgrey', linestyle='--')
#axes[2].plot([75,75], [-0.2,1.2], c='darkgrey', linestyle='--')
axes[2].plot(lag_AGo_bins, AGo_lag_y_1, c='mediumpurple', label='n=1')
axes[2].plot(lag_AGo_bins, AGo_lag_y_2, c='deepskyblue', label='n=2')
axes[2].plot(lag_AGo_bins, AGo_lag_y_3, c='lightgreen',label='n=3')
axes[2].plot(lag_AGo_bins, AGo_lag_y_4, c='orange', label='n=4')
axes[2].plot(lag_AGo_bins, AGo_lag_y_5, c='r', label='n=5')

axes[3].plot([0,0], [-0.2,1.2], c='darkgrey', linestyle='--')
#axes[3].plot([75,75], [-0.2,1.2], c='darkgrey', linestyle='--')
axes[3].plot(lag_RG_bins, RG_lag_y_1, c='mediumpurple', label='n=1')
axes[3].plot(lag_RG_bins, RG_lag_y_2, c='deepskyblue', label='n=2')
axes[3].plot(lag_RG_bins, RG_lag_y_3, c='lightgreen',label='n=3')
axes[3].plot(lag_RG_bins, RG_lag_y_4, c='orange', label='n=4')
axes[3].plot(lag_RG_bins, RG_lag_y_5, c='r', label='n=5')

for ax in axes:
	#ax.set_xlim(10, 160)
	#ax.set_xlabel(r'$\rm \sigma\ (km\ s^{-1})$', fontsize=13)
	ax.set_xlim(-100, 200)
	#ax.set_xlabel(r'v_{a} \ \rm(km\ s^{-1})', fontsize=13)
	ax.set_ylim(-0.1,1.1)
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
#axes[0].set_ylabel(r'$\rm N(< \sigma)$', fontsize=13)
axes[0].set_ylabel(r'$\rm N(< \it v_{a}\rm)$', fontsize=13)
axes[0].annotate('MS', xy=(177,0.7), horizontalalignment='right', fontsize=12) #177 for va, 150 for sigma
axes[1].annotate('young AGB', xy=(177,0), horizontalalignment='right', fontsize=12)
axes[2].annotate('older AGB', xy=(177,0), horizontalalignment='right', fontsize=12)
axes[3].annotate('RGB', xy=(177,0), horizontalalignment='right', fontsize=12)
nbins = len(axes[0].get_xticklabels())
axes[0].xaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='upper'))
axes[1].xaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='upper'))
axes[2].xaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='upper'))
axes[3].xaxis.set_major_locator(MaxNLocator(nbins=nbins))
plt.subplots_adjust(wspace=0, hspace=0)
axes[0].legend(frameon=False)
f.subplots_adjust(bottom=0.13)
f.text(0.5, 0.008, r'$\rmAsymmetric\ Drift:\ \it v_{a}\ \rm(km\ s^{-1})$', ha='center', fontsize=13)
plt.savefig('/Users/amandaquirk/Desktop/lag_hist.pdf', bbox_inches='tight')
plt.close()




















