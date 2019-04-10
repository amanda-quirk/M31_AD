import matplotlib.pyplot as plt 
from matplotlib import rc 
import numpy as np 

#age = [.03, .4, 2, 4] #Gyr
log_age = [np.log10(30e6), np.log10(.4e9), np.log10(2e9), np.log10(4e9)]
CO_AD = [-12.86, 0.18,36.99, 37.15] #km/s
HI_AD = [-8.15, 17.69, 50.43, 62.97] #km/s
HI_low_errors = [23.515 / np.sqrt(1009), 68.237/ np.sqrt(430), 32.248 / np.sqrt(880), 33.043 / np.sqrt(3129)]
HI_high_errors = [22.951 / np.sqrt(1009), 57.714 / np.sqrt(430), 37.272 / np.sqrt(880), 22.351 / np.sqrt(3129)]
CO_low_errors = [35.494 / np.sqrt(738), 117.362 / np.sqrt(338), 55.888 / np.sqrt(713), 38.107 / np.sqrt(2320)]
CO_high_errors = [23.334 / np.sqrt(738), 76.638 / np.sqrt(338), 41.397 / np.sqrt(713), 30.796/ np.sqrt(2320)]

IWM_HI_AD = [-13.54, 6.32, 36.45, 50.43] #km/s
IWM_HI_low_errors = [21.34 / np.sqrt(983), 36.66/ np.sqrt(420), 37.61 / np.sqrt(836), 25.31 / np.sqrt(2998)]
IWM_HI_high_errors = [27.44 / np.sqrt(983), 81.77 / np.sqrt(420), 41.98 / np.sqrt(836), 21.14 / np.sqrt(2998)]

HI_errors = np.row_stack((HI_low_errors, HI_high_errors))
CO_errors = np.row_stack((CO_low_errors, CO_high_errors))
IWM_HI_errors = np.row_stack((IWM_HI_low_errors, IWM_HI_high_errors))
 
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

elb1 = plt.errorbar(log_age, CO_AD, yerr=CO_errors, capsize=7, c='darkcyan', linewidth = 3, label=r'$\rm w.r.t.\ CO$')
elb2 = plt.errorbar(log_age, HI_AD, yerr=HI_errors, capsize=7,  c='darkgrey', linewidth = 3,linestyle='-', label=r'$\rm w.r.t.\ HI$')
#elb3 = plt.errorbar(log_age, IWM_HI_AD, yerr=IWM_HI_errors, capsize=5,  linestyle='--', c='darkgrey', label=r'$\rm w.r.t.\ IWM\ HI$')
# elb1[-1][0].set_linestyle('-.')
#elb3[-1][0].set_linestyle('--')
plt.legend(frameon=False, loc=4)
#plt.xscale('log')
#plt.xlabel(r'$\rm Average\ Stellar\ Age\ (Gyr)$', fontsize=13)
plt.xlabel(r'$\rm Average\ Stellar\ Age\ [log(yr)]$', fontsize=13)
plt.ylabel(r'$\rm Median\ Asymmetric\ Drift:\ \itv_{a}\ \rm(km\ s^{-1})$', fontsize=13)
plt.savefig('/Users/amandaquirk/Desktop/AD_transition.pdf', bbox_inches='tight')
#print('HI IWM')
#print(HI_errors)
#print('CO')
#print(CO_errors)
