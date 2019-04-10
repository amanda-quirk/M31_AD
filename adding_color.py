# from astropy.io import fits 
# from astropy.table import Table 
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt 
import numpy as np 
from matplotlib import rc,rcParams
from shapely.geometry import Point 
from shapely.geometry.polygon import Polygon

#import the data file
# hdulist= fits.open('/Users/amandaquirk/Documents/AsymmetricDrift/Data/subMasterSPLASH.fits', memmap=True)
# MS_ind= np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/MS_individual_chemin.txt', usecols=(8,), unpack=True)
# MS_ind_smooth= np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/MS_smoothed_chemin.txt', usecols=(9,), unpack=True)
# AGy_ind= np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/AGy_individual_chemin.txt', usecols=(8,), unpack=True)
# AGo_ind= np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/AGo_individual_chemin.txt', usecols=(8,), unpack=True)
# RG_ind= np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/RG_individual_chemin.txt', usecols=(8,), unpack=True)

# print('read in data')
# #puts data into an accessible table
# data=Table(hdulist[1].data)

# print('data in table')
# #creates arrays of the parameters than I need
# F336W=data['F336W']
# F814W=data['F814W']
# F475W=data['F475W']

# #adds the color data
# MS_475=[]
# MS_814=[]
# for i in range(len(MS_ind)):
# 	N=int(MS_ind[i])
# 	MS_475.append(F475W[N])
# 	MS_814.append(F814W[N])

# MSs_336=[]
# MSs_475=[]
# MSs_814=[]
# for i in range(len(MS_ind_smooth)):
# 	N=int(MS_ind_smooth[i])
# 	MSs_336.append(F336W[N])
# 	MSs_475.append(F475W[N])
# 	MSs_814.append(F814W[N])

# AGy_475=[]
# AGy_814=[]
# for i in range(len(AGy_ind)):
# 	N=int(AGy_ind[i])
# 	AGy_475.append(F475W[N])
# 	AGy_814.append(F814W[N])

# AGo_475=[]
# AGo_814=[]
# for i in range(len(AGo_ind)):
# 	N=int(AGo_ind[i])
# 	AGo_475.append(F475W[N])
# 	AGo_814.append(F814W[N])

# RG_475=[]
# RG_814=[]
# for i in range(len(RG_ind)):
# 	N=int(RG_ind[i])
# 	RG_475.append(F475W[N])
# 	RG_814.append(F814W[N])

MS_475, MS_814=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/MS_individual_color.txt', unpack=True)
AGy_475, AGy_814=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/AGy_individual_color.txt', unpack=True)
AGo_475, AGo_814=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/AGo_individual_color.txt', unpack=True)
RG_475, RG_814=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/RG_individual_color.txt', unpack=True)

def color(mag1, mag2):
	return [a-b for a,b in zip(mag1,mag2)]

MS_color=color(MS_475, MS_814)
#MS_color2=color(MSs_336, MSs_475)
AGy_color=color(AGy_475, AGy_814)
AGo_color=color(AGo_475, AGo_814)
RG_color=color(RG_475, RG_814)

color=MS_color+AGy_color+AGo_color+RG_color
print(len(color))
mag=list(MS_814) + list(AGy_814) + list(AGo_814) + list(RG_814)
print(len(mag))

MS_color=[]
MS_F8W=[]
for i in range(len(color)):
	if color[i]<=1:
                MS_color.append(color[i])
                MS_F8W.append(mag[i])

RG_color=[]
RG_F8W=[]
#polygon creates a shape that corresponds to the area of the CMD that contains RGB stars
polygon= Polygon([(2,23), (2.7, 20.4), (7.5,23)])
for i in range(len(color)):
        point= Point(color[i], mag[i])
        if polygon.contains(point)==True:
                RG_color.append(color[i])
                RG_F8W.append(mag[i])

#below tests to see if the star is in the old AGB
AGo_color=[]
AGo_F8W=[]
#polygon creates a shape that corresponds to the area of the CMD that contains old AGB stars
polygon= Polygon([(2.7,20.4), (7.5, 23), (8, 20.5)])
for i in range(len(color)):
        point= Point(color[i], mag[i])
        if polygon.contains(point)==True:
                AGo_color.append(color[i])
                AGo_F8W.append(mag[i])

#below tests to see if the star is in the young AGB
AGy_color=[]
AGy_F8W=[]
#polygon creates a shape that corresponds to the area of the CMD that contains young AGB stars
polygon= Polygon([(3.5,18), (8,18), (2.7, 20.4), (3.5,20.5), (8,20.5)])
for i in range(len(color)):
	point= Point(color[i], mag[i])
	if polygon.contains(point)==True:
		AGy_color.append(color[i])
		AGy_F8W.append(mag[i])

color1 = np.concatenate((MS_color, AGy_color), axis=None)
color2= np.concatenate((color1, AGo_color), axis=None)
color = np.concatenate((color2, RG_color), axis=None)
mag1 = np.concatenate((MS_F8W, AGy_F8W), axis=None)
mag2= np.concatenate((mag1, AGo_F8W), axis=None)
mag = np.concatenate((mag2, RG_F8W), axis=None)
#plotting CMD
# fig, ax=plt.subplots(1)
# ax.set_rasterization_zorder(1)
# for axis in ['top','bottom','left','right']:
#         ax.spines[axis].set_linewidth(2)
# ax.tick_params(axis='x',which='minor',bottom='off')
# plt.tick_params(which='both', width=2)
# plt.tick_params(which='major', length=7)
# plt.tick_params(which='minor', length=4)
# plt.tick_params(labelsize=10) 
# plt.minorticks_on()
rc('font', family = 'serif', size = 13)
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
plt.scatter(MS_color, MS_F8W, color='b', alpha=0.3, s=14, label='MS', zorder=0)
plt.scatter(AGy_color, AGy_F8W, color='m', alpha=0.3, s=14, label='young AGB',zorder=0)
plt.scatter(AGo_color, AGo_F8W, color='green', alpha=0.3, s=14, label='older AGB',zorder=0, marker='s')
plt.scatter(RG_color, RG_F8W, color='r', alpha=0.2, s=14, label='RGB',zorder=0)
#ax.annotate('MS', xy=(1.7,21), horizontalalignment='right', fontsize=13)
#ax.annotate('young AGB', xy=(5,19.3), horizontalalignment='right', fontsize=13)
#ax.annotate('older AGB', xy=(7.6,21), horizontalalignment='right', fontsize=13)
#ax.annotate('RGB', xy=(4,22.5), horizontalalignment='right', fontsize=13)
#im=ax.hist2d(color, mag, bins=900, cmap=plt.cm.binary,  norm=LogNorm())
#plt.colorbar(im[3], ax=ax)
plt.xlabel(r'$\rm F475W - F814W$', fontsize=14)
plt.xlim(-1,8)
plt.ylim(18,23)
plt.gca().invert_yaxis()
plt.ylabel(r'$\rm F814W$', fontsize=14)
#plt.legend(frameon=False, fontsize=10)
plt.savefig('/Users/amandaquirk/Desktop/CMD.pdf', bbox_inches='tight')
plt.close()

# plt.scatter(MS_color2, MSs_475)
# plt.xlabel('F336W - F475W (mag)', fontsize=16)
# plt.xlim(-2,2)
# plt.ylim(18,23)
# plt.gca().invert_yaxis()
# plt.ylabel('F475W (mag)', fontsize=16)
# plt.savefig('/Users/amandaquirk/Desktop/MS_CMD.png')

# #save files
# file=open('/Users/amandaquirk/Desktop/MS_individual_color.txt', 'w')
# file.write('#F475W, F814W\n')
# for i in range(len(MS_color)):
# 	file.write('{} {}\n'.format(MS_475[i], MS_814[i]))
# file.close()

# file=open('/Users/amandaquirk/Desktop/AGy_individual_color.txt', 'w')
# file.write('#F475W, F814W\n')
# for i in range(len(AGy_color)):
# 	file.write('{} {}\n'.format(AGy_475[i], AGy_814[i]))
# file.close()

# file=open('/Users/amandaquirk/Desktop/AGo_individual_color.txt', 'w')
# file.write('#F475W, F814W\n')
# for i in range(len(AGo_color)):
# 	file.write('{} {}\n'.format(AGo_475[i], AGo_814[i]))
# file.close()

# file=open('/Users/amandaquirk/Desktop/RG_individual_color.txt', 'w')
# file.write('#F475W, F814W\n')
# for i in range(len(RG_color)):
# 	file.write('{} {}\n'.format(RG_475[i], RG_814[i]))
# file.close()

# file=open('/Users/amandaquirk/Desktop/MS_smooth_color.txt', 'w')
# file.write('#F336W, F475W, F814W\n')
# for i in range(len(MS_color2)):
# 	file.write('{} {} {}\n'.format(MSs_336[i], MSs_475[i], MSs_814[i]))
# file.close()







