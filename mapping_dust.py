import numpy as np 
import matplotlib.pyplot as plt
import matplotlib

def mapping_dust(x,y,dust,age):
	plt.figure(figsize=(4, 6))
	cmap=plt.cm.RdYlBu
	norm = matplotlib.colors.BoundaryNorm(np.linspace(0,3.5,9), cmap.N)
	plt.scatter(x, y, c=dust, s=8,cmap=cmap,norm=norm,edgecolor='none')
	# clb=plt.colorbar(ticks=np.linspace(1,5,5))
	plt.xlim(2.5,12,5)
	plt.gca().invert_xaxis()
	plt.xlabel('xi (kpc)')
	plt.ylabel('eta (kpc)')
	plt.ylim(-2.5,15)
	clb=plt.colorbar()
	clb.set_label('Extinction')
	plt.savefig('/Users/amandaquirk/Desktop/mappying_dust_{}.png'.format(age))

MS_x, MS_y, MS_dust=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/MS_dust.txt', usecols=(0,1,2,), unpack=True)
AGy_x, AGy_y, AGy_dust=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/AGy_dust.txt', usecols=(0,1,2,), unpack=True)
AGo_x, AGo_y, AGo_dust=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/AGo_dust.txt', usecols=(0,1,2,), unpack=True)
RG_x, RG_y, RG_dust=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/RG_dust.txt', usecols=(0,1,2,), unpack=True)

mapping_dust(MS_x, MS_y, MS_dust,'MS')
mapping_dust(AGy_x, AGy_y, AGy_dust,'AGy')
mapping_dust(AGo_x, AGo_y, AGo_dust,'AGo')
mapping_dust(RG_x, RG_y, RG_dust,'RG')