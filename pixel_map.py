import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

#read in data-- eta and xi are in degrees
MS_xi, MS_eta, MS_v, MS_err, MS_n, MS_HImain, MS_HIclose=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/MS_individual_chemin.txt', usecols=(0,1,2,3,4,5,6,), unpack=True)

#convert ADJUSTED errors to weights
#function to calculate the weights
def weights(err):
	return 1 / (err**2)

#function to normalize the weights
def normed_weight(w):
	sum_weights=sum(w)
	return w / sum_weights  

MS_weights=weights(MS_err)

#function does the weighted meean
def weighted_mean(data,norm_w):
	return sum([a*b for a,b in zip(data, norm_w)])

#function does the weighted RMSE
def weighted_rMSe(norm_w, data, mean):
	differences=[x - mean for x in data]
	diff_sq=[d**2 for d in differences]
	products=[a*b for a,b in zip(diff_sq, norm_w)]
	return np.sqrt(sum(products))

#set up the grid
n=150/3600 #pixel size in deg
rows=np.linspace(min(MS_xi), max(MS_xi)+n, int((max(MS_xi)-min(MS_xi))/ (n) +1))
columns=np.linspace(min(MS_eta), max(MS_eta)+n, int((max(MS_eta)-min(MS_eta))/ (n) +1))  

#adding grid on top of original scatter plot for visual check
fig=plt.figure(figsize=(4, 6))
ax = fig.gca()
ax.set_xticks(rows)
ax.set_yticks(columns)
plt.scatter(MS_xi, MS_eta, s=3, c=MS_v, cmap='rainbow',vmin=-300,vmax=0)
plt.gca().invert_xaxis()
plt.xlabel('xi (deg)')
plt.ylabel('eta (deg)')
clb=plt.colorbar()
clb.set_label('v (km/s)')
plt.grid()
plt.savefig('/Users/amandaquirk/Desktop/MS_grid.png')
plt.close()

#makes 2D histogram
ret=stats.binned_statistic_2d(MS_xi, MS_eta, None, 'count', bins=[rows, columns])

#this gives the bin number that each data point is put into
N=ret.binnumber

#will have averaged data values for each grid point
#NOTE! right now I am averaging the positions and HI data-- I dont think this is ideal but is fine
MS_v_pix=[]
MS_dispersion=[]
MS_xi_pix=[] #kpc
MS_eta_pix=[] #kpc
MS_n_pix=[]
MS_HImain_pix=[]
MS_HIclose_pix=[]
MS_err_pix=[]

#goes through each bin number (numbering scheme is not apparent to me-- doesn't start at 0 or go sequentially)
i_start=min(N)
i=i_start
while i < max(N)+1:
	weights=[]
	cell_v=[]
	cell_xi=[]
	cell_eta=[]
	cell_n=[]
	cell_HImain=[]
	cell_HIclose=[]
	cell_err=[]
	for j in range(len(N)):
		if N[j]==i: #groups all of the data points with the same bin number together
			weights.append(MS_weights[j]) 
			cell_v.append(MS_v[j])
			cell_xi.append(MS_xi[j])
			cell_eta.append(MS_eta[j])
			cell_n.append(MS_n[j])
			cell_HImain.append(MS_HImain[j])
			cell_HIclose.append(MS_HIclose[j])
			cell_err.append(MS_err[j])
	if len(cell_v)>4: #eliminates grid points that don't have enough neighbors --should be 15 but that's too big
		weights_norm=normed_weight(weights)
		mean_v=weighted_mean(cell_v, weights_norm)
		MS_v_pix.append(mean_v)
		MS_xi_pix.append(weighted_mean(cell_xi, weights_norm)*13.67) #converts to kpc 
		MS_eta_pix.append(weighted_mean(cell_eta, weights_norm)*13.67) #converts to kpc
		MS_n_pix.append(weighted_mean(cell_n, weights_norm))
		MS_HImain_pix.append(weighted_mean(cell_HImain, weights_norm))
		MS_HIclose_pix.append(weighted_mean(cell_HIclose, weights_norm))
		MS_dispersion.append(weighted_rMSe(weights_norm, cell_v, mean_v))
		MS_err_pix.append(np.mean(cell_err))
	i=i+1

#examining how things are binned
#print(ret.statistic)
print(len(MS_v_pix))

#plotting final map
plt.figure(figsize=(4, 6))
plt.scatter(MS_xi_pix, MS_eta_pix, s=3, c=MS_v_pix, cmap='rainbow',vmin=-300,vmax=0)
plt.gca().invert_xaxis()
plt.xlabel('xi (kpc)')
plt.ylabel('eta (kpc)')
clb=plt.colorbar()
clb.set_label('v (km/s)')
plt.savefig('/Users/amandaquirk/Desktop/MS_pixelated_map.png')
plt.close()

#saving to file
file=open('/Users/amandaquirk/Documents/AsymmetricDrift/Data/MS_pixelated_map.txt', 'w')
file.write('#xi (kpc), eta (kpc), average v(km/s), verr, dispersion, n, HI main, HI close\n')
for i in range(len(MS_xi_pix)):
	file.write('{} {} {} {} {} {} {} {}\n'.format(MS_xi_pix[i],MS_eta_pix[i], MS_v_pix[i], MS_err_pix[i], MS_dispersion[i], MS_n_pix[i], MS_HImain_pix[i], MS_HIclose_pix[i]))
file.close()












