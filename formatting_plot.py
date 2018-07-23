import matplotlib.pyplot as plt 
import numpy as np

MS_eta, MS_xi= np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/MS_smoothed_goodcenters.txt', usecols=(0,1,), unpack=True)

plt.figure(figsize=(4, 6))
plt.scatter(MS_xi, MS_eta, s=1)#, c=MS_v, cmap='rainbow',vmin=-300,vmax=0)
plt.xlim(-.5,1)
plt.gca().invert_xaxis()
plt.ylim(-0.75,1.2)
plt.xlabel('xi (deg)')
plt.ylabel('eta (deg)')
#clb=plt.colorbar()
#clb.set_label('v (km/s)')
plt.savefig('/Users/amandaquirk/Desktop/MS_sizetesting.png')
plt.close()