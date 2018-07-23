import numpy as np 
import matplotlib.pyplot as plt 
from astropy.io import fits 
from astropy.table import Table 
from astropy import units as u
from astropy.coordinates import SkyCoord

img=plt.imread('/Users/amandaquirk/Documents/AsymmetricDrift/Plots/surveys_map/tilted_map_black1.jpeg')
hdulist= fits.open('/Users/amandaquirk/Documents/AsymmetricDrift/Data/subMasterSPLASH.fits', memmap=True)

#puts data into an accessible table
data=Table(hdulist[1].data)

#creates arrays of the parameters than I need
RA=data['RA']
Dec=data['DEC']

m31 = SkyCoord(10.6847083*u.deg, 41.26875*u.deg, frame='icrs')
c=SkyCoord(ra=RA, dec=Dec, frame='icrs', unit=(u.hourangle,u.deg))
c_inm31=c.transform_to(m31.skyoffset_frame())
xi, eta=c_inm31.lon, c_inm31.lat
xi=xi.degree
eta=eta.degree

#converts xi (kpc) and eta (kpc) to x in the shifted coordinates frame
def x_deg(xi, eta):
	sine=np.sin(float(45*np.pi) / 180)
	cosine=np.cos(float(45*np.pi) / 180)
	x=[(a*cosine)-(b*sine) for a,b in zip(xi, eta)]
	return x 

#converts xi (kpc) and eta (kpc) to y in the shifted coordinates frame
def y_deg(xi, eta):
	sine=np.sin(float(45*np.pi) / 180)
	cosine=np.cos(float(45*np.pi) / 180)
	y=[(a*cosine)+(b*sine) for a,b in zip(eta, xi)]
	return y 

# x=[a*13.67 +1 for a in x_deg(xi,eta)]
# y=[a*13.67 for a in y_deg(xi,eta)]

# x_good=[a for a,b in zip(x,y) if a<=7 and a>=-3 and b<=19 and b>=-19]
# y_good=[b for a,b in zip(x,y) if a<=7 and a>=-3 and b<=19 and b>=-19]

xi_kpc=[a*13.67 for a in xi]
eta_kpc=[a*13.67 for a in eta]

xi_good=[a +1.5 for a,b in zip(xi_kpc, eta_kpc) if b > (16/13) * a -12 and b < (25/15) * a + 15 - (25/3) and b<15]
eta_good=[b +0.7 for a,b in zip(xi_kpc, eta_kpc) if b > (16/13) * a -12 and b < (25/15) * a + 15 - (25/3) and b<15]
#x_HI=[42.377,42.377, -42.377, -42.377]
#y_HI=[38.276, -43.744, -43.744, 38.276]
x1=[.55, 15]
y1=[.2, 15]
x2=[.55, 12]
y2=[.2,17]

fig, ax=plt.subplots(1)
for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(2)
ax.tick_params(axis='x',which='minor',bottom='off')
plt.tick_params(which='both', width=2)
plt.tick_params(which='major', length=7)
plt.tick_params(which='minor', length=4)
plt.tick_params(labelsize=10) 
plt.minorticks_on()
plt.imshow(img, extent=[15, -15, -14.1, 16], zorder=0)#[max(xi_good),-max(xi_good),-max(eta_good), max(eta_good)])
plt.plot(x1,y1, c='aqua', linewidth=1, zorder=1)
plt.plot(x2,y2, c='aqua', linewidth=1, zorder=1)
plt.scatter(xi_good, eta_good, c='crimson', alpha=0.15, s=1, zorder=2)
#plt.plot(x_HI, y_HI, c='g')
plt.xlim(15,-15)
plt.ylim(-14.1,16)
plt.xlabel(r'$\xi\ (kpc)$', fontsize=13)
plt.ylabel(r'$\eta\ (kpc)$', fontsize=13)
plt.savefig('/Users/amandaquirk/Desktop/surveys_no_lines.pdf', bbox_inches='tight')




