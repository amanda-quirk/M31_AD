import numpy as np 
from astropy import units as u
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt 
from matplotlib import rc

RA1, RA2, RA3, Dec1, Dec2, Dec3, HI_IWM = np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/list-veliwm.ascii', usecols=(0, 1, 2, 3, 4, 5, 8,), unpack=True)

#reformat the coordinates
RAs = []
Decs = []
m31 = SkyCoord(10.6847083*u.deg, 41.26875*u.deg, frame='icrs')
for i in range(len(RA1)):
	RAs.append('{}:{}:{}'.format(int(RA1[i]), int(RA2[i]), RA3[i]))
	Decs.append('{}:{}:{}'.format(int(Dec1[i]), int(Dec2[i]), Dec3[i]))
	
#convert to xi and eta
c = SkyCoord(ra = RAs, dec = Decs, frame='icrs', unit=(u.hourangle,u.deg))
c_inm31 = c.transform_to(m31.skyoffset_frame())
xi, eta = c_inm31.lon, c_inm31.lat


def calculate_AD(agebin):

	HImain, HIclose, age = np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/velocityhi-coords-amanda.dat', usecols=(7, 8, 10,), unpack=True)

	#sort into age bins
	if agebin == 'MS':
		age_tag = 1
	elif agebin == 'AGy':
		age_tag = 2
	elif agebin == 'AGo':
		age_tag = 3
	elif agebin == 'RG':
		age_tag = 4

	agegroup = age == age_tag

	xis = xi[agegroup]
	etas = eta[agegroup]
	HI_avg = HI_IWM[agegroup]
	HImains = HImain[agegroup]
	HIcloses = HIclose[agegroup]

	xi0, eta0, HImain0, HIclose0 = np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/{}_smoothed_chemin.txt'.format(agebin), usecols=(0, 1, 6, 7), unpack=True) #kpc
	vrot_star = np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/{}_master_vrot.txt'.format(agebin), usecols=(2,), unpack=True) 

	#match stars I actually used in my analysis
	xi_keep = [] #deg
	eta_keep = []
	HI_keep = []
	vrot_star_keep = []
	for j in range(len(xi0)): 
		number = []
		for i in range(len(xis)):
			if (HImain0[j] - HImains[i] == 0) and (HIclose0[j] - HIcloses[i] ==0) and abs(xi0[j] - xis[i].value * 13.67) == 0 and abs(eta0[j] - etas[i].value * 13.67) == 0:
				number.append(j)
		if len(number) < 2: #what to do with the repeated values?
			xi_keep.append(xis.value[i])
			eta_keep.append(etas.value[i])
			HI_keep.append(HI_avg[i])
			vrot_star_keep.append(vrot_star[j])

	xi_keep = np.array((xi_keep))
	eta_keep = np.array((eta_keep))

	def shifted_coords(xi, eta): #deg
		sine = np.sin((37 * np.pi) / 180)
		cosine=np.cos((37* np.pi) / 180)
		x = (xi * cosine) - (eta * sine)
		y = (eta * cosine) + (xi * sine)
		return x, y

	shifted_x, shifted_y = shifted_coords(xi_keep, eta_keep)

	def distance(x, y):
		inclination_factor = np.cos((77 * np.pi) / 180)**2
		ang_distance_sq=[(a**2) + ((b**2) / inclination_factor) for a,b in zip(y,x)]
		ang_dist=[np.sqrt(a) for a in ang_distance_sq]
		dist=[a * 13.67 for a in ang_dist]
		return dist

	r = distance(shifted_x, shifted_y)

	def PA(x,y): #defined as west of north as positioned in maps (inverted x axis)
		deg = np.zeros_like(x)
		for i in range(len(x)):
			if x[i] > 0:
				rad=np.arctan((y[i]) / x[i])
				deg[i]=90 - ((rad * 180) / np.pi)
			else:
				rad=np.arctan((y[i]) / x[i])
				deg[i]=270 - ((rad*180) / np.pi)
		return deg + 37

	PAs = PA(shifted_x, shifted_y)

	HI_r, HI_PA, HI_i, HI_v=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/HI_PA_i_vrot.txt', unpack=True)
	
	#below defines a function to determine which HI ring a star is cloest to
	def find_nearest_ring(radius):
		idx=[]
		for i in range(len(HI_r)):
			idx.append(np.abs(HI_r[i] - radius))
			n=np.argmin(idx)
		return n #index of the radius closest to the star's radius
	
	#assigning a PA  and i to each star and HI component
	PA_ring=[]
	i_ring=[]
	for j in range(len(r)):
		N = find_nearest_ring(r[j])
		PA_ring.append(HI_PA[N])
		i_ring.append(HI_i[N])

	def Vrot_tilted_ring(v,PA_ring,PA_star, i_ring): 
		vsys = -300 #km/s, as defined in Claire's thesis
		A = [float(a-vsys)/np.sin(float(b*np.pi) / 180) for a, b in zip(v, i_ring)]
		B = [(np.tan(float((a-b)*np.pi) / 180))**2 for a,b in zip(PA_ring, PA_star)]
		C = [(np.cos(float(a*np.pi) / 180))**2 for a in i_ring]
		rotation_velocity= [a*np.sqrt(1+(float(b)/c)) for a,b,c in zip(A,B,C)]
		positive=[np.absolute(a) for a in rotation_velocity]
		return positive

	HI_vrot = Vrot_tilted_ring(HI_keep, PA_ring, PAs, i_ring)

	return np.array((HI_vrot)) - np.array((vrot_star_keep))


MS_AD = calculate_AD('MS')
AGy_AD = calculate_AD('AGy')
AGo_AD = calculate_AD('AGo')
RG_AD = calculate_AD('RG')

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
# plt.hist(MS_AD, bins=range(-200, 300, 20), label='MS={} stars'.format(len(MS_AD)),normed=1, histtype='step', linewidth=1.6,linestyle='--',stacked=True,fill=False, color='b')
# plt.hist(AGy_AD,  bins=range(-200, 300, 20),label='young AGB={} stars'.format(len(AGy_AD)),normed=1, histtype='step', linewidth=1.6,stacked=True,fill=False, color='m', hatch='//', alpha=0.4)
# plt.hist(AGo_AD,  bins=range(-200, 300, 20),label='older AGB={} stars'.format(len(AGo_AD)),normed=1,histtype='step', linewidth=1.6,stacked=True,fill=True, color='green', alpha=0.4)
# plt.hist(RG_AD, bins=range(-200, 300, 20), label='RGB={} stars'.format(len(RG_AD)),normed=1,histtype='step', linewidth=1.6,linestyle='-',stacked=True,fill=False, color='r')
# plt.legend(loc=1, frameon=False, fontsize=10)
# plt.xlim(-200,300)
# plt.xlabel(r'$ \rm Asymmetric\ Drift:\ \itv_{a}\ \rm(km\ s^{-1})$', fontsize=13)
# plt.savefig('/Users/amandaquirk/Desktop/IWM_asymmetric_drift_zoom.pdf')

def errors(AD):
	result = np.percentile(AD, [16, 50, 84])
	median = result[1]
	lower_error = result[1] - result[0]
	upper_error = result[2] - result[1]
	return median, lower_error, upper_error
	
print(errors(MS_AD))
print(errors(AGy_AD))
print(errors(AGo_AD))
print(errors(RG_AD))


