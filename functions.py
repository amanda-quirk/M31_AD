import matplotlib.pyplot as plt 
import numpy as np 

'''Geometry Stuff
==============================================================================================================================='''
#computes radial distance from shifted coordinates
def distance(x, y):
	inclination_factor=np.cos(float(77*np.pi) / 180)**2
	ang_distance_sq=[(a**2)+(float(b**2)/inclination_factor) for a,b in zip(y,x)]
	ang_dist=[np.sqrt(a) for a in ang_distance_sq]
	dist=[a*13.67 for a in ang_dist]
	return dist

#the function below will assign each star and gas cloud a position angle using their xi and eta coordinates
def PA(x,y): #defined as west of north as positioned in maps (inverted x axis)
	if x>0:
		rad=np.arctan(float(y)/x)
		deg=90-(float(rad*180)/np.pi)
	else:
		rad=np.arctan(float(y)/x)
		deg=270-(float(rad*180)/np.pi)
	return deg + 37

#converts xi (kpc) and eta (kpc) to x in the shifted coordinates frame
def x_kpc(xi, eta):
	xi_deg=[a/13.67 for a in xi]
	eta_deg=[a/13.67 for a in eta]
	sine=np.sin(float(37*np.pi) / 180)
	cosine=np.cos(float(37*np.pi) / 180)
	x=[(a*cosine)-(b*sine) for a,b in zip(xi_deg, eta_deg)]
	return x 

#converts xi (kpc) and eta (kpc) to y in the shifted coordinates frame
def y_kpc(xi, eta):
	xi_deg=[a/13.67 for a in xi]
	eta_deg=[a/13.67 for a in eta]
	sine=np.sin(float(37*np.pi) / 180)
	cosine=np.cos(float(37*np.pi) / 180)
	y=[(a*cosine)+(b*sine) for a,b in zip(eta_deg, xi_deg)]
	return y

'''Major Plotting
======================================================================================================================='''
#pretty single plot formatting
def single_plot():
	rc('font', family = 'serif')
	fig, ax=plt.subplots(1)
	for axis in ['top','bottom','left','right']:
	        ax.spines[axis].set_linewidth(2)
	ax.tick_params(axis='x',which='both',bottom='on',top='on', direction='out')
	ax.tick_params(axis='x',which='both',top='on', direction='in')
	ax.tick_params(axis='y',which='both',left='on',top='off', direction='out')
	ax.tick_params(axis='y',which='both',right='on', direction='in')
	plt.tick_params(which='both', width=2)
	plt.tick_params(which='major', length=7)
	plt.tick_params(which='minor', length=4)
	plt.tick_params(labelsize=12) 
	plt.minorticks_on()

#velocity dispersion maps
def dispersion_map(xi, eta, dispersion):
	plt.figure(figsize=(4, 6))
	plt.scatter(xi, eta, s=1, c=dispersion, cmap='copper',vmin=10,vmax=150)
	plt.xlim(2.5,12,5)
	plt.gca().invert_xaxis()
	plt.ylim(-2.5,15)
	plt.xlabel('xi (kpc)')
	plt.ylabel('eta (kpc)')
	clb=plt.colorbar()
	clb.set_label('v (km/s)')
	plt.close()

#RCs for dusty and/or clear LOSs
def dust_rotation_curve(Av, r, vrot, HI_vrot, age, dusty=True):
	if age=='MS':
		color='b'
	if age=='RG':
		color='r'
	if age!='MS' and age !='RG':
		color='m'
	r_dust=[]
	r_clear=[]
	vrot_dust=[]
	vrot_clear=[]
	HI_vrot_clear=[]
	HI_vrot_dust=[]
	for i in range(len(r)):
		if Av[i]<1:
			r_clear.append(r[i])
			vrot_clear.append(vrot[i])
			HI_vrot_clear.append(HI_vrot[i])
		else:
			r_dust.append(r[i])
			vrot_dust.append(vrot[i])
			HI_vrot_dust.append(HI_vrot[i])
	if dusty==True:
		description='dusty'
		x=r_dust
		star_y=vrot_dust
		HI_y=HI_vrot_dust
	if dusty==False:
		description='translucent'
		x=r_clear
		star_y=vrot_clear
		HI_y=HI_vrot_clear 

	plt.scatter(x, star_y, s=1, c='{}'.format(color))
	plt.scatter(x, HI_y, s=1, c='k')
	plt.xlim(4,20)
	plt.ylim(100,300)
	plt.xlabel('r [kpc]')
	plt.ylabel('rotational v [km/s]')
	plt.savefig('/Users/amandaquirk/Desktop/{}_rotationcurve_{}.png'.format(description,age))
	plt.close()

#general rotation curve
def rotation_curve(r, star_vrot, HI_vrot, age, description):
	if age=='MS':
		color='b'
	if age=='RG':
		color='r'
	if age!='MS' and age !='RG':
		color='m'
	plt.scatter(r, star_vrot, s=1, c='{}'.format(color)) #should put quotes around the bracket?
	plt.scatter(r, HI_vrot, s=1, c='k')
	plt.xlim(4,20)
	plt.ylim(100,300)
	plt.xlabel('r [kpc]')
	plt.ylabel('rotational v [km/s]')
	plt.savefig('/Users/amandaquirk/Desktop/{}_vrot_{}.png'.format(age, description))
	plt.close()

#velocity maps
def velocity_maps(xi, eta, velocity):
	plt.figure(figsize=(4, 6))
	plt.scatter(xi, eta, s=1, c=velocity, cmap='rainbow',vmin=-300,vmax=0)
	plt.xlim(2.5,12,5)
	plt.gca().invert_xaxis()
	plt.ylim(-2.5,15)
	plt.xlabel('xi (kpc)')
	plt.ylabel('eta (kpc)')
	clb=plt.colorbar()
	clb.set_label('v (km/s)')
	plt.close()

'''Misc.
========================================================================================================================'''
#calculates the asymmetric drift
def asymmetric_drift(v_star, v_gas):
	v_gas = np.array((v_gas))
	v_star = np.array((v_star))
	return v_gas-v_star


def color(mag1, mag2):
	return [a-b for a,b in zip(mag1,mag2)]

#finds median lines for rotation curves
def median_line(r, v):

	delta_r = 0.5 #kpc
	bins = np.arange(5, 20, delta_r)
	median_r = bins - delta_r / 2
	
	medians = np.zeros_like(bins)

	for i in range(len(bins)):
		data = [b for a,b in zip(r, v) if a < bins[i] and a > bins[i] - delta_r]
		if len(data) > 5:
			medians[i] = np.median(data)
		else:
			medians[i] = np.nan

	data = (np.isnan(medians) == False)

	good_r = median_r[data]
	good_v = medians[data]

	return good_r, good_v

'''Rotation Velocities
========================================================================================================================'''
#I define a functionto calculate the rotation speed without a tiled ring model
def Vrot(v, PA_star):
	sine=np.sin(float(77*np.pi) / 180)
	inclination_factor=(np.cos(float(77*np.pi) / 180))**2
	vsys= -300 #km/s, as defined in Claire's thesis
	A=[float(a-vsys)/sine for a in v]
	B= [(np.tan(float((37-b)*np.pi) / 180))**2 for b in PA_star]
	rotation_velocity= [a*np.sqrt(1+(float(b)/inclination_factor))for a,b in zip(A,B)]
	positive=[np.absolute(a) for a in rotation_velocity]
	return positive

#I define a function to caluclate the rotation speed using the tiled ring model
def Vrot_tilted_ring(v,PA_ring,PA_star, i_ring): 
	vsys= -300 #km/s, as defined in Claire's thesis
	A=[float(a-vsys)/np.sin(float(b*np.pi) / 180) for a, b in zip(v, i_ring)]
	B= [(np.tan(float((a-b)*np.pi) / 180))**2 for a,b in zip(PA_ring, PA_star)]
	C= [(np.cos(float(a*np.pi) / 180))**2 for a in i_ring]
	rotation_velocity= [a*np.sqrt(1+(float(b)/c)) for a,b,c in zip(A,B,C)]
	positive=[np.absolute(a) for a in rotation_velocity]
	return positive

'''Stats Stuff
==================================================================================================================='''
#function to normalize the weights
def normed_weight(w):
	sum_weights=sum(w)
	return w / sum_weights 

#function to calculate weights from errors
def weights(err):
	return 1 / (err**2)

#function does the weighted mean
def weighted_mean(data,norm_w):
	return sum([a*b for a,b in zip(data, norm_w)])

#function does the weighted RMSE
def weighted_rmse(norm_w, data, mean):
	differences=[x - mean for x in data]
	diff_sq=[d**2 for d in differences]
	products=[a*b for a,b in zip(diff_sq, norm_w)]
	return np.sqrt(sum(products))







