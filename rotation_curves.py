import numpy as np 
import matplotlib.pyplot as plt 
from astropy import units as u
from astropy.coordinates import SkyCoord

#first we will import the data of the star info for each age bin
MSs_eta, MSs_xi, MSs_v= np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/MS_smoothed_goodcenters.txt', unpack=True)
RGs_eta, RGs_xi, RGs_v= np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/RG_smoothed_goodcenters.txt', unpack=True)
AGos_eta, AGos_xi, AGos_v= np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/AGo_smoothed_goodcenters.txt', unpack=True)
AGys_eta, AGys_xi, AGys_v= np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/AGy_smoothed_goodcenters.txt', unpack=True)

#individual velocities
MS_eta, MS_xi, MS_v= np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/MS_individualv.txt', unpack=True) #eliminate the s here when done with testing
RG_eta, RG_xi, RG_v= np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/RG_individualv.txt', unpack=True)
AGo_eta, AGo_xi, AGo_v= np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/AGo_individualv.txt', unpack=True)
AGy_eta, AGy_xi, AGy_v= np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/AGy_individualv.txt', unpack=True)

#first we will calculate the distance of each star from the center

#below defines a function that converts eta and xi and then uses the shifted coordinates to find the deprojected radius. to calculate the deprojected radius from the center, we must shift the xi and eta according to the PA of the disk. (This is a simple shift of coordinate axes.) We will assume the disc has a PA of 37 deg and an inclination of 77 deg
def x(xi, eta):
	sine=np.sin(float(37*np.pi) / 180)
	cosine=np.cos(float(37*np.pi) / 180)
	inclination_factor=(np.cos(float(77*np.pi) / 180))**2
	x=(xi*cosine)-(eta*sine)
	return x 

def y(xi, eta):
	sine=np.sin(float(37*np.pi) / 180)
	cosine=np.cos(float(37*np.pi) / 180)
	inclination_factor=(np.cos(float(77*np.pi) / 180))**2
	y=(eta*cosine)+(xi*sine)
	return y

#below will contain x and y coordinates for all of the star data
#individual
MS_x=[]
MS_y=[]
for i in range(len(MS_xi)):
	x_coord=x(MS_xi[i], MS_eta[i])
	MS_x.append(x_coord)
	y_coord=y(MS_xi[i], MS_eta[i])
	MS_y.append(y_coord)
AGy_x=[]
AGy_y=[]
for i in range(len(AGy_xi)):
	x_coord=x(AGy_xi[i], AGy_eta[i])
	AGy_x.append(x_coord)
	y_coord=y(AGy_xi[i], AGy_eta[i])
	AGy_y.append(y_coord)
AGo_x=[]
AGo_y=[]
for i in range(len(AGo_xi)):
	x_coord=x(AGo_xi[i], AGo_eta[i])
	AGo_x.append(x_coord)
	y_coord=y(AGo_xi[i], AGo_eta[i])
	AGo_y.append(y_coord)
RG_x=[]
RG_y=[]
for i in range(len(RG_xi)):
	x_coord=x(RG_xi[i], RG_eta[i])
	RG_x.append(x_coord)
	y_coord=y(RG_xi[i], RG_eta[i])
	RG_y.append(y_coord)

#smoothed
MSs_x=[]
MSs_y=[]
for i in range(len(MSs_xi)):
	x_coord=x(MSs_xi[i], MSs_eta[i])
	MSs_x.append(x_coord)
	y_coord=y(MSs_xi[i], MSs_eta[i])
	MSs_y.append(y_coord)
AGys_x=[]
AGys_y=[]
for i in range(len(AGys_xi)):
	x_coord=x(AGys_xi[i], AGys_eta[i])
	AGys_x.append(x_coord)
	y_coord=y(AGys_xi[i], AGys_eta[i])
	AGys_y.append(y_coord)
AGos_x=[]
AGos_y=[]
for i in range(len(AGos_xi)):
	x_coord=x(AGos_xi[i], AGos_eta[i])
	AGos_x.append(x_coord)
	y_coord=y(AGos_xi[i], AGos_eta[i])
	AGos_y.append(y_coord)
RGs_x=[]
RGs_y=[]
for i in range(len(RGs_xi)):
	x_coord=x(RGs_xi[i], RGs_eta[i])
	RGs_x.append(x_coord)
	y_coord=y(RGs_xi[i], RGs_eta[i])
	RGs_y.append(y_coord)

def distance(x, y):
	inclination_factor=np.cos(float(77*np.pi) / 180)**2
	ang_distance_sq= (y**2)+(float(x**2)/inclination_factor)
	ang_dist=np.sqrt(ang_distance_sq) #this is a distance in degrees. we must now convert to kpc
	dist=ang_dist*13.67
	return dist

#below will contain the deprojected distances of each star
#individual
MS_r=[]
for i in range(len(MS_x)):
	true_dist=distance(MS_x[i], MS_y[i])
	MS_r.append(true_dist)

RG_r=[]
for i in range(len(RG_x)):
	true_dist=distance(RG_x[i], RG_y[i])
	RG_r.append(true_dist)

AGy_r=[]
for i in range(len(AGy_x)):
	true_dist=distance(AGy_x[i], AGy_y[i])
	AGy_r.append(true_dist)

AGo_r=[]
for i in range(len(AGo_x)):
	true_dist=distance(AGo_x[i], AGo_y[i])
	AGo_r.append(true_dist)

#smoothed
MSs_r=[]
for i in range(len(MSs_x)):
	true_dist=distance(MSs_x[i], MSs_y[i])
	MSs_r.append(true_dist)

RGs_r=[]
for i in range(len(RGs_x)):
	true_dist=distance(RGs_x[i], RGs_y[i])
	RGs_r.append(true_dist)

AGys_r=[]
for i in range(len(AGys_x)):
	true_dist=distance(AGys_x[i], AGys_y[i])
	AGys_r.append(true_dist)

AGos_r=[]
for i in range(len(AGos_x)):
	true_dist=distance(AGos_x[i], AGos_y[i])
	AGos_r.append(true_dist)

#the function below will assign each star a position angle using their xi and eta coordinates
def PA(x,y):
	if x>0:
		rad=np.arctan(float(y)/x)
		deg=90-float(rad*180)/np.pi
	else:
		rad=np.arctan(float(y)/x)
		deg=270-(float(rad*180)/np.pi)
	return deg + 37

#individual
MS_PA=[]
for i in range(len(MS_x)):
	pa=PA(MS_x[i], MS_y[i])
	MS_PA.append(pa)

RG_PA=[]
for i in range(len(RG_x)):
	pa=PA(RG_x[i], RG_y[i])
	RG_PA.append(pa)

AGy_PA=[]
for i in range(len(AGy_x)):
	pa=PA(AGy_x[i], AGy_y[i])
	AGy_PA.append(pa)

AGo_PA=[]
for i in range(len(AGo_x)):
	pa=PA(AGo_x[i], AGo_y[i])
	AGo_PA.append(pa)

#smoothed
MSs_PA=[]
for i in range(len(MSs_x)):
	pa=PA(MSs_x[i], MSs_y[i])
	MSs_PA.append(pa)

RGs_PA=[]
for i in range(len(RGs_x)):
	pa=PA(RGs_xi[i], RGs_eta[i])
	RGs_PA.append(pa)

AGys_PA=[]
for i in range(len(AGys_xi)):
	pa=PA(AGys_x[i], AGys_y[i])
	AGys_PA.append(pa)

AGos_PA=[]
for i in range(len(AGos_x)):
	pa=PA(AGos_x[i], AGos_y[i])
	AGos_PA.append(pa)

#now we need to assign each star a PA_ring and i_ring using HI data. we define these values by seeing which HI ring a star is cloest to and assigning it the corresponding PA_ring and i_ring value (HI_PA and HI_i respectively). we also will pair a HI velocity to each star to later calculate the asymmetric drift
#reads in the HI data
HI_r, HI_PA, HI_i, HI_v=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/HI_PA_i_vrot.txt', unpack=True)

#below defines a function to determine which HI ring a star is cloest to
def find_nearest_ring(radius):
	idx=[]
	for i in range(len(HI_r)):
		idx.append(np.abs(HI_r[i]-radius))
		n=np.argmin(idx)
	return n #index of the radius closest to the star's radius

#individual
MS_PA_ring=[]
MS_i=[]
for j in range(len(MS_r)):
	N=find_nearest_ring(MS_r[j])
	MS_PA_ring.append(HI_PA[N])
	MS_i.append(HI_i[N])

RG_PA_ring=[]
RG_i=[]
for j in range(len(RG_r)):
	N=find_nearest_ring(RG_r[j])
	RG_PA_ring.append(HI_PA[N])
	RG_i.append(HI_i[N])

AGy_PA_ring=[]
AGy_i=[]
for j in range(len(AGy_r)):
	N=find_nearest_ring(AGy_r[j])
	AGy_PA_ring.append(HI_PA[N])
	AGy_i.append(HI_i[N])

AGo_PA_ring=[]
AGo_i=[]
for j in range(len(AGo_r)):
	N=find_nearest_ring(AGo_r[j])
	AGo_PA_ring.append(HI_PA[N])
	AGo_i.append(HI_i[N])

#smoothed 
MSs_PA_ring=[]
MSs_i=[]
# MSs_HIv=[] #used later to calculate the asyemmtric drift
for j in range(len(MSs_r)):
	N=find_nearest_ring(MSs_r[j])
	MSs_PA_ring.append(HI_PA[N])
	MSs_i.append(HI_i[N])
	# MSs_HIv.append(HI_v[N])

RGs_PA_ring=[]
RGs_i=[]
# RGs_HIv=[]
for j in range(len(RGs_r)):
	N=find_nearest_ring(RGs_r[j])
	RGs_PA_ring.append(HI_PA[N])
	RGs_i.append(HI_i[N])
	# RGs_HIv.append(HI_v[N])

AGys_PA_ring=[]
AGys_i=[]
# AGys_HIv=[]
for j in range(len(AGys_r)):
	N=find_nearest_ring(AGys_r[j])
	AGys_PA_ring.append(HI_PA[N])
	AGys_i.append(HI_i[N])
	# AGys_HIv.append(HI_v[N])

AGos_PA_ring=[]
AGos_i=[]
# AGos_HIv=[]
for j in range(len(AGos_r)):
	N=find_nearest_ring(AGos_r[j])
	AGos_PA_ring.append(HI_PA[N])
	AGos_i.append(HI_i[N])
	# AGos_HIv.append(HI_v[N])

#i want to also smooth the PA in the same way velocity has been smoothed
MSs_sPA=[] #will contain the smoothed PA
for i in range(len(MSs_PA)):
	c1=SkyCoord(MSs_xi[i], MSs_eta[i], unit=(u.deg, u.deg)) #our first point
	positions=[]
	for j in range(len(MSs_PA)):
		c2=SkyCoord(MSs_xi[j], MSs_eta[j], unit=(u.deg, u.deg)) #second point
		sep=c1.separation(c2) #gets the separation of the two points
		if sep.arcsecond<200: #if the two points are close enough, the position is added to an array to later be averaged
			positions.append(MSs_PA[j])
	MSs_sPA.append(np.mean(positions))

RGs_sPA=[] #will contain the smoothed PA
for i in range(len(RGs_PA)):
	c1=SkyCoord(RGs_xi[i], RGs_eta[i], unit=(u.deg, u.deg)) #our first point
	positions=[]
	for j in range(len(RGs_PA)):
		c2=SkyCoord(RGs_xi[j], RGs_eta[j], unit=(u.deg, u.deg)) #second point
		sep=c1.separation(c2) #gets the separation of the two points
		if sep.arcsecond<200: #if the two points are close enough, the position is added to an array to later be averaged
			positions.append(RGs_PA[j])
	RGs_sPA.append(np.mean(positions))

AGos_sPA=[] #will contain the smoothed PA
for i in range(len(AGos_PA)):
	c1=SkyCoord(AGos_xi[i], AGos_eta[i], unit=(u.deg, u.deg)) #our first point
	positions=[]
	for j in range(len(AGos_PA)):
		c2=SkyCoord(AGos_xi[j], AGos_eta[j], unit=(u.deg, u.deg)) #second point
		sep=c1.separation(c2) #gets the separation of the two points
		if sep.arcsecond<275: #if the two points are close enough, the position is added to an array to later be averaged
			positions.append(AGos_PA[j])
	AGos_sPA.append(np.mean(positions))

AGys_sPA=[] #will contain the smoothed PA
for i in range(len(AGys_PA)):
	c1=SkyCoord(AGys_xi[i], AGys_eta[i], unit=(u.deg, u.deg)) #our first point
	positions=[]
	for j in range(len(AGys_PA)):
		c2=SkyCoord(AGys_xi[j], AGys_eta[j], unit=(u.deg, u.deg)) #second point
		sep=c1.separation(c2) #gets the separation of the two points
		if sep.arcsecond<275: #if the two points are close enough, the position is added to an array to later be averaged
			positions.append(AGys_PA[j])
	AGys_sPA.append(np.mean(positions))

#now we need to calculate the rotation speed of the stars using the PA and LOS velocity. I will use the equation in Claire's thesis on page 203 (but the corrected version I have in my evernote) using the PA and i value given by Table 5 in Chemin et al. 2009
#I define a functionto calculate the rotation speed without a tiled ring model
def Vrot(v, PA_star):
	sine=np.sin(float(77*np.pi) / 180)
	inclination_factor=(np.cos(float(77*np.pi) / 180))**2
	vsys= -300 #km/s, as defined in Claire's thesis
	A=float(v-vsys)/sine
	B= (np.tan(float((37-PA_star)*np.pi) / 180))**2 
	rotation_velocity= A*np.sqrt(1+(float(B)/inclination_factor))
	return np.absolute(rotation_velocity)

#I define a function to caluclate the rotation speed using the tiled ring model
def Vrot_tilted_ring(v,PA_ring,PA_star, i_ring): 
	vsys= -300 #km/s, as defined in Claire's thesis
	A=float(v-vsys)/np.sin(float(i_ring*np.pi) / 180)
	B= (np.tan(float((PA_ring-PA_star)*np.pi) / 180))**2
	C= (np.cos(float(i_ring*np.pi) / 180))**2
	rotation_velocity= A*np.sqrt(1+(float(B)/C))
	return np.absolute(rotation_velocity)

#caluclating the rotation velocity using the above defined function
#individual v and PA without tilted ring model
MS_ii_vrot=[]
for i in range(len(MS_r)):
	vrot=Vrot(MS_v[i],MS_PA[i])
	MS_ii_vrot.append(vrot)

RG_ii_vrot=[]
for i in range(len(RG_r)):
	vrot=Vrot(RG_v[i],RG_PA[i])
	RG_ii_vrot.append(vrot)

AGo_ii_vrot=[]
for i in range(len(AGo_r)):
	vrot=Vrot(AGo_v[i],AGo_PA[i])
	AGo_ii_vrot.append(vrot)

AGy_ii_vrot=[]
for i in range(len(AGy_r)):
	vrot=Vrot(AGy_v[i],AGy_PA[i])
	AGy_ii_vrot.append(vrot)

#individual v and PA with tilted ring model
MS_ii_vrot_tilt=[]
for i in range(len(MS_r)):
	vrot=Vrot_tilted_ring(MS_v[i],MS_PA_ring[i],MS_PA[i], MS_i[i])
	MS_ii_vrot_tilt.append(vrot)

RG_ii_vrot_tilt=[]
for i in range(len(RG_r)):
	vrot=Vrot_tilted_ring(RG_v[i],RG_PA_ring[i],RG_PA[i], RG_i[i])
	RG_ii_vrot_tilt.append(vrot)

AGo_ii_vrot_tilt=[]
for i in range(len(AGo_r)):
	vrot=Vrot_tilted_ring(AGo_v[i],AGo_PA_ring[i],AGo_PA[i], AGo_i[i])
	AGo_ii_vrot_tilt.append(vrot)

AGy_ii_vrot_tilt=[]
for i in range(len(AGy_r)):
	vrot=Vrot_tilted_ring(AGy_v[i],AGy_PA_ring[i],AGy_PA[i], AGy_i[i])
	AGy_ii_vrot_tilt.append(vrot)

#smoothed velocity and individual PA without tilted ring model
MS_si_vrot=[]
for i in range(len(MSs_r)):
	vrot=Vrot(MSs_v[i],MSs_PA[i])
	MS_si_vrot.append(vrot)

RG_si_vrot=[]
for i in range(len(RGs_r)):
	vrot=Vrot(RGs_v[i],RGs_PA[i])
	RG_si_vrot.append(vrot)

AGo_si_vrot=[]
for i in range(len(AGos_r)):
	vrot=Vrot(AGos_v[i],AGos_PA[i])
	AGo_si_vrot.append(vrot)

AGy_si_vrot=[]
for i in range(len(AGys_r)):
	vrot=Vrot(AGys_v[i],AGys_PA[i])
	AGy_si_vrot.append(vrot)

#smoothed v and individual PA with tilted ring model
MS_si_vrot_tilt=[]
for i in range(len(MSs_r)):
	vrot=Vrot_tilted_ring(MSs_v[i],MSs_PA_ring[i],MSs_PA[i], MSs_i[i])
	MS_si_vrot_tilt.append(vrot)

RG_si_vrot_tilt=[]
for i in range(len(RGs_r)):
	vrot=Vrot_tilted_ring(RGs_v[i],RGs_PA_ring[i],RGs_PA[i], RGs_i[i])
	RG_si_vrot_tilt.append(vrot)

AGo_si_vrot_tilt=[]
for i in range(len(AGos_r)):
	vrot=Vrot_tilted_ring(AGos_v[i],AGos_PA_ring[i],AGos_PA[i], AGos_i[i])
	AGo_si_vrot_tilt.append(vrot)

AGy_si_vrot_tilt=[]
for i in range(len(AGys_r)):
	vrot=Vrot_tilted_ring(AGys_v[i],AGys_PA_ring[i],AGys_PA[i], AGys_i[i])
	AGy_si_vrot_tilt.append(vrot)

#smoothed velocity and smoothed PA without tilted ring model
MS_ss_vrot=[]
for i in range(len(MSs_r)):
	vrot=Vrot(MSs_v[i],MSs_sPA[i])
	MS_ss_vrot.append(vrot)

RG_ss_vrot=[]
for i in range(len(RGs_r)):
	vrot=Vrot(RGs_v[i],RGs_sPA[i])
	RG_ss_vrot.append(vrot)

AGo_ss_vrot=[]
for i in range(len(AGos_r)):
	vrot=Vrot(AGos_v[i],AGos_sPA[i])
	AGo_ss_vrot.append(vrot)

AGy_ss_vrot=[]
for i in range(len(AGys_r)):
	vrot=Vrot(AGys_v[i],AGys_sPA[i])
	AGy_ss_vrot.append(vrot)

#smoothed v and smoothed PA with tilted ring model
MS_ss_vrot_tilt=[]
for i in range(len(MSs_r)):
	vrot=Vrot_tilted_ring(MSs_v[i],MSs_PA_ring[i],MSs_sPA[i], MSs_i[i])
	MS_ss_vrot_tilt.append(vrot)

RG_ss_vrot_tilt=[]
for i in range(len(RGs_r)):
	vrot=Vrot_tilted_ring(RGs_v[i],RGs_PA_ring[i],RGs_sPA[i], RGs_i[i])
	RG_ss_vrot_tilt.append(vrot)

AGo_ss_vrot_tilt=[]
for i in range(len(AGos_r)):
	vrot=Vrot_tilted_ring(AGos_v[i],AGos_PA_ring[i],AGos_sPA[i], AGos_i[i])
	AGo_ss_vrot_tilt.append(vrot)

AGy_ss_vrot_tilt=[]
for i in range(len(AGys_r)):
	vrot=Vrot_tilted_ring(AGys_v[i],AGys_PA_ring[i],AGys_sPA[i], AGys_i[i])
	AGy_ss_vrot_tilt.append(vrot)

#plotting rotation curves
#smoothed velocity and position without model
plt.scatter(MSs_r, MS_ss_vrot, s=1, c='b')
plt.scatter(HI_r, HI_v, s=5, c='k')
plt.xlim(4,20)
plt.ylim(100,300)
plt.savefig('/Users/amandaquirk/Desktop/MS_ss_vrot.png')
plt.close()

plt.scatter(RGs_r, RG_ss_vrot, s=1, c='r')
plt.scatter(HI_r, HI_v, s=5, c='k')
plt.xlim(4,20)
plt.ylim(100,300)
plt.savefig('/Users/amandaquirk/Desktop/RG_ss_vrot.png')
plt.close()

plt.scatter(AGos_r, AGo_ss_vrot, s=1, c='m')
plt.scatter(HI_r, HI_v, s=5, c='k')
plt.xlim(4,20)
plt.ylim(100,300)
plt.savefig('/Users/amandaquirk/Desktop/AGo_ss_vrot.png')
plt.close()

plt.scatter(AGys_r, AGy_ss_vrot, s=1, c='m')
plt.scatter(HI_r, HI_v, s=5, c='k')
plt.xlim(4,20)
plt.ylim(100,300)
plt.savefig('/Users/amandaquirk/Desktop/AGy_ss_vrot.png')
plt.close()

#smoothed velocity and positin with model
plt.scatter(MSs_r, MS_ss_vrot_tilt, s=1, c='b')
plt.scatter(HI_r, HI_v, s=5, c='k')
plt.xlim(4,20)
plt.ylim(100,300)
plt.savefig('/Users/amandaquirk/Desktop/MS_ss_vrot_tr.png')
plt.close()

plt.scatter(RGs_r, RG_ss_vrot_tilt, s=1, c='r')
plt.scatter(HI_r, HI_v, s=5, c='k')
plt.xlim(4,20)
plt.ylim(100,300)
plt.savefig('/Users/amandaquirk/Desktop/RG_ss_vrot_tr.png')
plt.close()

plt.scatter(AGos_r, AGo_ss_vrot_tilt, s=1, c='m')
plt.scatter(HI_r, HI_v, s=5, c='k')
plt.xlim(4,20)
plt.ylim(100,300)
plt.savefig('/Users/amandaquirk/Desktop/AGo_ss_vrot_tr.png')
plt.close()

plt.scatter(AGys_r, AGy_ss_vrot_tilt, s=1, c='m')
plt.scatter(HI_r, HI_v, s=5, c='k')
plt.xlim(4,20)
plt.ylim(100,300)
plt.savefig('/Users/amandaquirk/Desktop/AGy_ss_vrot_tr.png')
plt.close()

#smoothed velocity, individual position without model
plt.scatter(MSs_r, MS_si_vrot, s=1, c='b')
plt.scatter(HI_r, HI_v, s=5, c='k')
plt.xlim(4,20)
plt.ylim(100,300)
plt.savefig('/Users/amandaquirk/Desktop/MS_si_vrot.png')
plt.close()

plt.scatter(RGs_r, RG_si_vrot, s=1, c='r')
plt.scatter(HI_r, HI_v, s=5, c='k')
plt.xlim(4,20)
plt.ylim(100,300)
plt.savefig('/Users/amandaquirk/Desktop/RG_si_vrot.png')
plt.close()

plt.scatter(AGos_r, AGo_si_vrot, s=1, c='m')
plt.scatter(HI_r, HI_v, s=5, c='k')
plt.xlim(4,20)
plt.ylim(100,300)
plt.savefig('/Users/amandaquirk/Desktop/AGo_si_vrot.png')
plt.close()

plt.scatter(AGys_r, AGy_si_vrot, s=1, c='m')
plt.scatter(HI_r, HI_v, s=5, c='k')
plt.xlim(4,20)
plt.ylim(100,300)
plt.savefig('/Users/amandaquirk/Desktop/AGy_si_vrot.png')
plt.close()

#smoothed velocity, individual PA with model
plt.scatter(MSs_r, MS_si_vrot_tilt, s=1, c='b')
plt.scatter(HI_r, HI_v, s=5, c='k')
plt.xlim(4,20)
plt.ylim(100,300)
plt.savefig('/Users/amandaquirk/Desktop/MS_si_vrot_tr.png')
plt.close()

plt.scatter(RGs_r, RG_si_vrot_tilt, s=1, c='r')
plt.scatter(HI_r, HI_v, s=5, c='k')
plt.xlim(4,20)
plt.ylim(100,300)
plt.savefig('/Users/amandaquirk/Desktop/RG_si_vrot_tr.png')
plt.close()

plt.scatter(AGos_r, AGo_si_vrot_tilt, s=1, c='m')
plt.scatter(HI_r, HI_v, s=5, c='k')
plt.xlim(4,20)
plt.ylim(100,300)
plt.savefig('/Users/amandaquirk/Desktop/AGo_si_vrot_tr.png')
plt.close()

plt.scatter(AGys_r, AGy_si_vrot_tilt, s=1, c='m')
plt.scatter(HI_r, HI_v, s=5, c='k')
plt.xlim(4,20)
plt.ylim(100,300)
plt.savefig('/Users/amandaquirk/Desktop/AGy_si_vrot_tr.png')
plt.close()

#individual velocity and position without model
plt.scatter(MS_r, MS_ii_vrot, s=1, c='b')
plt.scatter(HI_r, HI_v, s=5, c='k')
plt.xlim(4,20)
plt.ylim(100,300)
plt.savefig('/Users/amandaquirk/Desktop/MS_ii_vrot.png')
plt.close()

plt.scatter(RG_r, RG_ii_vrot, s=1, c='r')
plt.scatter(HI_r, HI_v, s=5, c='k')
plt.xlim(4,20)
plt.ylim(100,300)
plt.savefig('/Users/amandaquirk/Desktop/RG_ii_vrot.png')
plt.close()

plt.scatter(AGo_r, AGo_ii_vrot, s=1, c='m')
plt.scatter(HI_r, HI_v, s=5, c='k')
plt.xlim(4,20)
plt.ylim(100,300)
plt.savefig('/Users/amandaquirk/Desktop/AGo_ii_vrot.png')
plt.close()

plt.scatter(AGy_r, AGy_ii_vrot, s=1, c='m')
plt.scatter(HI_r, HI_v, s=5, c='k')
plt.xlim(4,20)
plt.ylim(100,300)
plt.savefig('/Users/amandaquirk/Desktop/AGy_ii_vrot.png')
plt.close()

#individual velocity and position with model
plt.scatter(MS_r, MS_ii_vrot_tilt, s=1, c='b')
plt.scatter(HI_r, HI_v, s=5, c='k')
plt.xlim(4,20)
plt.ylim(100,300)
plt.savefig('/Users/amandaquirk/Desktop/MS_ii_vrot_tr.png')
plt.close()

plt.scatter(RG_r, RG_ii_vrot_tilt, s=1, c='r')
plt.scatter(HI_r, HI_v, s=5, c='k')
plt.xlim(4,20)
plt.ylim(100,300)
plt.savefig('/Users/amandaquirk/Desktop/RG_ii_vrot_tr.png')
plt.close()

plt.scatter(AGo_r, AGo_ii_vrot_tilt, s=1, c='m')
plt.scatter(HI_r, HI_v, s=5, c='k')
plt.xlim(4,20)
plt.ylim(100,300)
plt.savefig('/Users/amandaquirk/Desktop/AGo_ii_vrot_tr.png')
plt.close()

plt.scatter(AGy_r, AGy_ii_vrot_tilt, s=1, c='m')
plt.scatter(HI_r, HI_v, s=5, c='k')
plt.xlim(4,20)
plt.ylim(100,300)
plt.savefig('/Users/amandaquirk/Desktop/AGy_ii_vrot_tr.png')
plt.close()

#files
#individual without model
file=open('/Users/amandaquirk/Desktop/MS_ii_vrot.txt', 'w')
file.write('#r (kpc), v_rot (km/s)\n')
for i in range(len(MS_r)):
	file.write('{} {}\n'.format(MS_r[i],MS_ii_vrot[i]))
file.close()

file=open('/Users/amandaquirk/Desktop/AGy_ii_vrot.txt', 'w')
file.write('#r (kpc), v_rot (km/s)\n')
for i in range(len(AGy_r)):
	file.write('{} {}\n'.format(AGy_r[i],AGy_ii_vrot[i]))
file.close()

file=open('/Users/amandaquirk/Desktop/AGo_ii_vrot.txt', 'w')
file.write('#r (kpc), v_rot (km/s)\n')
for i in range(len(AGo_r)):
	file.write('{} {}\n'.format(AGo_r[i],AGo_ii_vrot[i]))
file.close()

file=open('/Users/amandaquirk/Desktop/RG_ii_vrot.txt', 'w')
file.write('#r (kpc), v_rot (km/s)\n')
for i in range(len(RG_r)):
	file.write('{} {}\n'.format(RG_r[i],RG_ii_vrot[i]))
file.close()

#individual with model
file=open('/Users/amandaquirk/Desktop/MS_ii_vrot_tr.txt', 'w')
file.write('#r (kpc), v_rot (km/s)\n')
for i in range(len(MS_r)):
	file.write('{} {}\n'.format(MS_r[i],MS_ii_vrot_tilt[i]))
file.close()

file=open('/Users/amandaquirk/Desktop/AGy_ii_vrot_tr.txt', 'w')
file.write('#r (kpc), v_rot (km/s)\n')
for i in range(len(AGy_r)):
	file.write('{} {}\n'.format(AGy_r[i],AGy_ii_vrot_tilt[i]))
file.close()

file=open('/Users/amandaquirk/Desktop/AGo_ii_vrot_tr.txt', 'w')
file.write('#r (kpc), v_rot (km/s)\n')
for i in range(len(AGo_r)):
	file.write('{} {}\n'.format(AGo_r[i],AGo_ii_vrot_tilt[i]))
file.close()

file=open('/Users/amandaquirk/Desktop/RG_ii_vrot_tr.txt', 'w')
file.write('#r (kpc), v_rot (km/s)\n')
for i in range(len(RG_r)):
	file.write('{} {}\n'.format(RG_r[i],RG_ii_vrot_tilt[i]))
file.close()

#smoothed v, individual PA without model
file=open('/Users/amandaquirk/Desktop/MS_si_vrot.txt', 'w')
file.write('#r (kpc), v_rot (km/s)\n')
for i in range(len(MSs_r)):
	file.write('{} {}\n'.format(MSs_r[i],MS_si_vrot[i]))
file.close()

file=open('/Users/amandaquirk/Desktop/AGy_si_vrot.txt', 'w')
file.write('#r (kpc), v_rot (km/s)\n')
for i in range(len(AGys_r)):
	file.write('{} {}\n'.format(AGys_r[i],AGy_si_vrot[i]))
file.close()

file=open('/Users/amandaquirk/Desktop/AGo_si_vrot.txt', 'w')
file.write('#r (kpc), v_rot (km/s)\n')
for i in range(len(AGos_r)):
	file.write('{} {}\n'.format(AGos_r[i],AGo_si_vrot[i]))
file.close()

file=open('/Users/amandaquirk/Desktop/RG_si_vrot.txt', 'w')
file.write('#r (kpc), v_rot (km/s)\n')
for i in range(len(RGs_r)):
	file.write('{} {}\n'.format(RGs_r[i],RG_si_vrot[i]))
file.close()

#smoothed v, individual PA with model
file=open('/Users/amandaquirk/Desktop/MS_si_vrot_tr.txt', 'w')
file.write('#r (kpc), v_rot (km/s)\n')
for i in range(len(MSs_r)):
	file.write('{} {}\n'.format(MSs_r[i],MS_si_vrot_tilt[i]))
file.close()

file=open('/Users/amandaquirk/Desktop/AGy_si_vrot_tr.txt', 'w')
file.write('#r (kpc), v_rot (km/s)\n')
for i in range(len(AGys_r)):
	file.write('{} {}\n'.format(AGys_r[i],AGy_si_vrot_tilt[i]))
file.close()

file=open('/Users/amandaquirk/Desktop/AGo_si_vrot_tr.txt', 'w')
file.write('#r (kpc), v_rot (km/s)\n')
for i in range(len(AGos_r)):
	file.write('{} {}\n'.format(AGos_r[i],AGo_si_vrot_tilt[i]))
file.close()

file=open('/Users/amandaquirk/Desktop/RG_si_vrot_tr.txt', 'w')
file.write('#r (kpc), v_rot (km/s)\n')
for i in range(len(RGs_r)):
	file.write('{} {}\n'.format(RGs_r[i],RG_si_vrot_tilt[i]))
file.close()

#smoothed v and PA without model
file=open('/Users/amandaquirk/Desktop/MS_ss_vrot.txt', 'w')
file.write('#r (kpc), v_rot (km/s)\n')
for i in range(len(MSs_r)):
	file.write('{} {}\n'.format(MSs_r[i],MS_ss_vrot[i]))
file.close()

file=open('/Users/amandaquirk/Desktop/AGy_ss_vrot.txt', 'w')
file.write('#r (kpc), v_rot (km/s)\n')
for i in range(len(AGys_r)):
	file.write('{} {}\n'.format(AGys_r[i],AGy_ss_vrot[i]))
file.close()

file=open('/Users/amandaquirk/Desktop/AGo_ss_vrot.txt', 'w')
file.write('#r (kpc), v_rot (km/s)\n')
for i in range(len(AGos_r)):
	file.write('{} {}\n'.format(AGos_r[i],AGo_ss_vrot[i]))
file.close()

file=open('/Users/amandaquirk/Desktop/RG_ss_vrot.txt', 'w')
file.write('#r (kpc), v_rot (km/s)\n')
for i in range(len(RGs_r)):
	file.write('{} {}\n'.format(RGs_r[i],RG_ss_vrot[i]))
file.close()

#smoothed v and individual PA with model
file=open('/Users/amandaquirk/Desktop/MS_ss_vrot_tr.txt', 'w')
file.write('#r (kpc), v_rot (km/s)\n')
for i in range(len(MSs_r)):
	file.write('{} {}\n'.format(MSs_r[i],MS_ss_vrot_tilt[i]))
file.close()

file=open('/Users/amandaquirk/Desktop/AGy_ss_vrot_tr.txt', 'w')
file.write('#r (kpc), v_rot (km/s)\n')
for i in range(len(AGys_r)):
	file.write('{} {}\n'.format(AGys_r[i],AGy_ss_vrot_tilt[i]))
file.close()

file=open('/Users/amandaquirk/Desktop/AGo_ss_vrot_tr.txt', 'w')
file.write('#r (kpc), v_rot (km/s)\n')
for i in range(len(AGos_r)):
	file.write('{} {}\n'.format(AGos_r[i],AGo_ss_vrot_tilt[i]))
file.close()

file=open('/Users/amandaquirk/Desktop/RG_ss_vrot_tr.txt', 'w')
file.write('#r (kpc), v_rot (km/s)\n')
for i in range(len(RGs_r)):
	file.write('{} {}\n'.format(RGs_r[i],RG_ss_vrot_tilt[i]))
file.close()
#now i calculate the asymmetric drift 
def asymmetric_drift(v_HI, v_star):
	return v_HI-v_star

# #smoothed position, smoothed v, tilted ring model
# MS_ad=[]
# for i in range(len(MSs_HIv)):
# 	MS_ad.append(asymmetric_drift(MSs_HIv[i], MS_ss_vrot_tilt[i]))

# RG_ad=[]
# for i in range(len(RGs_HIv)):
# 	RG_ad.append(asymmetric_drift(RGs_HIv[i], RG_ss_vrot_tilt[i]))

# AGy_ad=[]
# for i in range(len(AGys_HIv)):
# 	AGy_ad.append(asymmetric_drift(AGys_HIv[i], AGy_ss_vrot_tilt[i]))

# AGo_ad=[]
# for i in range(len(AGos_HIv)):
# 	AGo_ad.append(asymmetric_drift(AGos_HIv[i], AGo_ss_vrot_tilt[i]))

# plt.hist(MS_ad)
# plt.savefig('/Users/amandaquirk/Desktop/MS_ss_tilt_ad.png')
# plt.close()

# plt.hist(RG_ad)
# plt.savefig('/Users/amandaquirk/Desktop/RG_ss_tilt_ad.png')
# plt.close()

# plt.hist(AGy_ad)
# plt.savefig('/Users/amandaquirk/Desktop/AGy_ss_tilt_ad.png')
# plt.close()

# plt.hist(AGo_ad)
# plt.savefig('/Users/amandaquirk/Desktop/AGo_ss_tilt_ad.png')
# plt.close()




