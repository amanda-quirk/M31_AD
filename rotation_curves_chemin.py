import numpy as np 
import matplotlib.pyplot as plt 
from astropy import units as u
from astropy.coordinates import SkyCoord

#importing all of the data
#xi (kpc), eta (kpc), average v(km/s), v err,var, n, HI main, HI close <-- header of data file to be read in 
MS_xi, MS_eta, MS_v, MS_n, MS_HImain, MS_HIclose=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/MS_smoothed_chemin.txt', usecols=(0,1,2,5,6,7,), unpack=True)
AGy_xi, AGy_eta, AGy_v, AGy_n, AGy_HImain, AGy_HIclose=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/AGy_smoothed_chemin.txt', usecols=(0,1,2,5,6,7,), unpack=True)
AGo_xi, AGo_eta, AGo_v, AGo_n, AGo_HImain, AGo_HIclose=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/AGo_smoothed_chemin.txt', usecols=(0,1,2,5,6,7,), unpack=True)
RG_xi, RG_eta, RG_v, RG_n, RG_HImain, RG_HIclose=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/RG_smoothed_chemin.txt', usecols=(0,1,2,5,6,7,), unpack=True)

#below defines a function that converts eta and xi and then uses the shifted coordinates to find the deprojected radius. to calculate the deprojected radius from the center, we must shift the xi and eta according to the PA of the disk. (This is a simple shift of coordinate axes.) We will assume the disk has a PA of 37 deg and an inclination of 77 deg 
def x(xi, eta):
	xi_deg=xi/13.67
	eta_deg=eta/13.67
	sine=np.sin(float(37*np.pi) / 180)
	cosine=np.cos(float(37*np.pi) / 180)
	x=(xi_deg*cosine)-(eta_deg*sine)
	return x 

def y(xi, eta):
	xi_deg=xi/13.67
	eta_deg=eta/13.67
	sine=np.sin(float(37*np.pi) / 180)
	cosine=np.cos(float(37*np.pi) / 180)
	y=(eta_deg*cosine)+(xi_deg*sine)
	return y

#getting coordinates in shifted coordinate system (centered on M31)
MS_x=x(MS_xi, MS_eta) 
MS_y=y(MS_xi, MS_eta)

AGy_x=x(AGy_xi, AGy_eta) 
AGy_y=y(AGy_xi, AGy_eta)

AGo_x=x(AGo_xi, AGo_eta) 
AGo_y=y(AGo_xi, AGo_eta)

RG_x=x(RG_xi, RG_eta) 
RG_y=y(RG_xi, RG_eta)

def distance(x, y):
	inclination_factor=np.cos(float(77*np.pi) / 180)**2
	ang_distance_sq=[(a**2)+(float(b**2)/inclination_factor) for a,b in zip(y,x)]
	ang_dist=[np.sqrt(a) for a in ang_distance_sq]
	dist=[a*13.67 for a in ang_dist]
	return dist

#below will contain the deprojected distances of each star and gas cloud
MS_r=distance(MS_x,MS_y)
AGy_r=distance(AGy_x,AGy_y)
AGo_r=distance(AGo_x,AGo_y)
RG_r=distance(RG_x,RG_y)

#the function below will assign each star and gas cloud a position angle using their xi and eta coordinates
def PA(x,y): #defined as west of north as positioned in maps (inverted x axis)
	if x>0:
		rad=np.arctan(float(y)/x)
		deg=90-(float(rad*180)/np.pi)
	else:
		rad=np.arctan(float(y)/x)
		deg=270-(float(rad*180)/np.pi)
	return deg + 37

MS_PA=[]
for i in range(len(MS_r)):
	MS_PA.append(PA(MS_x[i], MS_y[i]))
AGy_PA=[]
for i in range(len(AGy_r)):
	AGy_PA.append(PA(AGy_x[i], AGy_y[i]))
AGo_PA=[]
for i in range(len(AGo_r)):
	AGo_PA.append(PA(AGo_x[i], AGo_y[i]))
RG_PA=[]
for i in range(len(RG_r)):
	RG_PA.append(PA(RG_x[i], RG_y[i]))

#now we need to assign each star a PA_ring and i_ring using HI data. we define these values by seeing which HI ring a star is cloest to and assigning it the corresponding PA_ring and i_ring value (HI_PA and HI_i respectively). we also will pair a HI velocity to each star to later calculate the asymmetric drift
#reads in the HI data
HI_r, HI_PA, HI_i, HI_v=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/HI_PA_i_vrot.txt', unpack=True)
HI_PA_avg = np.mean(HI_PA)
HI_i_avg = np.mean(HI_i)

#below defines a function to determine which HI ring a star is cloest to
# def find_nearest_ring(radius):
# 	idx=[]
# 	for i in range(len(HI_r)):
# 		idx.append(np.abs(HI_r[i]-radius))
# 		n=np.argmin(idx)
# 	return n #index of the radius closest to the star's radius

# #assigning a PA  and i to each star and HI component
# MS_PA_ring=[]
# MS_i=[]
# for j in range(len(MS_r)):
# 	N=find_nearest_ring(MS_r[j])
# 	MS_PA_ring.append(HI_PA[N])
# 	MS_i.append(HI_i[N])

# RG_PA_ring=[]
# RG_i=[]
# for j in range(len(RG_r)):
# 	N=find_nearest_ring(RG_r[j])
# 	RG_PA_ring.append(HI_PA[N])
# 	RG_i.append(HI_i[N])

# AGy_PA_ring=[]
# AGy_i=[]
# for j in range(len(AGy_r)):
# 	N=find_nearest_ring(AGy_r[j])
# 	AGy_PA_ring.append(HI_PA[N])
# 	AGy_i.append(HI_i[N])

# AGo_PA_ring=[]
# AGo_i=[]
# for j in range(len(AGo_r)):
# 	N=find_nearest_ring(AGo_r[j])
# 	AGo_PA_ring.append(HI_PA[N])
# 	AGo_i.append(HI_i[N])

# #now we need to calculate the rotation speed of the stars using the PA and LOS velocity. I will use the equation in Claire's thesis on page 203 (but the corrected version I have in my evernote) using the PA and i value given by Table 5 in Chemin et al. 2009
# #I define a functionto calculate the rotation speed without a tiled ring model
# def Vrot(v, PA_star):
# 	sine=np.sin(float(77*np.pi) / 180)
# 	inclination_factor=(np.cos(float(77*np.pi) / 180))**2
# 	vsys= -300 #km/s, as defined in Claire's thesis
# 	A=[float(a-vsys)/sine for a in v]
# 	B= [(np.tan(float((37-b)*np.pi) / 180))**2 for b in PA_star]
# 	rotation_velocity= [a*np.sqrt(1+(float(b)/inclination_factor))for a,b in zip(A,B)]
# 	positive=[np.absolute(a) for a in rotation_velocity]
# 	return positive

#I define a function to caluclate the rotation speed using the tiled ring model
def Vrot_tilted_ring(v,PA_ring,PA_star, i_ring): 
	vsys= -300 #km/s, as defined in Claire's thesis
	A=[float(a-vsys)/np.sin(float(b*np.pi) / 180) for a, b in zip(v, i_ring)]
	B= [(np.tan(float((a-b)*np.pi) / 180))**2 for a,b in zip(PA_ring, PA_star)]
	C= [(np.cos(float(a*np.pi) / 180))**2 for a in i_ring]
	rotation_velocity= [a*np.sqrt(1+(float(b)/c)) for a,b,c in zip(A,B,C)]
	positive=[np.absolute(a) for a in rotation_velocity]
	return positive

#caluclating the rotation velocity using the above defined function
# #smoothed v and individual PA without tilted ring model
# MS_si_vrot=Vrot(MS_v,MS_PA)
# RG_si_vrot=Vrot(RG_v,RG_PA)
# AGo_si_vrot=Vrot(AGo_v,AGo_PA)
# AGy_si_vrot=Vrot(AGy_v,AGy_PA)

# # #smoothed v and individual PA with tilted ring model
# MS_si_vrot_tilt=Vrot_tilted_ring(MS_v,MS_PA_ring,MS_PA, MS_i)
# RG_si_vrot_tilt=Vrot_tilted_ring(RG_v,RG_PA_ring,RG_PA, RG_i)
# AGo_si_vrot_tilt=Vrot_tilted_ring(AGo_v,AGo_PA_ring,AGo_PA, AGo_i)
# AGy_si_vrot_tilt=Vrot_tilted_ring(AGy_v,AGy_PA_ring,AGy_PA, AGy_i)

# #i want to also smooth the PA in the same way velocity has been smoothed
# MS_sPA=[] #will contain the smoothed PA
# for i in range(len(MS_PA)):
# 	c1=SkyCoord(MS_xi[i], MS_eta[i], unit=(u.deg, u.deg)) #our first point
# 	positions=[]
# 	for j in range(len(MS_PA)):
# 		c2=SkyCoord(MS_xi[j], MS_eta[j], unit=(u.deg, u.deg)) #second point
# 		sep=c1.separation(c2) #gets the separation of the two points
# 		if sep.arcsecond<200: #if the two points are close enough, the position is added to an array to later be averaged
# 			positions.append(MS_PA[j])
# 	MS_sPA.append(np.mean(positions))

# AGy_sPA=[] #will contain the smoothed PA
# for i in range(len(AGy_PA)):
# 	c1=SkyCoord(AGy_xi[i], AGy_eta[i], unit=(u.deg, u.deg)) #our first point
# 	positions=[]
# 	for j in range(len(AGy_PA)):
# 		c2=SkyCoord(AGy_xi[j], AGy_eta[j], unit=(u.deg, u.deg)) #second point
# 		sep=c1.separation(c2) #gets the separation of the two points
# 		if sep.arcsecond<275: #if the two points are close enough, the position is added to an array to later be averaged
# 			positions.append(AGy_PA[j])
# 	AGy_sPA.append(np.mean(positions))

# AGo_sPA=[] #will contain the smoothed PA
# for i in range(len(AGo_PA)):
# 	c1=SkyCoord(AGo_xi[i], AGo_eta[i], unit=(u.deg, u.deg)) #our first point
# 	positions=[]
# 	for j in range(len(AGo_PA)):
# 		c2=SkyCoord(AGo_xi[j], AGo_eta[j], unit=(u.deg, u.deg)) #second point
# 		sep=c1.separation(c2) #gets the separation of the two points
# 		if sep.arcsecond<275: #if the two points are close enough, the position is added to an array to later be averaged
# 			positions.append(AGo_PA[j])
# 	AGo_sPA.append(np.mean(positions))

# RG_sPA=[] #will contain the smoothed PA
# for i in range(len(RG_PA)):
# 	c1=SkyCoord(RG_xi[i], RG_eta[i], unit=(u.deg, u.deg)) #our first point
# 	positions=[]
# 	for j in range(len(RG_PA)):
# 		c2=SkyCoord(RG_xi[j], RG_eta[j], unit=(u.deg, u.deg)) #second point
# 		sep=c1.separation(c2) #gets the separation of the two points
# 		if sep.arcsecond<200: #if the two points are close enough, the position is added to an array to later be averaged
# 			positions.append(RG_PA[j])
# 	RG_sPA.append(np.mean(positions))

# #smoothed v and PA without tilted ring model
# MS_ss_vrot=Vrot(MS_v,MS_sPA)
# RG_ss_vrot=Vrot(RG_v,RG_sPA)
# AGo_ss_vrot=Vrot(AGo_v,AGo_sPA)
# AGy_ss_vrot=Vrot(AGy_v,AGy_sPA)

# #smoothed v and PA with tilted ring model
# MS_ss_vrot_tilt=Vrot_tilted_ring(MS_v,MS_PA_ring,MS_sPA, MS_i)
# RG_ss_vrot_tilt=Vrot_tilted_ring(RG_v,RG_PA_ring,RG_sPA, RG_i)
# AGo_ss_vrot_tilt=Vrot_tilted_ring(AGo_v,AGo_PA_ring,AGo_sPA, AGo_i)
# AGy_ss_vrot_tilt=Vrot_tilted_ring(AGy_v,AGy_PA_ring,AGy_sPA, AGy_i)

# # #HI without model
# MS_HImain_vrot=Vrot(MS_HImain, MS_PA) #not using smoothed PA
# MS_HIclose_vrot=Vrot(MS_HIclose,MS_PA)
# AGy_HImain_vrot=Vrot(AGy_HImain, AGy_PA) #not using smoothed PA
# AGy_HIclose_vrot=Vrot(AGy_HIclose,AGy_PA)
# AGo_HImain_vrot=Vrot(AGo_HImain, AGo_PA) #not using smoothed PA
# AGo_HIclose_vrot=Vrot(AGo_HIclose,AGo_PA)
# RG_HImain_vrot=Vrot(RG_HImain, RG_PA) #not using smoothed PA
# RG_HIclose_vrot=Vrot(RG_HIclose,RG_PA)

# # #HI with model
# MS_HImain_vrot_tilt=Vrot_tilted_ring(MS_HImain,MS_PA_ring,MS_PA, MS_i)
# MS_HIclose_vrot_tilt=Vrot_tilted_ring(MS_HIclose,MS_PA_ring,MS_PA, MS_i)
# AGy_HImain_vrot_tilt=Vrot_tilted_ring(AGy_HImain,AGy_PA_ring,AGy_PA, AGy_i)
# AGy_HIclose_vrot_tilt=Vrot_tilted_ring(AGy_HIclose,AGy_PA_ring,AGy_PA, AGy_i)
# AGo_HImain_vrot_tilt=Vrot_tilted_ring(AGo_HImain,AGo_PA_ring,AGo_PA, AGo_i)
# AGo_HIclose_vrot_tilt=Vrot_tilted_ring(AGo_HIclose,AGo_PA_ring,AGo_PA, AGo_i)
# RG_HImain_vrot_tilt=Vrot_tilted_ring(RG_HImain,RG_PA_ring,RG_PA, RG_i)
# RG_HIclose_vrot_tilt=Vrot_tilted_ring(RG_HIclose,RG_PA_ring,RG_PA, RG_i)

# #plotting
# #smooth v and individual PA without model
# plt.scatter(MS_r, MS_si_vrot, s=1, c='b')
# plt.scatter(MS_r, MS_HImain_vrot, s=1, c='k')
# plt.xlim(4,20)
# plt.ylim(100,300)
# plt.xlabel('r [kpc]')
# plt.ylabel('rotational v [km/s]')
# plt.savefig('/Users/amandaquirk/Desktop/MS_si_vrot_HImain.png')
# plt.close()

# plt.scatter(RG_r, RG_si_vrot, s=1, c='r')
# plt.scatter(RG_r, RG_HImain_vrot, s=1, c='k')
# plt.xlim(4,20)
# plt.ylim(100,300)
# plt.xlabel('r [kpc]')
# plt.ylabel('rotational v [km/s]')
# plt.savefig('/Users/amandaquirk/Desktop/RG_si_vrot_HImain.png')
# plt.close()

# plt.scatter(AGo_r, AGo_si_vrot, s=1, c='m')
# plt.scatter(AGo_r, AGo_HImain_vrot, s=1, c='k')
# plt.xlim(4,20)
# plt.ylim(100,300)
# plt.xlabel('r [kpc]')
# plt.ylabel('rotational v [km/s]')
# plt.savefig('/Users/amandaquirk/Desktop/AGo_si_vrot_HImain.png')
# plt.close()

# plt.scatter(AGy_r, AGy_si_vrot, s=1, c='m')
# plt.scatter(AGy_r, AGy_HImain_vrot, s=1, c='k')
# plt.xlim(4,20)
# plt.ylim(100,300)
# plt.xlabel('r [kpc]')
# plt.ylabel('rotational v [km/s]')
# plt.savefig('/Users/amandaquirk/Desktop/AGy_si_vrot_HImain.png')
# plt.close()

# plt.scatter(MS_r, MS_si_vrot, s=1, c='b')
# plt.scatter(MS_r, MS_HIclose_vrot, s=1, c='k')
# plt.xlim(4,20)
# plt.ylim(100,300)
# plt.xlabel('r [kpc]')
# plt.ylabel('rotational v [km/s]')
# plt.savefig('/Users/amandaquirk/Desktop/MS_si_vrot_HIclose.png')
# plt.close()

# plt.scatter(RG_r, RG_si_vrot, s=1, c='r')
# plt.scatter(RG_r, RG_HIclose_vrot, s=1, c='k')
# plt.xlim(4,20)
# plt.ylim(100,300)
# plt.xlabel('r [kpc]')
# plt.ylabel('rotational v [km/s]')
# plt.savefig('/Users/amandaquirk/Desktop/RG_si_vrot_HIclose.png')
# plt.close()

# plt.scatter(AGo_r, AGo_si_vrot, s=1, c='m')
# plt.scatter(AGo_r, AGo_HIclose_vrot, s=1, c='k')
# plt.xlim(4,20)
# plt.ylim(100,300)
# plt.xlabel('r [kpc]')
# plt.ylabel('rotational v [km/s]')
# plt.savefig('/Users/amandaquirk/Desktop/AGo_si_vrot_HIclose.png')
# plt.close()

# plt.scatter(AGy_r, AGy_si_vrot, s=1, c='m')
# plt.scatter(AGy_r, AGy_HIclose_vrot, s=1, c='k')
# plt.xlim(4,20)
# plt.ylim(100,300)
# plt.xlabel('r [kpc]')
# plt.ylabel('rotational v [km/s]')
# plt.savefig('/Users/amandaquirk/Desktop/AGy_si_vrot_HIclose.png')
# plt.close()

# #smooth v and PA without model
# plt.scatter(MS_r, MS_ss_vrot, s=1, c='b')
# plt.scatter(MS_r, MS_HImain_vrot, s=1, c='k')
# plt.xlim(4,20)
# plt.ylim(100,300)
# plt.xlabel('r [kpc]')
# plt.ylabel('rotational v [km/s]')
# plt.savefig('/Users/amandaquirk/Desktop/MS_ss_vrot_HImain.png')
# plt.close()

# plt.scatter(RG_r, RG_ss_vrot, s=1, c='r')
# plt.scatter(RG_r, RG_HImain_vrot, s=1, c='k')
# plt.xlim(4,20)
# plt.ylim(100,300)
# plt.xlabel('r [kpc]')
# plt.ylabel('rotational v [km/s]')
# plt.savefig('/Users/amandaquirk/Desktop/RG_ss_vrot_HImain.png')
# plt.close()

# plt.scatter(AGo_r, AGo_ss_vrot, s=1, c='m')
# plt.scatter(AGo_r, AGo_HImain_vrot, s=1, c='k')
# plt.xlim(4,20)
# plt.ylim(100,300)
# plt.xlabel('r [kpc]')
# plt.ylabel('rotational v [km/s]')
# plt.savefig('/Users/amandaquirk/Desktop/AGo_ss_vrot_HImain.png')
# plt.close()

# plt.scatter(AGy_r, AGy_ss_vrot, s=1, c='m')
# plt.scatter(AGy_r, AGy_HImain_vrot, s=1, c='k')
# plt.xlim(4,20)
# plt.ylim(100,300)
# plt.xlabel('r [kpc]')
# plt.ylabel('rotational v [km/s]')
# plt.savefig('/Users/amandaquirk/Desktop/AGy_ss_vrot_HImain.png')
# plt.close()

# plt.scatter(MS_r, MS_ss_vrot, s=1, c='b')
# plt.scatter(MS_r, MS_HIclose_vrot, s=1, c='k')
# plt.xlim(4,20)
# plt.ylim(100,300)
# plt.xlabel('r [kpc]')
# plt.ylabel('rotational v [km/s]')
# plt.savefig('/Users/amandaquirk/Desktop/MS_ss_vrot_HIclose.png')
# plt.close()

# plt.scatter(RG_r, RG_ss_vrot, s=1, c='r')
# plt.scatter(RG_r, RG_HIclose_vrot, s=1, c='k')
# plt.xlim(4,20)
# plt.ylim(100,300)
# plt.xlabel('r [kpc]')
# plt.ylabel('rotational v [km/s]')
# plt.savefig('/Users/amandaquirk/Desktop/RG_ss_vrot_HIclose.png')
# plt.close()

# plt.scatter(AGo_r, AGo_ss_vrot, s=1, c='m')
# plt.scatter(AGo_r, AGo_HIclose_vrot, s=1, c='k')
# plt.xlim(4,20)
# plt.ylim(100,300)
# plt.xlabel('r [kpc]')
# plt.ylabel('rotational v [km/s]')
# plt.savefig('/Users/amandaquirk/Desktop/AGo_ss_vrot_HIclose.png')
# plt.close()

# plt.scatter(AGy_r, AGy_ss_vrot, s=1, c='m')
# plt.scatter(AGy_r, AGy_HIclose_vrot, s=1, c='k')
# plt.xlim(4,20)
# plt.ylim(100,300)
# plt.xlabel('r [kpc]')
# plt.ylabel('rotational v [km/s]')
# plt.savefig('/Users/amandaquirk/Desktop/AGy_ss_vrot_HIclose.png')
# plt.close()

# #smooth v and individual PA with model
# plt.scatter(MS_r, MS_si_vrot_tilt, s=1, c='b')
# plt.scatter(MS_r, MS_HImain_vrot_tilt, s=1, c='k')
# plt.xlim(4,20)
# plt.ylim(100,300)
# plt.xlabel('r [kpc]')
# plt.ylabel('rotational v [km/s]')
# plt.savefig('/Users/amandaquirk/Desktop/MS_si_vrot_HImain_tilt.png')
# plt.close()

# plt.scatter(RG_r, RG_si_vrot_tilt, s=1, c='r')
# plt.scatter(RG_r, RG_HImain_vrot_tilt, s=1, c='k')
# plt.xlim(4,20)
# plt.ylim(100,300)
# plt.xlabel('r [kpc]')
# plt.ylabel('rotational v [km/s]')
# plt.savefig('/Users/amandaquirk/Desktop/RG_si_vrot_HImain_tilt.png')
# plt.close()

# plt.scatter(AGo_r, AGo_si_vrot_tilt, s=1, c='m')
# plt.scatter(AGo_r, AGo_HImain_vrot_tilt, s=1, c='k')
# plt.xlim(4,20)
# plt.ylim(100,300)
# plt.xlabel('r [kpc]')
# plt.ylabel('rotational v [km/s]')
# plt.savefig('/Users/amandaquirk/Desktop/AGo_si_vrot_HImain_tilt.png')
# plt.close()

# plt.scatter(AGy_r, AGy_si_vrot_tilt, s=1, c='m')
# plt.scatter(AGy_r, AGy_HImain_vrot_tilt, s=1, c='k')
# plt.xlim(4,20)
# plt.ylim(100,300)
# plt.xlabel('r [kpc]')
# plt.ylabel('rotational v [km/s]')
# plt.savefig('/Users/amandaquirk/Desktop/AGy_si_vrot_HImain_tilt.png')
# plt.close()

# plt.scatter(MS_r, MS_si_vrot_tilt, s=1, c='b')
# plt.scatter(MS_r, MS_HIclose_vrot_tilt, s=1, c='k')
# plt.xlim(4,20)
# plt.ylim(100,300)
# plt.xlabel('r [kpc]')
# plt.ylabel('rotational v [km/s]')
# plt.savefig('/Users/amandaquirk/Desktop/MS_si_vrot_HIclose_tilt.png')
# plt.close()

# plt.scatter(RG_r, RG_si_vrot_tilt, s=1, c='r')
# plt.scatter(RG_r, RG_HIclose_vrot_tilt, s=1, c='k')
# plt.xlim(4,20)
# plt.ylim(100,300)
# plt.xlabel('r [kpc]')
# plt.ylabel('rotational v [km/s]')
# plt.savefig('/Users/amandaquirk/Desktop/RG_si_vrot_HIclose_tilt.png')
# plt.close()

# plt.scatter(AGo_r, AGo_si_vrot_tilt, s=1, c='m')
# plt.scatter(AGo_r, AGo_HIclose_vrot_tilt, s=1, c='k')
# plt.xlim(4,20)
# plt.ylim(100,300)
# plt.xlabel('r [kpc]')
# plt.ylabel('rotational v [km/s]')
# plt.savefig('/Users/amandaquirk/Desktop/AGo_si_vrot_HIclose_tilt.png')
# plt.close()

# plt.scatter(AGy_r, AGy_si_vrot_tilt, s=1, c='m')
# plt.scatter(AGy_r, AGy_HIclose_vrot_tilt, s=1, c='k')
# plt.xlim(4,20)
# plt.ylim(100,300)
# plt.xlabel('r [kpc]')
# plt.ylabel('rotational v [km/s]')
# plt.savefig('/Users/amandaquirk/Desktop/AGy_si_vrot_HIclose_tilt.png')
# plt.close()

# #smooth v and PA without model
# plt.scatter(MS_r, MS_ss_vrot_tilt, s=1, c='b')
# plt.scatter(MS_r, MS_HImain_vrot_tilt, s=1, c='k')
# plt.xlim(4,20)
# plt.ylim(100,300)
# plt.xlabel('r [kpc]')
# plt.ylabel('rotational v [km/s]')
# plt.savefig('/Users/amandaquirk/Desktop/MS_ss_vrot_HImain_tilt.png')
# plt.close()

# plt.scatter(RG_r, RG_ss_vrot_tilt, s=1, c='r')
# plt.scatter(RG_r, RG_HImain_vrot_tilt, s=1, c='k')
# plt.xlim(4,20)
# plt.ylim(100,300)
# plt.xlabel('r [kpc]')
# plt.ylabel('rotational v [km/s]')
# plt.savefig('/Users/amandaquirk/Desktop/RG_ss_vrot_HImain_tilt.png')
# plt.close()

# plt.scatter(AGo_r, AGo_ss_vrot_tilt, s=1, c='m')
# plt.scatter(AGo_r, AGo_HImain_vrot_tilt, s=1, c='k')
# plt.xlim(4,20)
# plt.ylim(100,300)
# plt.xlabel('r [kpc]')
# plt.ylabel('rotational v [km/s]')
# plt.savefig('/Users/amandaquirk/Desktop/AGo_ss_vrot_HImain_tilt.png')
# plt.close()

# plt.scatter(AGy_r, AGy_ss_vrot_tilt, s=1, c='m')
# plt.scatter(AGy_r, AGy_HImain_vrot_tilt, s=1, c='k')
# plt.xlim(4,20)
# plt.ylim(100,300)
# plt.xlabel('r [kpc]')
# plt.ylabel('rotational v [km/s]')
# plt.savefig('/Users/amandaquirk/Desktop/AGy_ss_vrot_HImain_tilt.png')
# plt.close()

# plt.scatter(MS_r, MS_ss_vrot_tilt, s=1, c='b')
# plt.scatter(MS_r, MS_HIclose_vrot_tilt, s=1, c='k')
# plt.xlim(4,20)
# plt.ylim(100,300)
# plt.xlabel('r [kpc]')
# plt.ylabel('rotational v [km/s]')
# plt.savefig('/Users/amandaquirk/Desktop/MS_ss_vrot_HIclose_tilt.png')
# plt.close()

# plt.scatter(RG_r, RG_ss_vrot_tilt, s=1, c='r')
# plt.scatter(RG_r, RG_HIclose_vrot_tilt, s=1, c='k')
# plt.xlim(4,20)
# plt.ylim(100,300)
# plt.xlabel('r [kpc]')
# plt.ylabel('rotational v [km/s]')
# plt.savefig('/Users/amandaquirk/Desktop/RG_ss_vrot_HIclose_tilt.png')
# plt.close()

# plt.scatter(AGo_r, AGo_ss_vrot_tilt, s=1, c='m')
# plt.scatter(AGo_r, AGo_HIclose_vrot_tilt, s=1, c='k')
# plt.xlim(4,20)
# plt.ylim(100,300)
# plt.xlabel('r [kpc]')
# plt.ylabel('rotational v [km/s]')
# plt.savefig('/Users/amandaquirk/Desktop/AGo_ss_vrot_HIclose_tilt.png')
# plt.close()

# plt.scatter(AGy_r, AGy_ss_vrot_tilt, s=1, c='m')
# plt.scatter(AGy_r, AGy_HIclose_vrot_tilt, s=1, c='k')
# plt.xlim(4,20)
# plt.ylim(100,300)
# plt.xlabel('r [kpc]')
# plt.ylabel('rotational v [km/s]')
# plt.savefig('/Users/amandaquirk/Desktop/AGy_ss_vrot_HIclose_tilt.png')
# plt.close()

# #writing data to files
# file=open('/Users/amandaquirk/Desktop/MS_master_vrot.txt', 'w')
# file.write('#r (kpc), v_rot no model (km/s), v_rot tr model, smoothed PA, sPA v_rot no model, sPA v_rot tr model, n, HI main no model, HI main tr, HI close no model, HI close tr\n')
# for i in range(len(MS_xi)):
# 	file.write('{} {} {} {} {} {} {} {} {} {} {}\n'.format(MS_r[i],MS_si_vrot[i], MS_si_vrot_tilt[i],MS_sPA[i], MS_ss_vrot[i], MS_ss_vrot_tilt[i], MS_n[i], MS_HImain_vrot[i], MS_HImain_vrot_tilt[i], MS_HIclose_vrot[i], MS_HIclose_vrot_tilt[i]))
# file.close()

# file=open('/Users/amandaquirk/Desktop/AGy_master_vrot.txt', 'w')
# file.write('#r (kpc), v_rot no model (km/s), v_rot tr model, smoothed PA, sPA v_rot no model, sPA v_rot tr model, n, HI main no model, HI main tr, HI close no model, HI close tr\n')
# for i in range(len(AGy_xi)):
# 	file.write('{} {} {} {} {} {} {} {} {} {} {}\n'.format(AGy_r[i],AGy_si_vrot[i], AGy_si_vrot_tilt[i],AGy_sPA[i], AGy_ss_vrot[i], AGy_ss_vrot_tilt[i], AGy_n[i], AGy_HImain_vrot[i], AGy_HImain_vrot_tilt[i], AGy_HIclose_vrot[i], AGy_HIclose_vrot_tilt[i]))
# file.close()

# file=open('/Users/amandaquirk/Desktop/AGo_master_vrot.txt', 'w')
# file.write('#r (kpc), v_rot no model (km/s), v_rot tr model, smoothed PA, sPA v_rot no model, sPA v_rot tr model, n, HI main no model, HI main tr, HI close no model, HI close tr\n')
# for i in range(len(AGo_xi)):
# 	file.write('{} {} {} {} {} {} {} {} {} {} {}\n'.format(AGo_r[i],AGo_si_vrot[i], AGo_si_vrot_tilt[i],AGo_sPA[i], AGo_ss_vrot[i], AGo_ss_vrot_tilt[i], AGo_n[i], AGo_HImain_vrot[i], AGo_HImain_vrot_tilt[i], AGo_HIclose_vrot[i], AGo_HIclose_vrot_tilt[i]))
# file.close()

# file=open('/Users/amandaquirk/Desktop/RG_master_vrot.txt', 'w')
# file.write('#r (kpc), v_rot no model (km/s), v_rot tr model, smoothed PA, sPA v_rot no model, sPA v_rot tr model, n, HI main no model, HI main tr, HI close no model, HI close tr\n')
# for i in range(len(RG_xi)):
# 	file.write('{} {} {} {} {} {} {} {} {} {} {}\n'.format(RG_r[i],RG_si_vrot[i], RG_si_vrot_tilt[i],RG_sPA[i], RG_ss_vrot[i], RG_ss_vrot_tilt[i], RG_n[i], RG_HImain_vrot[i], RG_HImain_vrot_tilt[i], RG_HIclose_vrot[i], RG_HIclose_vrot_tilt[i]))
# file.close()


#infinite disk model
def Vrot_tilted_ring(v,PA_star): 
	i_ring = HI_i_avg
	PA_ring = HI_PA_avg
	vsys= -300 #km/s, as defined in Claire's thesis
	A=[float(a-vsys)/np.sin(float(i_ring*np.pi) / 180) for a in v]
	B= [(np.tan(float((PA_ring-b)*np.pi) / 180))**2 for b in PA_star]
	C= (np.cos(float(i_ring*np.pi) / 180))**2
	rotation_velocity= [a*np.sqrt(1+(float(b)/C)) for a,b in zip(A,B)]
	positive=[np.absolute(a) for a in rotation_velocity]
	return positive

MS_if_vrot_tilt=Vrot_tilted_ring(MS_v, MS_PA)
RG_if_vrot_tilt=Vrot_tilted_ring(RG_v, RG_PA)
AGo_if_vrot_tilt=Vrot_tilted_ring(AGo_v, AGo_PA)
AGy_if_vrot_tilt=Vrot_tilted_ring(AGy_v, AGy_PA)

MS_HImain_vrot_tr=Vrot_tilted_ring(MS_HImain,MS_PA)
AGy_HImain_vrot_tr=Vrot_tilted_ring(AGy_HImain,AGy_PA)
AGo_HImain_vrot_tr=Vrot_tilted_ring(AGo_HImain,AGo_PA)
RG_HImain_vrot_tr=Vrot_tilted_ring(RG_HImain,RG_PA)

from matplotlib import rc 
from matplotlib.ticker import MaxNLocator

rc('font', family = 'serif')
f, axes= plt.subplots(4,1, sharey=False, sharex=True, figsize=(4,9.8))
axes[0].scatter(MS_r, MS_HImain_vrot_tr, s=2, c='darkgray')
axes[0].scatter(MS_r, MS_if_vrot_tilt, s=2, c='b', alpha=0.4)
axes[1].scatter(AGy_r, AGy_HImain_vrot_tr, s=2, c='darkgray')
axes[1].scatter(AGy_r, AGy_if_vrot_tilt, s=2, c='m', alpha=0.4)
axes[2].scatter(AGo_r, AGo_HImain_vrot_tr, s=2, c='darkgray')
axes[2].scatter(AGo_r, AGo_if_vrot_tilt, s=2, c='k', alpha=0.4)
axes[3].scatter(RG_r, RG_HImain_vrot_tr, s=2, c='darkgray')
axes[3].scatter(RG_r, RG_if_vrot_tilt, s=2, c='r', alpha=0.4)
axes[0].annotate('MS', xy=(19,115), horizontalalignment='right', fontsize=12)
axes[1].annotate('young AGB', xy=(19,115), horizontalalignment='right', fontsize=12)
axes[2].annotate('older AGB', xy=(19,115), horizontalalignment='right', fontsize=12)
axes[3].annotate('RGB', xy=(19,115), horizontalalignment='right', fontsize=12)

for ax in axes:
	ax.set_xlim(4, 20)
	#ax.set_ylabel(r'$\rm v_{\rm rot} \ (km\ s^{-1})$', fontsize=13)
	ax.set_ylim(100,300)
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
axes[3].set_xlabel(r'$\rm Radial\ Distance:\ \it r \ \rm(kpc)$', fontsize=13)
nbins = len(axes[0].get_yticklabels())-1
axes[3].yaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='upper'))
axes[1].yaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='upper'))
axes[2].yaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='upper'))
f.subplots_adjust(left=0.17)
f.text(0.008, 0.5, r'$\rm Rotation\ Velocity:\ \it v_{\rm rot}\ \rm(km\ s^{-1})$', va='center', rotation='vertical', fontsize=13)
plt.subplots_adjust(wspace=0, hspace=0)
plt.savefig('/Users/amandaquirk/Desktop/rotation_curves_infinite_disk.pdf', bbox_inches='tight')
