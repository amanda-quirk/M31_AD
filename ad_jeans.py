import numpy as np 
from scipy.optimize import fsolve
import matplotlib.pyplot as plt 

'''
purpose: use Jean's equation and estimated velocity ellipsoid shape to predict the AD
'''

#import data here
HI_r, HI_PA, HI_i, HI_v=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/HI_PA_i_vrot.txt', unpack=True)
RG_r = np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/RG_master_vrot.txt', usecols=(0,), unpack=True)
RG_var = np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/RG_smoothed_chemin.txt', usecols=(4), unpack=True)

#below defines a function to determine which HI ring a star is cloest to
def find_nearest_ring(radius):
	idx=[]
	for i in range(len(HI_r)):
		idx.append(np.abs(HI_r[i]-radius))
		n=np.argmin(idx)
	return n #index of the radius closest to the star's radius

def jeans_va(r, sigmaLOS): #kpc and km/s
	# k = 1
	# PA = 38
	# vc = 260
	# inc = 77
	k = np.zeros_like(r)
	vc = np.zeros_like(r)
	PA = np.zeros_like(r)
	inc = np.zeros_like(r)
	for i in range(len(k)):
		if r[i] < 10:
			k[i] = 2
		else: 
			k[i] = 1 

		ind = find_nearest_ring(r[i])
		print(ind)
		vc[i] = HI_v[ind]
		PA[i] = HI_PA[ind]
		inc[i] = HI_i[ind]

	def sine(angle): #degrees -- convert to radians
		return np.sin(angle * np.pi / 180)

	def cosine(angle): #degrees -- convert to radians
		return np.cos(angle * np.pi / 180)

	sigmaR2 = sigmaLOS**2 / ((((sine(PA)**2 + (0.64 * cosine(PA)**2)) * sine(inc)**2)) + (0.2 * cosine(inc)**2))

	#solve quad formula
	a = -1
	b = 2 * vc  
	c = (sigmaR2 / (2 * vc)) * ((k * r / 5.76) - 0.76)

	va_solution_p = (-b + np.sqrt(b**2 - (4 * a * c))) / (2 * a) 
	va_solution_n = (-b - np.sqrt(b**2 - (4 * a * c))) / (2 * a) 

	sigmaPHI2 = (0.8 * np.sqrt(sigmaR2))**2

	#solving using fsolve
	# func = lambda va: -va**2 + (2 * va * vc) + sigmaR2 / (2 * vc) * ((k * r) / 5.76 - 0.76)
	# initial_guess = sigmaPHI2
	# va_solution = fsolve(func, initial_guess)

	return va_solution_p, va_solution_n, sigmaPHI2

yp, yn, sigmaPHI2 = jeans_va(RG_r, RG_var)
# plt.scatter(sigmaPHI2, yp, c='r')
plt.scatter(sigmaPHI2, np.sqrt(yn), c='b')
# # plt.xlim(0, 20000)
# # plt.ylim(-20, 140)
plt.show()
# RG_ind = []
# for i in range(len(RG_r)):
# 	RG_ind.append(find_nearest_ring(RG_r[i]))
# plt.scatter(RG_r, RG_ind)
# plt.show()

