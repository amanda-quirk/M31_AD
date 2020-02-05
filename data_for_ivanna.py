import numpy as np 
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord, EarthLocation
from astropy import constants as const
from astropy.time import Time
import matplotlib.pyplot as plt
from astropy.table import Table

'''
goal: give Ivanna the data of star's in Claire's 2012 and 2013 MCT + PHAT masks
'''

#original data ================================================================================================================================
hdu = fits.open('/Users/amandaquirk/Documents/AsymmetricDrift/Data/subMasterSPLASH.fits', memmap = True)
data = hdu[1].data
ra = data['RA']
dec = data['Dec']
redshift = data['Z']
time = data['MJD']
aband = data['ABAND']
SN = data['SSNR']
mask = data['MASK']
name = data['OBJNAME']
grating = data['GRATE']
flux = data['SPEC']
wavelength = data['LBIN']
F275W = data['F275W']
F336W = data['F336W']
F475W = data['F475W']
F814W = data['F814W']
F110W = data['F110W']
F160W = data['F160W']
print('read in subMasterSPLASH')

'''
old version before I read the directions carefully :)
#convert to eta and xi and do a location cut==================================================================================================
#shift to frame centered on M31
# m31 = SkyCoord(10.6847083*u.deg, 41.26875*u.deg, frame='icrs')
# c = SkyCoord(ra=ra, dec=dec, frame='icrs', unit=(u.hourangle, u.deg))
# c_inm31 = c.transform_to(m31.skyoffset_frame())
# xi, eta = c_inm31.lon, c_inm31.lat
# xi = xi.deg
# eta = eta.deg

# #cut to approximately the MCT + PHAT masks
# window = (0.2 < eta) & (eta < 1.2) & (0.1 < xi) & (xi < 1)
# ra = ra[window]
# dec = dec[window]
# redshift = redshift[window]
# time = time[window]
# aband = aband[window]
# SN = SN[window]
# mask = mask[window]
# name = name[window]
# grating = grating[window]
# flux = flux[window]
# wavelength = wavelength[window]
# F275W = F275W[window]
# F336W = F336W[window]
# F475W = F475W[window]
# F814W = F814W[window]
# F110W = F110W[window]
# F160W = F160W[window]

# #only want targets with S/N > 15==============================================================================================================
# high_SN = SN > 15
# ra = ra[high_SN]
# dec = dec[high_SN]
# redshift = redshift[high_SN]
# time = time[high_SN]
# aband = aband[high_SN]
# SN = SN[high_SN]
# mask = mask[high_SN]
# name = name[high_SN]
# grating = grating[high_SN]
# flux = flux[high_SN]
# wavelength = wavelength[high_SN]
# F275W = F275W[high_SN]
# F336W = F336W[high_SN]
# F475W = F475W[high_SN]
# F814W = F814W[high_SN]
# F110W = F110W[high_SN]
# F160W = F160W[high_SN]

# #only want stars with the 1200G line grating==================================================================================================
# G1200 = grating == 1200
# ra = ra[G1200]
# dec = dec[G1200]
# redshift = redshift[G1200]
# time = time[G1200]
# aband = aband[G1200]
# SN = SN[G1200]
# mask = mask[G1200]
# name = name[G1200]
# grating = grating[G1200]
# flux = flux[G1200]
# wavelength = wavelength[G1200]
# F275W = F275W[G1200]
# F336W = F336W[G1200]
# F475W = F475W[G1200]
# F814W = F814W[G1200]
# F110W = F110W[G1200]
# F160W = F160W[G1200]

# print('The mask names are', mask)
# #plotting check
# c = SkyCoord(ra=ra, dec=dec, frame='icrs', unit=(u.hourangle, u.deg))
# c_inm31 = c.transform_to(m31.skyoffset_frame())
# xi, eta = c_inm31.lon, c_inm31.lat
# xi = xi.deg
# eta = eta.deg
# plt.scatter(eta, xi)
# plt.xlim(1,0)
# plt.show()
'''

'''
new version where I actually read the paper and am going by mask name instead wow oh wow
'''

#MCT = ((mask=='mctA5') | (mask=='mctB4') | (mask=='mctC3') | (mask=='mctD3') | (mask=='mctE3')) & (np.isnan(F475W) == False) & (name.startswith('serendip') == False) #looking at mask names in the paper but want to remove serendips and also stars that don't have PHAT photometry; this is the version just from Claire's 2012 paper

#looking at all the masks that have PHAT photometry in the area Ivanna is looking into
masks_1200 = ['mctG', 'mctA5', 'mct15p', 'mct06p', 'mctJ', 'mctE3', 'mct07p', 'mct05p', 'mct12p', 'mct16p', 'mctL', 'mct09p', 'mct04p', 'mctC3', 'mct10p', 'mctK', 'mctF', 'mctD3', 'mctB4', 'mct13p']
masks_600 = ['mct6U', 'mct6K', 'mct6H', 'mct6R', 'mct6G', 'mct6O', 'mct6V', 'mct6X', 'mct6P', 'mct6T', 'mct6L', 'mct6D', 'mct6M', 'mct6F', 'mct6C', 'mct6Q', 'mct6W', 'mct6S']

MCT_1200 = (np.isin(mask, masks_1200)) & (np.isnan(F475W) == False) & (name.startswith('serendip') == False) 
MCT_600 = (np.isin(mask, masks_600)) & (np.isnan(F475W) == False) & (name.startswith('serendip') == False)


ra = ra[MCT_1200]
dec = dec[MCT_1200]
redshift = redshift[MCT_1200]
time = time[MCT_1200]
aband = aband[MCT_1200]
SN = SN[MCT_1200]
mask = mask[MCT_1200]
name = name[MCT_1200]
grating = grating[MCT_1200]
flux = flux[MCT_1200]
wavelength = wavelength[MCT_1200]
F275W = F275W[MCT_1200]
F336W = F336W[MCT_1200]
F475W = F475W[MCT_1200]
F814W = F814W[MCT_1200]
F110W = F110W[MCT_1200]
F160W = F160W[MCT_1200]
print(sum(MCT_1200))

#correct the velocities ======================================================================================================================
def correct_vel(ra, dec, time, redshift, aband):
	sc = SkyCoord(ra=ra, dec=dec, unit=(u.hourangle, u.deg))

	keck = EarthLocation.from_geodetic(lat=19.8283*u.deg, lon=-155.4783*u.deg, height=4160*u.m)

	heliocorr = sc.radial_velocity_correction('heliocentric', obstime=Time(time, format='mjd'), location=keck) 
	heliocorr_km_s = heliocorr.to(u.km/u.s) 
	vraw = redshift * const.c.to(u.km/u.s)
	vcorr = vraw + heliocorr_km_s - aband * const.c.to(u.km/u.s)

	return vcorr.value #km/s

velocity = correct_vel(ra, dec, time, redshift, aband)

#save the data================================================================================================================================
np.savetxt('/Users/amandaquirk/Desktop/MCT_PHAT_stars_velocity_1200_all.txt', np.c_[name, velocity], header='object name, helio + aband corrected vel (km/s)', fmt='%s')

#saving the spectra probably ===============================================================================================================
MCT_data = data[MCT_1200]

t = Table(MCT_data)
t.write('/Users/amandaquirk/Desktop/MCT_PHAT_stars_catalogue_spectra_1200_all.fits', format='fits')




