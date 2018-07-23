from astropy.io import fits 
from astropy.table import Table 
from astropy import units as u
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt 
import numpy as np 
from shapely.geometry import Point 
from shapely.geometry.polygon import Polygon

#import the data file
hdulist= fits.open('/Users/amandaquirk/Documents/AsymmetricDrift/Data/subMasterSPLASH.fits', memmap=True)

#puts data into an accessible table
data=Table(hdulist[1].data)

#creates arrays of the parameters than I need
zqual=data['ZQUAL'] #quality of data point
L=data['LIKELIHOOD'] #if negative, MWFG star
RA=data['RA']
Dec=data['DEC']
F814W=data['F814W']
F475W=data['F475W']
z=data['Z'] #redshift

#these arrays will contain data from stars with good zqual values and that are not MWFGstars
z_good=[]
RA_good=[]
Dec_good=[]
F8W=[]
F4W=[]
for i in range(len(zqual)):
        #below selects the candidate stars
        if zqual[i]>=3 and L[i]>0 or np.isnan(L[i])==True:
                z_good.append(z[i])
                RA_good.append(RA[i])
                Dec_good.append(Dec[i])
                F8W.append(F814W[i])
                F4W.append(F475W[i])

#the next step is to separate the stars into age groups: MS, RGB, young AGB, and old AGB based on a color magnitude diagram. first i create a CDM

#below will contain the color of every good star
color=[]
for i in range(len(F8W)):
        color.append(F4W[i]-F8W[i])

#print(len(color))

#below will contain position and redshift values for MS stars
MS_z=[]
MS_RA=[]
MS_Dec=[]
MS_tag=[]
for i in range(len(color)):
        if color[i]<=1:
                MS_z.append(z_good[i])
                MS_RA.append(RA_good[i])
                MS_Dec.append(Dec_good[i])
                MS_tag.append(1)

#print(len(MS_z), len(MS_RA), len(MS_Dec))

#below tests to see if the star is in the RGB
RG_z=[]
RG_RA=[]
RG_Dec=[]
RG_tag=[]
#polygon creates a shape that corresponds to the area of the CMD that contains RGB stars
polygon= Polygon([(2,23), (2.7, 20.4), (7.5,23)])
for i in range(len(color)):
        point= Point(color[i], F8W[i])
        if polygon.contains(point)==True:
                RG_z.append(z_good[i])
                RG_RA.append(RA_good[i])
                RG_Dec.append(Dec_good[i])
                RG_tag.append(4)

#print(len(RG_z), len(RG_RA), len(RG_Dec))

#below tests to see if the star is in the old AGB
AGo_z=[]
AGo_RA=[]
AGo_Dec=[]
AGo_tag=[]
#polygon creates a shape that corresponds to the area of the CMD that contains old AGB stars
polygon= Polygon([(2.7,20.4), (7.5, 23), (8, 20.5)])
for i in range(len(color)):
        point= Point(color[i], F8W[i])
        if polygon.contains(point)==True:
                AGo_z.append(z_good[i])
                AGo_RA.append(RA_good[i])
                AGo_Dec.append(Dec_good[i])
                AGo_tag.append(3)

#print(len(AGo_z), len(AGo_RA), len(AGo_Dec))

#below tests to see if the star is in the young AGB
AGy_z=[]
AGy_RA=[]
AGy_Dec=[]
AGy_tag=[]
#polygon creates a shape that corresponds to the area of the CMD that contains young AGB stars
polygon= Polygon([(3.5,18), (8,18), (2.7, 20.4), (3.5,20.5), (8,20.5)])
for i in range(len(color)):
        point= Point(color[i], F8W[i])
        if polygon.contains(point)==True:
                AGy_z.append(z_good[i])
                AGy_RA.append(RA_good[i])
                AGy_Dec.append(Dec_good[i])
                AGy_tag.append(2)

#print(len(AGy_z), len(AGy_RA), len(AGy_Dec))

#eliminates the outliers in the young AG stars
AGy_z_good=[]
AGy_RA_good=[]
AGy_Dec_good=[]
AGy_tag_good=[]
for i in range(len(AGy_RA)):
        if AGy_z[i]<(200/300000):
                AGy_z_good.append(AGy_z[i])
                AGy_RA_good.append(AGy_RA[i])
                AGy_Dec_good.append(AGy_Dec[i])
                AGy_tag_good.append(2)   

#combine arrays
RA_final=MS_RA+RG_RA+AGy_RA_good+AGo_RA
Dec_final=MS_Dec+RG_Dec+AGy_Dec_good+AGo_Dec
z_final=MS_z+RG_z+AGy_z_good+AGo_z
tag=MS_tag+RG_tag+AGy_tag_good+AGo_tag

#convert redshift to velocity
velocity=[]
for i in range(len(z_final)):
        velocity.append(z_final[i]*3*(10**5))

#print(len(RA_final), len(Dec_final), len(velocity), len(tag))

#create a fits file
c1=fits.Column(name='RA', array=RA_final, format='11A')
c2=fits.Column(name='Dec', array=Dec_final, format='11A')
c3=fits.Column(name='Velocity', array=velocity, format='D')
c4=fits.Column(name='Age', array=tag, format='K')
t= fits.BinTableHDU.from_columns([c1, c2, c3,c4])
t.writeto('/Users/amandaquirk/Desktop/submasterSPLASH_coordinates_agebins.fits')




