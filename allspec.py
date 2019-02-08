#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  8 11:29:04 2019

@author: Christina
"""

#Libraries for plotting
import matplotlib
import matplotlib.pyplot as plt

#Standard numerical and scientific libraries
import numpy as np
import scipy as sp

#Units, Constants, Fits reading libraries
from astropy import units as U
from astropy import constants as C
from astropy.io import fits
from astropy.io import ascii

#Automatic spectra reading tools
from specutils.io import read_fits

#Astronomy related packages
from PyAstronomy import pyasl
from specutils import extinction

import sfdmap  #for the SFD dust maps 

from decimal import Decimal #convert to decimals

#import pyfits
from string import rstrip

from matplotlib import rc   #for Latex text

#Load files:
mags = np.loadtxt("/home/christina/Desktop/QSOs/mags.list")
comp_file=np.loadtxt('/home/christina/Desktop/QSOs/compoM.data')
zab = np.loadtxt("/home/christina/Desktop/QSOs/zab1.list")
names = np.loadtxt("/home/christina/Desktop/QSOs/names.list", dtype=np.str)

#Read the coordinates:
f = open("/home/christina/Desktop/QSOs/coord.list", "r")
e=open("/home/christina/Desktop/QSOs/ebv.list", "w")    #open file to write 
rah,decd = [], []
for l in f:
    row = l.split()
    rah = row[0]
    decd = row[1]
    coor= np.array([rah,decd])
    a=' '.join(coor)
    ra, dec = pyasl.coordsSexaToDeg(a)
    m = sfdmap.SFDMap('/home/christina/Desktop/QSOs/Dust/')
    ebv = m.ebv(ra, dec, frame='icrs', interpolate=True)
    b=(''.join(str(ebv)))
    e.write("%s\n" % b)       #write ebv values to list file
    for row in str(ebv):
        b=(''.join(str(ebv)))
e.close()

ebv = np.loadtxt('/home/christina/Desktop/QSOs/ebv.list')
AB=zab[:,1]
z=zab[:,0]

wlr=comp_file[:,0]/1e4   #wavelength of composite spectrum in microns  (rest frame)
flx=comp_file[:,1]   #flux of composite spectrum

#Parameters for the Pei-parametrization of the SMC extinction curve
ai = np.array([185.,27.,0.005,0.010,0.012,0.030])
wli = np.array([0.042,0.08,0.22,9.7,18.,25.])
bi = np.array([90.,5.50,-1.95,-1.95,-1.80,0.0])
ni = np.array([2.0,4.0,2.0,2.0,2.0,2.0])
Ki = np.array([2.89,0.91,0.02,1.55,1.72,1.89])

### Central wavelengths for broad band filters (Aangstrom)
lambdau=3540.
lambdag=4750.
lambdar=6220.
lambdai=7630.
lambdaz=9050.
lambday=10310.
lambdaj=12480.
lambdah=16310.
lambdak=22010.
#Build an array with the wavelengths
wavar=np.array([lambdau, lambdag, lambdar, lambdai, lambdaz, lambday, lambdaj,lambdah, lambdak])  

#The magnitudes
u = mags[:,0]
g = mags[:,1]
r = mags[:,2]
i = mags[:,3]
zs = mags[:,4] #zs in order to not confuse with the redshift
y = mags[:,5]
j = mags[:,6]
h = mags[:,7]
k = mags[:,8]
#Build an array with the magnitudes
magsar=np.array([u,g,r,i,zs,y,j,h,k])
            
#Now run through the Nspec spectra
Nspec=18

list = open('/home/christina/Desktop/QSOs/redspec1.list')

wavelength = np.zeros((Nspec),float)
flux = np.zeros((Nspec),float)

for s in range(0,Nspec):
    spec = fits.open(rstrip(list.readline()))
    Header = spec[0].header
    flux = spec[0].data*1e17
    lambda_central = Header['CRVAL1']
    delta_lambda = Header['CD1_1']
    Naxis1 = spec[0].header['NAXIS1']
    wavelength = delta_lambda*np.arange(0,Naxis1)+lambda_central
    fluxe=flux[3]   #error specrum flux
    #print fluxe
    flux=flux[0]
    
    #Convert mags to flux for each filter:
    fluxu = 10.**(0.4*(-48.60-u[s]))*3.e18/lambdau**2    
    fluxg = 10.**(0.4*(-48.60-g[s]))*3.e18/lambdag**2
    fluxr = 10.**(0.4*(-48.60-r[s]))*3.e18/lambdar**2
    fluxi = 10.**(0.4*(-48.60-i[s]))*3.e18/lambdai**2
    fluxz = 10.**(0.4*(-48.60-zs[s]))*3.e18/lambdaz**2
    fluxy = 10.**(0.4*(-48.60-y[s]))*3.e18/lambday**2
    fluxj = 10.**(0.4*(-48.60-j[s]))*3.e18/lambdaj**2
    fluxh = 10.**(0.4*(-48.60-h[s]))*3.e18/lambdah**2
    fluxk = 10.**(0.4*(-48.60-k[s]))*3.e18/lambdak**2
    #convert to array
    fluxar=np.array([fluxu, fluxg, fluxr, fluxi, fluxz, fluxy, fluxj,fluxh, fluxk])  
    AB=zab[:,1]
    ABn= AB[s]
    
    #Correct for galactic extinction
    funred = pyasl.unred(wavelength, flux[0], ebv[s])
    deredr = pyasl.unred(wavar, fluxar, ebv[s])
    spec = spec[0]
    spec = funred
    filt = np.nonzero((wavelength*U.Angstrom > (5700*U.Angstrom)) & (wavelength*U.Angstrom < (6600*U.Angstrom)))
    norm = np.mean(spec[filt])
    specr = funred/norm*deredr[2]*1.e17    #Scaling spectrum tor the r-band--the index 2 refers to the r-band

    z=zab[:,0]
    zn=z[s]
    
    wl=comp_file[:,0]*(zn+1.)   #wavelength of composite spectrum
    flx=comp_file[:,1]    #flux of composite spectrum

    filtc = np.nonzero((wl*U.Angstrom > (5700*U.Angstrom)) & (wl*U.Angstrom < (6600*U.Angstrom)))
    normc=np.mean(flx[filtc])       #normalized flux of comp spec without reddening
    
    #Correct the reddened composite spectrum for SMC reddening: 
    Alambda = flx*0.      #Reset flux to zero
    for e in range(len(ai)):
        Alambda=Alambda+ai[e]/((wlr/wli[e])**ni[e]+(wli[e]/wlr)**ni[e]+bi[e])   #fitting function covering all wavelengths in microns
    Alambda = Alambda*ABn
    model = 10**(-0.4*Alambda)*flx    #use the flux of the composite spectrum
   # wl=comp_file[:,0]*(z+1.)   #wavelength of the composite spectrum
    flt = np.nonzero((wl*U.Angstrom > (5700*U.Angstrom)) & (wl*U.Angstrom < (6600*U.Angstrom)))
    normm = np.mean(model[flt])

    #Error Spectrum
    
    funrede = pyasl.unred(wavelength, fluxe[0], ebv[s])
    #spece = spec[0]
    spece = funrede
    spece = fluxe/norm*deredr[2]*1e17   #Scale error spectrum to r-band
    spece=spece[0]
    
    ########## PLOT ############
    
    fig = plt.figure(figsize=(20,10))
    ax = fig.add_subplot(111)
    ax1=fig.add_subplot(111)
    ax2=fig.add_subplot(111) 
    ax3=fig.add_subplot(111) 
    ax4=fig.add_subplot(111) 
    ax5=fig.add_subplot(111) 

    ax.set_xlabel("Wavelength (Angstrom)", fontsize=20)
    ax.set_ylabel("Flux (10$^-17$ erg/cm$^2$/s/Angstrom)",fontsize=20)
    ax.set_xticklabels(['5000','10000','20000'])
    #plt.xticks(range(3), ('5000', '10000', '20000'))

    #ax.semilogx(wavelength,specr,'k')
    #ax1.semilogx(wavelength,spece,'b')    #error spectrum

    
    if names[s] == str('BrightJ2103-0043'):
        plt.text(15000,100,names[s], fontsize=20)
        plt.text(15000,94, 'z=', fontsize=20)
        plt.text(15800,94,zn, fontsize=20)
        plt.text(17300,94,r'$A_B$=', fontsize=20)
        plt.text(18600,94,ABn, fontsize=20)
        ax.axis([3800,22000,0,110])
        ax.semilogx(wavelength,specr,'k')
        ax1.semilogx(wavelength,spece,color='grey')    #error spectrum
        ax3.semilogx(wl,flx/normc*deredr[2]*1.e17,'g')   #composite spectrum without SMC reddening
        ax4.semilogx(wavar, deredr*1e17, 'ro', markersize=14)   #photometric points 
        ax5.semilogx(wl, model/normm*deredr[2]*1.e17,'r')    #composite spectrum with SMC reddening
    elif names[s] == str('RedJ2113-0028'):
        plt.text(15000,9,names[s], fontsize=20)
        plt.text(15000,7, 'z=', fontsize=20)
        plt.text(15800,7,zn, fontsize=20)
        plt.text(17300,7,r'$A_B$=', fontsize=20)
        plt.text(18600,7,ABn, fontsize=20)
        ax.axis([3800,22000,0,10])
        ax.semilogx(wavelength,specr,'k')
        ax1.semilogx(wavelength,spece,color='grey')    #error spectrum
        ax3.semilogx(wl,flx/normc*deredr[2]*1.e17,'g')   #composite spectrum without SMC reddening
        ax4.semilogx(wavar, deredr*1e17, 'ro', markersize=14)   #photometric points 
        ax5.semilogx(wl, model/normm*deredr[2]*1.e17,'r')    #composite spectrum with SMC reddening
    elif names[s] == str('RedJ2045-0022'):
        plt.text(15000,30,names[s], fontsize=20)
        plt.text(15000,28, 'z=', fontsize=20)
        plt.text(15800,28,zn, fontsize=20)
        plt.text(17300,28,r'$A_B$=', fontsize=20)
        plt.text(18600,28,ABn, fontsize=20)
        ax.axis([3800,22000,0,32])
        ax.semilogx(wavelength,specr,'k')
        ax1.semilogx(wavelength,spece,color='grey')    #error spectrum
        ax3.semilogx(wl,flx/normc*deredr[2]*1.e17,'g')   #composite spectrum without SMC reddening
        ax4.semilogx(wavar, deredr*1e17, 'ro', markersize=14)   #photometric points 
        ax5.semilogx(wl, model/normm*deredr[2]*1.e17,'r')    #composite spectrum with SMC reddening
    elif names[s] == str('RedJ0312+0032'):
        plt.text(15000,13,names[s], fontsize=20)
        plt.text(15000,11, 'z=', fontsize=20)
        plt.text(15800,11,zn, fontsize=20)
        plt.text(17300,11,r'$A_B$=', fontsize=20)
        plt.text(18600,11,ABn, fontsize=20)
        ax.axis([3800,22000,0,15])
        ax.semilogx(wavelength,specr,'k')
        ax1.semilogx(wavelength,spece,color='grey')    #error spectrum
        ax3.semilogx(wl,flx/normc*deredr[2]*1.e17,'g')   #composite spectrum without SMC reddening
        ax4.semilogx(wavar, deredr*1e17, 'ro', markersize=14)   #photometric points 
        ax5.semilogx(wl, model/normm*deredr[2]*1.e17,'r')    #composite spectrum with SMC reddening
    elif names[s] == str('RedJ0312+0035'):
        plt.text(15000,18,names[s], fontsize=20)
        plt.text(15000,16, 'z=', fontsize=20)
        plt.text(15800,16,zn, fontsize=20)
        plt.text(17300,16,r'$A_B$=', fontsize=20)
        plt.text(18600,16,ABn, fontsize=20)
        ax.axis([3800,22000,0,20])   
        ax.semilogx(wavelength,specr,'k')
        ax1.semilogx(wavelength,spece,color='grey')    #error spectrum
        ax3.semilogx(wl,flx/normc*deredr[2]*1.e17,'g')   #composite spectrum without SMC reddening
        ax4.semilogx(wavar, deredr*1e17, 'ro', markersize=14)   #photometric points 
        ax5.semilogx(wl, model/normm*deredr[2]*1.e17,'r')    #composite spectrum with SMC reddening
    elif names[s] == str('BrightJ2241-0103'):
        plt.text(15000,56,names[s], fontsize=20)
        plt.text(15000,53, 'z=', fontsize=20)
        plt.text(15800,53,zn, fontsize=20)
        plt.text(17300,53,r'$A_B$=', fontsize=20)
        plt.text(18600,53,ABn, fontsize=20)
        ax.axis([3800,22000,0,60])
        ax.semilogx(wavelength,specr,'k')
        ax1.semilogx(wavelength,spece,color='grey')    #error spectrum
        ax3.semilogx(wl,flx/normc*deredr[2]*1.e17,'g')   #composite spectrum without SMC reddening
        ax4.semilogx(wavar, deredr*1e17, 'ro', markersize=14)   #photometric points 
        ax5.semilogx(wl, model/normm*deredr[2]*1.e17,'r')    #composite spectrum with SMC reddening
    elif names[s] == str('BrightJ2259+0103'):
        plt.text(15000,56,names[s], fontsize=20)
        plt.text(15000,53, 'z=', fontsize=20)
        plt.text(15800,53,zn, fontsize=20)
        plt.text(17300,53,r'$A_B$=', fontsize=20)
        plt.text(18600,53,ABn, fontsize=20)
        ax.axis([3800,22000,0,60])
        ax.semilogx(wavelength,specr,'k')
        ax1.semilogx(wavelength,spece,color='grey')    #error spectrum
        ax3.semilogx(wl,flx/normc*deredr[2]*1.e17,'g')   #composite spectrum without SMC reddening
        ax4.semilogx(wavar, deredr*1e17, 'ro', markersize=14)   #photometric points 
        ax5.semilogx(wl, model/normm*deredr[2]*1.e17,'r')    #composite spectrum with SMC reddening
    elif names[s] == str('RedJ2320+0028'):
        plt.text(15000,6,names[s], fontsize=20)
        plt.text(15000,5, 'z=', fontsize=20)
        plt.text(15800,5,zn, fontsize=20)
        plt.text(17300,5,r'$A_B$=', fontsize=20)
        plt.text(18600,5,ABn, fontsize=20)
        ax.axis([3800,22000,0,7])
        ax.semilogx(wavelength,specr,'k')
        ax1.semilogx(wavelength,spece,color='grey')    #error spectrum
        ax4.semilogx(wavar, deredr*1e17, 'ro', markersize=14)   #photometric points 
    elif names[s] == str('RedJ2047+0045'):
        plt.text(15000,48,names[s], fontsize=20)
        plt.text(15000,45, 'z=', fontsize=20)
        plt.text(15800,45,zn, fontsize=20)
        plt.text(17300,45,r'$A_B$=', fontsize=20)
        plt.text(18600,45,ABn, fontsize=20)
        ax.axis([3800,22000,0,50])
        ax.semilogx(wavelength,specr,'k')
        ax1.semilogx(wavelength,spece,color='grey')    #error spectrum
        ax3.semilogx(wl,flx/normc*deredr[2]*1.e17,'g')   #composite spectrum without SMC reddening
        ax4.semilogx(wavar, deredr*1e17, 'ro', markersize=14)   #photometric points 
        ax5.semilogx(wl, model/normm*deredr[2]*1.e17,'r')    #composite spectrum with SMC reddening
    elif names[s] == str('BrightJ2043+0045'):
        plt.text(15000,48,names[s], fontsize=20)
        plt.text(15000,45, 'z=', fontsize=20)
        plt.text(15800,45,zn, fontsize=20)
        plt.text(17300,45,r'$A_B$=', fontsize=20)
        plt.text(18600,45,ABn, fontsize=20)
        ax.axis([3800,22000,0,52])
        ax.semilogx(wavelength,specr,'k')
        ax1.semilogx(wavelength,spece,color='grey')    #error spectrum
        ax3.semilogx(wl,flx/normc*deredr[2]*1.e17,'g')   #composite spectrum without SMC reddening
        ax4.semilogx(wavar, deredr*1e17, 'ro', markersize=14)   #photometric points 
        ax5.semilogx(wl, model/normm*deredr[2]*1.e17,'r')    #composite spectrum with SMC reddening
    #elif names[s] == str('RedJ2048+0056'):
        #plt.text(15000,30,names[s], fontsize=20)
        #plt.text(5223,12, 'CIV', fontsize='18')
        #plt.text(4200,15, r'Ly-$\alpha$', fontsize='18')
        #plt.text(6436,15, r'CIII', fontsize='18')
        #plt.text(9422,15, r'MgII', fontsize='18')
        #plt.text(4720,15, r'SiIV', fontsize='18')
        #plt.text(15000,28, 'z=', fontsize=20)
        #plt.text(15800,28,zn, fontsize=20)
        #plt.text(17300,28,r'$A_B$=', fontsize=20)
        #plt.text(18600,28,ABn, fontsize=20)
        #ax.axis([3800,22000,0,32])
        #ax.semilogx(wavelength,specr,'k')
        #ax1.semilogx(wavelength,spece,color='grey')    #error spectrum
        #ax3.semilogx(wl,flx/normc*deredr[2]*1.e17,'g')   #composite spectrum without SMC reddening
        #ax4.semilogx(wavar, deredr*1e17, 'ro', markersize=14)   #photometric points 
        #ax5.semilogx(wl, model/normm*deredr[2]*1.e17,'r')    #composite spectrum with SMC reddening
        #l1=[5214,5214]
        #l2=[6.5,32]
        #l3=[4092,4092]
        #l4=[6.5,32]
        #l5=[6429,6429]
        #l6=[6.5,32]
        #l7=[9417,9417]
        #l8=[6.5,32]
        #l9=[4715,4715]
        #l10=[6.5,32]
        #plt.plot(l1,l2,'--', color = 'grey')
        #plt.plot(l3,l4,'--', color = 'grey')
        #plt.plot(l5,l6,'--', color = 'grey')
        #plt.plot(l7,l8,'--', color = 'grey')
        #plt.plot(l9,l10,'--', color = 'grey')
        #ax3.semilogx(wl,flx/normc*deredr[2]*1.e17,'g')   #composite spectrum without SMC reddening
    else:
        plt.text(15000,30,names[s], fontsize=20)
        plt.text(15000,28, 'z=', fontsize=20)
        plt.text(15800,28, zn, fontsize=20)
        plt.text(17300,28,r'$A_B$=', fontsize=20)
        plt.text(18600,28,ABn, fontsize=20)
        ax.axis([3800,22000,0,32])
        ax.semilogx(wavelength,specr,'k')
        ax1.semilogx(wavelength,spece,color='grey')    #error spectrum
        ax3.semilogx(wl,flx/normc*deredr[2]*1.e17,'g')   #composite spectrum without SMC reddening
        ax4.semilogx(wavar, deredr*1e17, 'ro', markersize=14)   #photometric points 
        ax5.semilogx(wl, model/normm*deredr[2]*1.e17,'r')    #composite spectrum with SMC reddening
    
    print(names[s],z[s],AB[s],ebv[s],deredr[4]*1.e17)
 
   #Composite spectrum
    ax2=ax.twiny()
    ax2.set_xlabel("Rest Wavelength (Angstrom)", fontsize=20)
    #plt.setp(ax,xticks=[5000, 10000, 20000], xticklabels=['5000', '10000', '20000']) 
    #plt.xticks(np.arange(3), ('1000', '2000', '5000'))
    #ax3.semilogx(wl,flx/normc*deredr[2]*1.e17,'g')   #composite spectrum without SMC reddening
    #ax4.semilogx(wavar, deredr*1e17, 'ro', markersize=14)   #photometric points 
    #ax5.semilogx(wl, model/normm*deredr[2]*1.e17,'r')    #composite spectrum with SMC reddening
  
    plt.show()
list.close()
