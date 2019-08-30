#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 12:00:47 2019

@author: seth
"""

import numpy as np
from astropy.table import Table, join, vstack, Column
import matplotlib.pyplot as plt
import pdb

#Solar luminosites from http://mips.as.arizona.edu/~cnaw/sun.html
msolv=4.81
msolr=4.65 #SDSS AB
msolb=5.44

s17=Table.read('spengler17_tab8.dat',format='ascii')
p11_1=Table.read('paudel11_tab1.tex')
p11_2=Table.read('paudel11_tab2.tex')
p11=join(p11_1,p11_2,keys='galaxy')
k09=Table.read('koleva09_tab5.tex')
gallazzi05=Table.read('gallazzi05_tab2.dat',format='ascii')
spname=['M54','NGC5102','NGC5206']
splogmass=[np.log10(3.0e8),np.log10(6.0e9),np.log10(2.4e9)]
spfeh=[-1.49,-0.38,-0.16]
# M54: mass from Niederste-Ostholt 2010/12, metallicity from Harris, 
#NGC5102 and 5206 stellar masses from Nguyen+ 2018, metallicities from Kacharov+ 2018
bc03=Table.read('bc03.fits')
bc03['feh']=np.log10(bc03['z']/0.019)

#now get galaxy masses for Paudel
p11['mlv']=0.0 #for storing the M/LV
p11['V-r']=0.0
for i in range(0,len(p11)):
#interpolate the M/L for Paudel galaxies using BC03 SSP models of closest metallicity
    argmin=np.argmin(np.abs(p11['fehgal'][i]-bc03['feh']))
    use=(bc03['feh'] == bc03['feh'][argmin])
#    print(p11['agegal'][i],np.interp(np.log10(p11['agegal'][i]*1.e9), bc03['log-age-yr'][use], bc03['M*/Lv'][use]))
    p11['mlv'][i]=np.interp(np.log10(p11['agegal'][i]*1.e9), bc03['log-age-yr'][use], bc03['M*/Lv'][use])
    p11['V-r'][i]=np.interp(np.log10(p11['agegal'][i]*1.e9), bc03['log-age-yr'][use], bc03['V-R'][use])-0.30
    #note that for old populations > 1 Gyr, Padova models show JohnsonR-SDSSr=-0.30 with a pretty ~0.03 mag variations
p11['mlr']=p11['mlv']+10**(-0.4*(p11['V-r']-(msolv-msolr)))
p11['logmassgal']=np.log10(10**(-0.4*(p11['mrgal']-msolr))*p11['mlr'])

#get galaxy masses for Koleva
k09['mlb']=0.0
for i in range(0,len(k09)):
    argmin=np.argmin(np.abs(k09['fehgal'][i]-bc03['feh']))
    use=(bc03['feh'] == bc03['feh'][argmin])
    k09['mlb'][i]=np.interp(np.log10(k09['agegal'][i]*1.e9), bc03['log-age-yr'][use], bc03['M*/Lv'][use])
k09['logmassgal']=np.log10(10**(-0.4*(k09['mb']-k09['distmod']-msolb))*k09['mlb'])

#Mass-metallicity relation from Kirby+ 2013: https://ui.adsabs.harvard.edu/abs/2013ApJ...779..102K/abstract
logmassarr=np.arange(8.0,9.2,0.1)
kirby13feh=(logmassarr-6.)*0.30-1.69

ufontsize=24
plt.rc('text', usetex=True)
plt.rc('font', family='sans-serif')
plt.rc('xtick', labelsize=ufontsize) 
plt.rc('ytick', labelsize=ufontsize) 
plt.figure(figsize=(6,6))
#   plt.xlim(7.5,10.4)
plt.ylim(-1.7,0.7)
plt.xlabel('log($M_{\star}$)',fontsize=ufontsize)
plt.ylabel('[Fe/H]',fontsize=ufontsize)

plt.plot(logmassarr,kirby13feh,linestyle='-',color='black')
plt.plot(gallazzi05['logmass'],gallazzi05['feh'],linestyle='-',color='black')
h1=plt.scatter(s17['logmstargal'],s17['fehnuc'],color='black',alpha=0.2)
h2=plt.scatter(p11['logmassgal'],p11['fehnuc'],color='r')
plt.scatter(k09['logmassgal'],k09['feh1'],color='r')
plt.scatter(splogmass[1:],spfeh[1:],color='r')
plt.scatter([splogmass[0]],[spfeh[0]],color='red',alpha=1.0,marker="*",s=200)
plt.text(splogmass[0]+0.02,spfeh[0]+0.02,'M54',fontsize=ufontsize/1.4)
plt.legend(labels=['Spec [Fe/H]','Phot [Fe/H]'],handles=[h2,h1],fontsize=ufontsize/1.6)
plt.savefig('mass_metallicity_nsc.pdf',bbox_inches='tight')
plt.show()


plt.figure(figsize=(6,6))
plt.plot([7.8,10.4],[0,0],'--',color='black')
plt.xlim(7.9,10.35)
h1=plt.scatter(s17['logmstargal'],s17['fehnuc']-s17['fehgal'],color='black',alpha=0.2)
h2=plt.scatter(p11['logmassgal'],p11['fehnuc']-p11['fehgal'],color='r')
plt.scatter(k09['logmassgal'],k09['feh1']-k09['fehgal'],color='r')

plt.xlabel('log($M_{\star}$)',fontsize=ufontsize)
plt.ylabel('[Fe/H]$_{nuc}$ - [Fe/H]$_{gal}$',fontsize=ufontsize)
plt.savefig('mass_metallicity_diff.pdf',bbox_inches='tight')

plt.show()

