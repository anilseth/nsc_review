#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 25 16:14:17 2019

@author: seth
"""

import numpy as np
from astropy.table import Table, join, vstack
import matplotlib.pyplot as plt

ufontsize=24
fontProperties = {'family':'sans-serif','sans-serif':['Helvetica'],
                  'weight' : 'normal', 'size' : ufontsize}
#a = plt.gca()
#a.set_xticklabels(a.get_xticks(), fontProperties)
#a.set_yticklabels(a.get_yticks(), fontProperties)
plt.rc('text', usetex=True)
plt.rc('font', family='sans-serif')
plt.rc('xtick', labelsize=ufontsize) 
plt.rc('ytick', labelsize=ufontsize) 
plt.rc('text.latex', preamble=r'\usepackage{cmbright}')




sp=Table.read('spengler17_privcomm.dat',format='ascii')
extra=Table.read('mass_ellip_extra.dat',format='ascii')

g16=Table.read('georgiev16.fits')
g16['Name']=g16['Gal']
g14full=Table.read('ge14ge16_joint.fits',2)
g16['Ell']=g14full['Ell']
good=((g16['reffNSC'] > g16['e_reffNSC']) & (g14full['MOD'] < 31.5)) #radius > than error

g16['oute']=np.nan

#filters to use 1/606, 2/814, 4/555, 6/439
#cut between 
t1=((np.abs(g14full['pa2']-g14full['pa4']) < 20.) & (g14full['e2']/g14full['e4'] > 0.6666) & (g14full['e2']/g14full['e4'] < 1.5) & (g16['reffNSC'] > g16['e_reffNSC']))
g16['oute'][t1]=(g14full['e2'][t1]+g14full['e4'][t1])/2.


#t2=((np.abs(g14full['pa1']-g14full['pa2']) < 20.) & (g14full['e1']/g14full['e2'] > 0.6666) & (g14full['e1']/g14full['e2'] < 1.5) & (g16['reffNSC'] > g16['e_reffNSC']))

t2=((np.abs(g14full['pa1']-g14full['pa2']) < 20.) & (g14full['e1']/g14full['e2'] > 0.6666) & (g14full['e1']/g14full['e2'] < 1.5) & (g16['reffNSC'] > g16['e_reffNSC']))

g16['oute'][t2]=(g14full['e1'][t2]+g14full['e2'][t2])/2.

plt.figure(figsize=(9,6))
plt.scatter(sp['nscmass'],1-sp['ba'],color='red', s=80*(0.3+np.log10(sp['re'])),alpha=0.5)
plt.scatter(np.log10(g16['MNSC']*1.e4),g16['oute'],color='blue',s=80*(0.3+np.log10(g16['reffNSC'])),alpha=0.5)
plt.scatter(extra['ncmass'],extra['ell'],color='blue',s=80*(0.3+np.log10(extra['reff'])),alpha=0.5,marker='s')


plt.legend(['Early-Type','Late-Type','Edge-On Late-Type'],fontsize=ufontsize/1.4)

plt.xlabel('log(NSC Mass)',fontsize=ufontsize)
plt.ylabel('Ellipticity',fontsize=ufontsize)
plt.savefig('mass_ellipticity.pdf',bbox_inches='tight')
plt.show()

print("Minimum and Maximum Sizes")
t=(np.isfinite(g16['oute']))
print(np.min(g16['reffNSC'][t]),np.max(sp['re']))



plt.scatter(sp['galmass'],1-sp['ba'],color='red', s=80*(0.3+np.log10(sp['re'])),alpha=0.5)
plt.scatter(np.log10(g16['Mgal_']*1.e7),g16['oute'],color='blue',s=80*(0.3+np.log10(g16['reffNSC'])),alpha=0.5)
plt.xlabel('log(Galaxy Mass)',fontsize=ufontsize)
plt.ylabel('Ellipticity',fontsize=ufontsize)
plt.show()

#Table 2
lowearly=(sp['galmass'] < 9.0)
highearly=(sp['galmass'] >= 9.0)
print(np.percentile(1-sp['ba'][lowearly],(16,50,84)))
print(np.percentile(1-sp['ba'][highearly],(16,50,84)))

lowlate=((g16['Mgal_'] < 100) & (np.isfinite(g16['oute'])))
highlate=((g16['Mgal_'] > 100) & (np.isfinite(g16['oute'])))
print(np.percentile(g16['oute'][lowlate],(16,50,84)))
print(np.percentile(g16['oute'][highlate],(16,50,84)))


#plt.scatter(g14full['Ell'][good],g14full['e2'][good])
#plt.scatter(g14full['Ell'][good],g14full['e1'][good])
#plt.scatter(g14full['Ell'][good],g14full['e4'][good])
#plt.xlabel('Galaxy Ellipticity',fontsize=ufontsize)
#plt.ylabel('NSC Ellipticity',fontsize=ufontsize)
#plt.show()
