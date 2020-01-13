#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 20:26:59 2019

@author: seth
"""


import numpy as np
from astropy.table import Table, join, vstack, Column
import matplotlib.pyplot as plt
import pdb

j07=Table.read('jordan09.fits')
s17=Table.read('spengler17_tab8.dat',format='ascii')
c06=Table.read('cote06_table1.dat',format='ascii')

s17.sort('logmstargal')
plt.figure(figsize=(6,6))    

#fig, ax = plt.subplots()
ufontsize=24
plt.gca().invert_yaxis()
plt.xlabel('log($M_{\star}$)',fontsize=ufontsize)
plt.ylabel('Absolute $g$ magnitude',fontsize=ufontsize)
plt.rc('text.latex', preamble=r'\usepackage{cmbright}')

for i in range(0,len(s17)):
    ind=((j07['VCC'] == s17['VCC'][i]) & (j07['gmag']-31.15 > -12))
    cind=(c06['VCC'] == s17['VCC'][i])
    galaxymass=np.zeros(np.sum(ind))+i#s17['logmstargal'][i]
    plt.scatter(galaxymass,j07['gmag'][ind]-31.15,color='black',alpha=0.1)
#    plt.scatter([s17['logmstargal'][i]],[c06['g_AB'][cind]],marker="*",s=200,color='r')
    plt.scatter([i],[c06['g_AB'][cind]-31.15],marker="*",s=300*(2+np.log10(c06['r_hg'][cind])),color='r',alpha=0.5)
    ngcs=np.sum(ind)
    minpos=np.argmin(j07['gmag'][ind])
    mincolor=np.array(j07['gmag'][ind][minpos]-j07['zmag'][ind][minpos])
    dmag=np.min(j07['gmag'][ind]-31.15)-(c06['g_AB'][cind]-31.15)
    ra=np.array(j07['RAJ2000'][ind][minpos])
    dec=np.array(j07['DEJ2000'][ind][minpos])
    print(np.array(s17['VCC'][i]),ngcs,np.array(np.min(j07['gmag'][ind]-31.15)),
          np.array(dmag),np.array(s17['logmstargal'][i]), np.array(c06['r_hg'][cind]*80.))



#plt.ylim(27,15)
#plt.xlim(7.9,10.3)
plt.savefig('nscs_special_mass_nox.pdf',bbox_inches='tight')
plt.show()

plt.figure(figsize=(6,6))    

for i in range(0,len(s17)):
    ind=((j07['VCC'] == s17['VCC'][i]) & (j07['gmag']-31.15 > -12))
    cind=(c06['VCC'] == s17['VCC'][i])
    galaxymass=np.zeros(np.sum(ind))+i#s17['logmstargal'][i]
#    plt.scatter(np.log10(j07['GDist'][ind]*80.),j07['gmag'][ind],color='black',alpha=0.1)
#    plt.scatter([np.log10((c06['offset'][cind]+0.001)*80.)],[c06['g_AB'][cind]],marker="*",s=200,color='r')  
    plt.scatter((j07['GDist'][ind]*80./1000.),j07['gmag'][ind]-31.15,color='black',alpha=0.1)
    plt.scatter([((c06['offset'][cind]+0.001)*80.)/1000.],[c06['g_AB'][cind]-31.15],marker="*",s=300*(2+np.log10(c06['r_hg'][cind])),color='r',alpha=0.5)  
#    print(c06['offset'][cind],np.log10((c06['offset'][cind]+0.001)*80.))
plt.gca().invert_yaxis()
plt.xlabel('Galactocentric Radius [kpc]',fontsize=ufontsize)
plt.ylabel('Absolute $g$ magnitude',fontsize=ufontsize)

plt.legend(['Virgo GCs','Virgo NSCs'],fontsize=ufontsize/1)

plt.savefig('nscs_special_radial.pdf',bbox_inches='tight')
plt.show()