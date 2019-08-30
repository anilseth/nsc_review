#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 11:03:44 2019

@author: seth
"""

import numpy as np
from astropy.table import Table, join, vstack
import matplotlib.pyplot as plt

plt.close('all')

ufontsize=24
plt.rc('text', usetex=True)
plt.rc('font', family='sans-serif')
plt.rc('xtick', labelsize=ufontsize) 
plt.rc('ytick', labelsize=ufontsize) 


g16=Table.read('georgiev16.fits')
s17=Table.read('spengler17_tab8.dat',format='ascii')
sjgal=Table.read('sanchez-janssen19_tab4.tex')
sjnuc=Table.read('sanchez-janssen19_tab5.tex')
sj=join(sjnuc,sjgal,keys='id')
harris=Table.read('mw_gc_cat.fits')
#ob=Table.read('ordenes-briceno18_tab1.tex') not easily used, no masses
logmnuc=np.concatenate((np.array(sj['logmstar_1']),np.array(s17['logmstarnuc']),np.log10(g16['MNSC']*1.e4)))
logmgal=np.concatenate((np.array(sj['logmstar_2']),np.array(s17['logmstargal']),np.log10(g16['Mgal_']*1.e7)))

highmass=(logmgal > 9.0)
#earlylogm=np.concatenate((np.array(sjnuc['logmstar']),np.array(s17['logmstarnuc'],np.log10(g16['MNSC']*1.e4))))


harrismv=np.array([i for i in harris['M_V,t'] if '-' in i],dtype=float)
harrislogmass=np.log10(10**(-0.4*(harrismv-4.8))*2.)

      

plt.figure(figsize=(9,6))
plt.xlim(4,9)

weights=np.zeros(len(logmnuc))+1./len(logmnuc)*(np.max(logmnuc)-np.min(logmnuc))
plt.hist(harrislogmass,color='black',alpha=0.3,bins=12,density=True)
#plt.hist(np.log10(g16['MNSC']*1.e4),color='b',alpha=0.5,density=True,bins=12)
(thishist,usebins,usepatches)=plt.hist(logmnuc,color='g',alpha=1,bins=15,histtype='step',weights=weights,linewidth=2)
plt.hist(logmnuc[highmass],color='g',alpha=0.4,bins=usebins,weights=weights[highmass])

plt.xlabel('log($M_{NSC,GC}$)',fontsize=ufontsize)
plt.ylabel('Normalized Number',fontsize=ufontsize)
plt.legend(labels=['All NSCs','Milky Way GCs','NSCs in $M_\star > 10^9$ M$_\odot$'],fontsize=ufontsize/1.4)
plt.savefig('mass_histogram.pdf',bbox_inches='tight')
plt.show()
