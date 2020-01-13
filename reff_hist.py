#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 10:41:45 2019

@author: seth
"""

import numpy as np
from astropy.table import Table, join, vstack
import matplotlib.pyplot as plt


ufontsize=24
#plt.rc('text', usetex=True)
plt.rc('font', family='sans-serif')
plt.rc('xtick', labelsize=ufontsize) 
plt.rc('ytick', labelsize=ufontsize) 
plt.rc('text.latex', preamble=r'\usepackage{cmbright}')


g16=Table.read('georgiev16.fits')
c06=Table.read('cote06_table1.dat',format='ascii')
rh_cote=c06['r_hg']/3600./180.*3.1415*16.52e6
cotegood=(rh_cote > 0.0)
g16good=(g16['e_reffNSC'] < g16['reffNSC']/1.1) #radius is more than 3x the lower error (NOT USING???)




#plt.hist(g16['reffNSC'][early],color='red')
plt.hist(np.log10(g16['reffNSC']),color='blue',histtype='step',linewidth=2)
plt.hist(np.log10(g16['reffNSC'][g16good]),color='blue',alpha=0.5)
#plt.hist(np.log10(rh_cote[cotegood]),color='red',density=True,alpha=0.4)
plt.show()

plt.figure(figsize=(9,6))
plt.xlim(-1,50)
(n,usebins,patches)=plt.hist(g16['reffNSC'][g16good],color='blue',density=True,alpha=0.5,bins=25)
plt.hist(rh_cote[cotegood],color='red',density=True,alpha=0.4,bins=usebins)

plt.xlabel('$r_{eff,NSC}$ [pc]',fontsize=ufontsize)
plt.ylabel('Normalized Number',fontsize=ufontsize)
plt.legend([r'Late-Type (Georgiev+ 2016)',r'Early-Type (C\^ot\'e+ 2006)'],fontsize=ufontsize/1.2)
plt.savefig('reff_histogram.pdf',bbox_inches='tight')
#plt.plot([np.log10(2.),np.log10(2)], [0,1.4], color='red')

plt.show()

allr=np.concatenate((np.array(g16['reffNSC'][g16good]),np.array(rh_cote[cotegood])))
print(np.percentile(allr,(16,50,84)))