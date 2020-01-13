#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 26 09:08:15 2019

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



vs=Table.read('v_sig_eps_table.txt',format='ascii')
kam=Table.read('GCs_Kamann.txt',format='ascii')




# prep isotropic klines
eps_int=np.arange(0.0,0.8,0.001)
inc=3.1415/2.
ell = np.sqrt(1-(1-eps_int)**2)
Omega_e = (0.5*(np.arcsin(ell)/(np.sqrt(1-ell**2))-ell)) / (ell - np.arcsin(ell)* np.sqrt(1-ell**2))

#Assuming alpha=0.15 (Cappellari 2016)
edgeon_vsig=np.sqrt((Omega_e -1) / (0.15*Omega_e+1))
edgeon_eps = 1- np.sqrt(1+eps_int * (eps_int -2) * np.sin(inc)**2)

allinc=np.arange(0,3.1415/2.,0.01)
useind=60

#Edge-on view: i=90deg
#Projections at different inclinations:
inc_vsig = edgeon_vsig[useind] * np.sin(allinc)  # since anisotropy/delta=0
inc_eps=1- np.sqrt(1+eps_int[useind] * (eps_int[useind] -2) * np.sin(allinc)**2)


fig, ax = plt.subplots(figsize=(7,6))

#plt.figure()


latetype=np.arange(11,13,1)
earlytype=np.append((np.arange(0,11,1)),13)
label=np.array((6,11,12))

ax.plot(edgeon_eps,edgeon_vsig,linestyle='--',color='black')
#plt.plot(inc_eps,inc_vsig,linestyle='--',color='black')

#plt.scatter(vs['epsE'],vs['v_sig_e'])
#ax.scatter(vs['epsE'],vs['v_sig_e'], color='black',s=100)
ldata = ax.scatter(vs['epsE'][latetype],vs['v_sig_e'][latetype], color='blue', s=100)
edata = ax.scatter(vs['epsE'][earlytype],vs['v_sig_e'][earlytype], color='red', s=100)
#data = ax.scatter(vs['epsE'],vs['v_sig_e'], c=vs['log_M'], edgecolors='face', cmap='plasma', s=60)
#cbar = fig.colorbar(data, ax=ax)#,ticks=np.arange(8.5,11.0,0.5))
#bar.set_label('log\,M$_\star$', rotation=270,fontsize=ufontsize,labelpad=20)
ax.scatter(kam['ellip'],kam['vsigma'],facecolors='none',edgecolors='black',s=100)
for i in np.arange(0,3):
    ax.text(vs['epsE'][label][i]-0.01,vs['v_sig_e'][label][i],vs['Name'][label][i],fontsize=ufontsize/1.2,horizontalalignment='right')

ax.set_ylabel('Rotational Support ($V_r / \sigma$)$_e$',fontsize=ufontsize)
ax.set_xlabel('Ellipticity $\epsilon_e$',fontsize=ufontsize)
ax.set_xlim(0,0.6)
ax.set_ylim(0,1)

plt.savefig('vsig_eps.pdf',bbox_inches='tight')


plt.show()