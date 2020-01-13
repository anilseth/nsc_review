#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 18 12:08:23 2019

@author: seth
"""

import numpy as np
from astropy.table import Table, join, vstack
import matplotlib.pyplot as plt
import matplotlib.colors as colors
#from matplotlib.mlab import bivariate_normal

plt.close('all')

ufontsize=20
plt.rc('text', usetex=True)
plt.rc('font', family='sans-serif')
plt.rc('xtick', labelsize=ufontsize) 
plt.rc('ytick', labelsize=ufontsize) 
plt.rc('text.latex', preamble=r'\usepackage{cmbright}')

bntab=Table.read('bh_nsc_galmass.csv')
dets=(bntab['bhulimit'] == 0) & (bntab['nsculimit'] == 0)
latedets=(bntab['bhulimit'] == 0) & (bntab['nsculimit'] == 0) & (bntab['galtype'] == 1)
earlydets=(bntab['bhulimit'] == 0) & (bntab['nsculimit'] == 0) & (bntab['galtype'] == 2)
ucddets=(bntab['bhulimit'] == 0) & (bntab['nsculimit'] == 0) & (bntab['galtype'] == 0)
latebhu=(bntab['bhulimit'] == 1) & (bntab['galtype'] == 1)
earlybhu=(bntab['bhulimit'] == 1) & (bntab['galtype'] == 2)
latenscu=(bntab['nsculimit'] == 1) & (bntab['galtype'] == 1)
earlynscu=(bntab['nsculimit'] == 1) & (bntab['galtype'] == 2)

plt.figure(figsize=(6,6))
plt.xlim(5.5,8.5)
plt.ylim(-4,4)
plt.xlabel('log($M_{NSC}$)',fontsize=ufontsize)
plt.ylabel('log($M_{BH}/M_{NSC}$)',fontsize=ufontsize)

bntab['xdirect']=0
bntab['ydirect']=-1
plt.plot((5.5,8.5),(0,0),linestyle='--',color='black',alpha=0.3)

plt.quiver(bntab['lognscmass'][latebhu],(bntab['logbhmass']-bntab['lognscmass'])[latebhu],bntab['xdirect'][latebhu],bntab['ydirect'][latebhu],color='b',alpha=0.3)
plt.quiver(bntab['lognscmass'][earlybhu],(bntab['logbhmass']-bntab['lognscmass'])[earlybhu],bntab['xdirect'][earlybhu],bntab['ydirect'][earlybhu],color='r',alpha=0.3)
bntab['xdirectn']=-1
bntab['ydirectn']=1
plt.quiver(bntab['lognscmass'][latenscu],(bntab['logbhmass']-bntab['lognscmass'])[latenscu],bntab['xdirectn'][latenscu],bntab['ydirectn'][latenscu],color='b',alpha=0.3)
plt.quiver(bntab['lognscmass'][earlynscu],(bntab['logbhmass']-bntab['lognscmass'])[earlynscu],bntab['xdirectn'][earlynscu],bntab['ydirectn'][earlynscu],color='r',alpha=0.3)

x1=plt.scatter(bntab['lognscmass'][latedets],(bntab['logbhmass']-bntab['lognscmass'])[latedets],color='b',s=100)
x2=plt.scatter(bntab['lognscmass'][earlydets],(bntab['logbhmass']-bntab['lognscmass'])[earlydets],color='r', s=100)
x3=plt.scatter(bntab['lognscmass'][ucddets],(bntab['logbhmass']-bntab['lognscmass'])[ucddets],color='m',s=100)

plt.legend(labels=['Late-Type','Early-Type','UCDs'],handles=[x1,x2,x3],fontsize=ufontsize/1.2)
plt.savefig('bhnsc_v_mnsc.pdf',bbox_inches='tight')
plt.show()



plt.figure(figsize=(6,6))
plt.xlim(8.6,11.9)
plt.ylim(-4,4)
plt.xlabel('log($M_\star$)',fontsize=ufontsize)
plt.ylabel('log($M_{BH}/M_{NSC}$)',fontsize=ufontsize)
bntab['xdirect']=0
bntab['ydirect']=-1
plt.plot((8.6,11.9),(0,0),linestyle='--',color='black',alpha=0.3)
plt.quiver(bntab['logmstar'][latebhu],(bntab['logbhmass']-bntab['lognscmass'])[latebhu],bntab['xdirect'][latebhu],bntab['ydirect'][latebhu],color='b',alpha=0.3)
plt.quiver(bntab['logmstar'][earlybhu],(bntab['logbhmass']-bntab['lognscmass'])[earlybhu],bntab['xdirect'][earlybhu],bntab['ydirect'][earlybhu],color='r',alpha=0.3)
bntab['xdirectn']=0
bntab['ydirectn']=1
plt.quiver(bntab['logmstar'][latenscu],(bntab['logbhmass']-bntab['lognscmass'])[latenscu],bntab['xdirectn'][latenscu],bntab['ydirectn'][latenscu],color='b',alpha=0.3)
plt.quiver(bntab['logmstar'][earlynscu],(bntab['logbhmass']-bntab['lognscmass'])[earlynscu],bntab['xdirectn'][earlynscu],bntab['ydirectn'][earlynscu],color='r',alpha=0.3)

x1=plt.scatter(bntab['logmstar'][latedets],(bntab['logbhmass']-bntab['lognscmass'])[latedets],color='b',s=100)
x2=plt.scatter(bntab['logmstar'][earlydets],(bntab['logbhmass']-bntab['lognscmass'])[earlydets],color='r', s=100)
#x3=plt.scatter(bntab['lognscmass'][ucddets],(bntab['logbhmass']-bntab['lognscmass'])[ucddets],color='m',s=100)

#plt.legend(labels=['Late-Type','Early-Type'],handles=[x1,x2],fontsize=ufontsize/1.2)
plt.savefig('bhnsc_v_mstar.pdf',bbox_inches='tight')
plt.show()


