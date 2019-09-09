#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  9 11:45:04 2019

@author: seth
"""

import numpy as np
from astropy.table import Table, join, vstack
import matplotlib.pyplot as plt

t2=Table.read("lauer05_tab2.tex")
t3=Table.read("lauer05_tab3.tex")
t7=Table.read("lauer05_tab7.tex")
t8=Table.read("lauer05_tab8.tex")
t23=join(t2,t3,keys='galaxy',join_type='left')
t237=join(t23,t7,keys='galaxy',join_type='left')
allt=join(t237,t8,keys='galaxy',join_type='left')
#Use Jordi+ 2006 relation to get V-I to g-i
#vigal is V-I color at 1" radius
allt['g-i']=1.481*allt['vigal']+2.*allt['dvidr']-0.536
#my own fit to I-i=-0.370282*(V-I)-0.161448 -- not a great fit
allt['sloani']=(np.array(allt['MV'],dtype='float')-allt['vigal'])-(-0.370282*allt['vigal']-0.161448)
allt['mli']=10.**(0.979*allt['g-i']-0.831)
allt['logmstar']=np.log10(allt['mli']*10**(-0.4*(allt['sloani']-4.53)))
allt['source']='L05'

#reset NSC flags for only those objects with spectype eq 'a'
emission=(allt['spectype'] != 'a')
allt['nucflag'][emission]=0



nuc=(allt['nucflag'] == 1)
nonuc=(allt['nucflag'] == 0)
plt.scatter(allt['logmstar'][nonuc],allt['g-i'][nonuc],facecolors='none', edgecolors='black')
plt.scatter(allt['logmstar'][nuc],allt['g-i'][nuc],facecolors='black', edgecolors='white',alpha=0.5)
allt.write('lauer05_alltab.fits',overwrite=True)
