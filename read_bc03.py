#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 10:40:13 2019

@author: seth
"""

import numpy as np
from astropy.table import Table, join, vstack, Column
import matplotlib.pyplot as plt
import pdb


basedir='/Users/seth/programs/bc03/models/Padova2000/chabrier/proc/'
# 0.0001(m22), 0.0004(m32), 0.004(m42), 0.008(m52), 0.02(m62), 0.05(m72)
metaltags=['m122','m132','m142','m152','m162','m172']
zs=[0.0001,0.0004,0.004,0.008,0.02,0.05]



for i in range(len(zs)):
    color2file=basedir+'bc2003_hr_' + metaltags[i] + '_chab_ssp.2color'
    color4file=basedir+'bc2003_hr_' + metaltags[i] + '_chab_ssp.4color'
    color2=Table.read(color2file,format='ascii')
    color4=Table.read(color4file,format='ascii')
    both=join(color2,color4,keys='log-age-yr')
    print(color2file,color2.colnames,color4.colnames)
    both['z']=zs[i]
    if (i == 0):
        final=both
    else:
        final=vstack([final,both],join_type='inner')

final.write('bc03.fits')