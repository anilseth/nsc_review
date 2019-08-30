#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 21 06:29:21 2019

@author: seth
"""

import numpy as np
from astropy.table import Table, join, vstack
import matplotlib.pyplot as plt

sb205=Table.read('NGC205_SBP_Valluri_2005.txt',format='ascii')
vc205=Table.read('ngc205_valluri_clicks.txt',format='ascii')
bh300=Table.read('ngc0300_Bland-Hawthorn.txt',format='ascii')
bo300=Table.read('ngc0300_Boeker.txt',format='ascii')
ki300=Table.read('ngc0300_kim2004.txt',format='ascii')

plt.rcParams['font.family'] = "sans-serif"




ufontsize=19
fontProperties = {'family':'sans-serif','sans-serif':['Helvetica'],
                  'weight' : 'normal', 'size' : ufontsize}
#a = plt.gca()
#a.set_xticklabels(a.get_xticks(), fontProperties)
#a.set_yticklabels(a.get_yticks(), fontProperties)
plt.rc('text', usetex=False)
plt.rc('font', family='sans-serif')
plt.rc('xtick', labelsize=ufontsize) 
plt.rc('ytick', labelsize=ufontsize) 
plt.rc('text.latex', preamble=r'\usepackage{cmbright}')

#first do NGC300
indbh=bh300['r'] > 300
indki=ki300['r'] > 10

ax1=plt.subplot(1,1,1)
ax1.scatter(bh300['r'][indbh],bh300['sb'][indbh],facecolors='none', edgecolors='black')
ax1.scatter(ki300['r'][indki][::3],ki300['sb'][indki][::3],facecolors='none', edgecolors='black')
ax1.scatter(bo300['r'],bo300['sb']+0.58,facecolors='none', edgecolors='black')

radius=np.arange(0.1,2.e3,0.1)
sl=135. 
expdisk=-2.5*np.log10(np.exp(-radius/sl))+2.5*np.log10(np.exp(-np.min(radius)/sl))+19.97
ax1.plot(radius,expdisk,linestyle='--',color='black')

ax1.set_xscale('log')
ax1.set_xlim(0.1,1500)
ax1.set_ylim(31.5,14)

ax1.set_xlabel('Radius ["]',fontsize=ufontsize*1.2)
ax1.set_ylabel('Surface Brightness',fontsize=ufontsize*1.2)
#ax1.set_yticks(fontname = "Helvetica")

distance=2.15 #Karachentsev 2013
newlabel = [1,10,100,1000,10000] # labels of the xticklabels: the position in the new x-axis
arclabel= np.array(newlabel)/2.15e6*3600*180./3.1415
#k2degc = lambda t: t-273.15 # convert function: from Kelvin to Degree Celsius
#arcsectodist = lambda t: t-273.15 # convert function: from Kelvin to Degree Celsius


ax2 = ax1.twiny()
ax2.set_xscale('log')
#newpos   = [k2degc(t) for t in newlabel]   # position of the xticklabels in the old #x-axis
ax2.set_xticks(arclabel)
ax2.set_xticklabels(newlabel)
ax2.set_xlabel('Radius [pc]',fontsize=ufontsize*1.2)
ax2.set_xlim(ax1.get_xlim())
plt.savefig('sb300.pdf',bbox_inches='tight')

plt.show()

vcind=(vc205['r'] < 2)
vcind2=(vc205['r'] > 2)

ax1=plt.subplot(1,1,1)
radius=np.arange(0.02,300,0.02)
re=120*2.2 #re from Nguyen+ 2018
n=1.4
bn=1.9992*n-0.3271
mue=21 # guessing to fit
sersic=mue+2.5*bn/np.log(10)*((radius/re)**(1/n)-1)
ax1.plot(radius,sersic,linestyle='--',color='black')
ax1.scatter(sb205['r'],sb205['m_I'],facecolors='none', edgecolors='black')
ax1.scatter(vc205['r'][vcind][::2],vc205['m_I'][vcind][::2],facecolors='none', edgecolors='black')
ax1.scatter(vc205['r'][vcind2],vc205['m_I'][vcind2],facecolors='none', edgecolors='black')

ax1.set_xscale('log')
ax1.set_xlim(0.02,300)
ax1.set_ylim(22,12)
ax1.set_xlabel('Radius ["]',fontsize=ufontsize*1.2)
ax1.set_ylabel('Surface Brightness',fontsize=ufontsize*1.2)
distance=0.82e6 #Karachentsev 2013
newlabel = [1,10,100,1000] # labels of the xticklabels: the position in the new x-axis
arclabel= np.array(newlabel)/distance*3600*180./3.1415
#k2degc = lambda t: t-273.15 # convert function: from Kelvin to Degree Celsius
#arcsectodist = lambda t: t-273.15 # convert function: from Kelvin to Degree Celsius


ax2 = ax1.twiny()
ax2.set_xscale('log')
#newpos   = [k2degc(t) for t in newlabel]   # position of the xticklabels in the old #x-axis
ax2.set_xticks(arclabel)
ax2.set_xticklabels(newlabel)
ax2.set_xlabel('Radius [pc]',fontsize=ufontsize*1.2)
ax2.set_xlim(ax1.get_xlim())
plt.savefig('sb205.pdf',bbox_inches='tight')

plt.show()
