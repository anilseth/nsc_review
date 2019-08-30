#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 10 16:08:34 2019

@author: seth
"""

import numpy as np
from astropy.table import Table, join, vstack, Column
import matplotlib.pyplot as plt
import pdb


plt.close('all')

ufontsize=20
plt.rc('text', usetex=True)
plt.rc('font', family='sans-serif')
plt.rc('xtick', labelsize=ufontsize) 
plt.rc('ytick', labelsize=ufontsize) 


#NGVS
sjgal=Table.read('sanchez-janssen19_tab4.tex')
print("sjgal Columns: ", sjgal.colnames)
sjnuc=Table.read('sanchez-janssen19_tab5.tex')
#print("sjnuc Columns: ", sjnuc.colnames)
#Masses determined by Roediger relation fits to all available photometry
sjgal['g-i']=sjgal['gmag']-sjgal['imag']
sjgal['source']='NGVS'
sj=join(sjnuc,sjgal,keys='id')
sj['logmnuc']=sj['logmstar_1']
sj['logmstargal']=sj['logmstar_2']
sj['T']=-1

#NGFS
#OB18 only has 61 lines, excludes objects where no nuclear data was determined
ob=Table.read('ordenes-briceno18_tab1.tex')
ob['NGFS'] = [x[4:-1] for x in ob['nucleus']]


#print("OB18 Columns: ", ob.colnames)
#75 nucleated galaxies, should use this for galaxy data
ei=Table.read('eigenthaler18_tab1.fits')
print("Ei18 Columns: ", ei.colnames)
#Uses Bell 2003 relations for stellar masses
ei['g-i']=ei['__g_-i__0']
ei['sloani']=ei['g_mag_lc']-ei['g-i']
ei['mli']=10.**(0.979*ei['g-i']-0.831)
ei['logmstargal']=np.log10(ei['mli']*10**(-0.4*(ei['sloani']-4.53)))
ei['nucflag']=np.zeros(len(ei),dtype=int)
ei['nucflag'][(ei['Type'] == 'o')]=1
ei['source']='NGFS'
ngfs=join(ob,ei,keys='NGFS')
ngfs['T']=-1
ngfs['logmnuc'][2]='0.00       '  #fix one missing mass
ngfs['logmnuc']=np.array([x[:5] for x in ngfs['logmnuc']],dtype='float')

#Georgiev2014/16 sample
g16=Table.read('georgiev16.fits')
g16['Name']=g16['Gal']
g14_gal=Table.read('georgiev14_tab1plus.fits')
print("g14_gal Columns: ", g14_gal.colnames)
#168 of 228 galaxies with Imag
g14_nonuc=Table.read('georgiev14_tab2.fits')
print("g14_nonuc Columns: ", g14_nonuc.colnames)
#71 of 95 with Imag
g14_gal['nucflag']=1
g14_nonuc['nucflag']=0
#remove all sources that don't have NoNSC in it
#test=[i for i in g14['Com'] if 'NoNSC' not in i]

g14=vstack([g14_gal,g14_nonuc],join_type='outer')
g14['source']='G14'

#prep Georgiev14
g14['g-i']=0.7451*(g14['Bmag']-g14['Imag'])-0.4388 #very good fit using Padova SSPs
g14['sloani']=g14['Imag']+0.3887+0.0831*(g14['Bmag']-g14['Imag'])  #poor fit using Padova SSPs, max residual ~0.1 mags.  
#Use roediger relation to get to g-i
g14['mli']=10.**(0.979*g14['g-i']-0.831)
g14['logmstargal']=np.log10(g14['mli']*10**(-0.4*(g14['sloani']-g14['MOD']-4.53)))
g14good=np.isfinite(g14['Imag'])
g14use=g14[g14good]
ge=join(g14use,g16,keys='Name')
ge['logmnuc']=np.log10(ge['MNSC']*1e4)
ge['T']=ge['t']

#ACSVCS Sample
s17=Table.read('spengler17_tab8.dat',format='ascii')
s17['source']="S17"
s17['logmnuc']=s17['logmstarnuc']
s17['T']=-1

e12=Table.read('erwin12_tab2.tex')
e12['source']="E12"

allnuc=vstack([sj,ngfs,ge,s17,e12],join_type='inner')
allnuc['dmass']=allnuc['logmnuc']-allnuc['logmstargal']


plt.figure(figsize=(6,6))
earlyind=(allnuc['T'] < 1)
lateind=(allnuc['T'] > 1)
plt.scatter(allnuc['logmstargal'][lateind],allnuc['logmnuc'][lateind],color='b',alpha=0.5)
plt.scatter(allnuc['logmstargal'][earlyind],allnuc['logmnuc'][earlyind],color='r',alpha=0.5)
e12earlyind=(e12['T'] < 1)
e12lateind=(e12['T'] > 1)
plt.scatter(e12['logmstargal'][e12lateind],e12['logmnuc'][e12lateind],color='b',alpha=1.0,marker="*",s=200)
plt.scatter(e12['logmstargal'][e12earlyind],e12['logmnuc'][e12earlyind],color='r',alpha=1.0,marker="*",s=200)

plt.xlim(5.5,11.2)
plt.ylim(4.5,9)
plt.xlabel('log($M_{\star}$)',fontsize=ufontsize)
plt.ylabel('log($M_{NSC}$)',fontsize=ufontsize)
plt.legend(['Late-Type','Early-Type','Good Masses'],fontsize=ufontsize/1.4)
plt.savefig('mnuc_mgal.pdf',bbox_inches='tight')
plt.show()



late=allnuc[lateind]
early=allnuc[earlyind]
usebins=np.arange(5.6,12.3,0.4)
plotbins=usebins[:-1]+0.25
latepercentiles=np.zeros((3,len(usebins)-1))
earlypercentiles=np.zeros((3,len(usebins)-1))
for i in range(len(usebins)-1):
    einbin=early[np.where((early['logmstargal'] > usebins[i]) & (early['logmstargal'] < usebins[i+1]))]
    linbin=late[np.where((late['logmstargal'] > usebins[i]) & (late['logmstargal'] < usebins[i+1]))]
    if len(einbin) > 5:
#        pdb.set_trace() 
        earlypercentiles[:,i]=np.percentile(einbin['dmass'],[25,50,75])
    if len(linbin) > 5:
        latepercentiles[:,i]=np.percentile(linbin['dmass'],[25,50,75])
        
    #    earlypercentiles[:,i]=np.percentile(earlyboth,[10,50,90])
#    latelow=late[(late['logmstargal'] > usebins[i])]
#    lateboth=latelow[(latelow['logmstargal'] < usebins[i+1])]
#    latepercentiles[:,i]=np.percentile(lateboth,[10,50,90])
lind=np.where(latepercentiles[1,:] < 0)
eind=np.where(earlypercentiles[1,:] < 0)
print(plotbins[eind])
print(plotbins[lind])

plt.figure(figsize=(6,6))
plt.xlim(5.5,11.2)
plt.ylim(-4,0)
plt.xlabel('log($M_{\star}$)',fontsize=ufontsize)
plt.ylabel('log($M_{NSC}/M_{\star}$)',fontsize=ufontsize)

#latepercentiles=10.**latepercentiles
#earlypercentiles=10.**earlypercentiles

plt.fill_between(plotbins[lind],latepercentiles[0][lind],latepercentiles[2][lind],color='blue',alpha=0.2)
plt.fill_between(plotbins[eind],earlypercentiles[0][eind],earlypercentiles[2][eind],color='red',alpha=0.2)
h1=plt.plot(plotbins[lind],latepercentiles[1][lind],color='b')
h2=plt.plot(plotbins[eind],earlypercentiles[1][eind],color='r')
plt.legend(['Late-Type','Early-Type'],fontsize=ufontsize/1.2)
#legend=plt.legend(labels=['no NSC detected','w/ NSC'],handles=[h1,h2],fontsize=ufontsize/1.2,frameon=True)
plt.savefig('frac_mnsc_mgal.pdf',bbox_inches='tight')

plt.show()