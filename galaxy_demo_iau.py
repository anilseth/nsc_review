#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 22 16:26:38 2019

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
nsa=Table.read('/Users/seth/astro/nsa/nsa_v0_1_2.fits')
nsalogmass=np.log10(nsa['MASS'])
nsagi=nsa['ABSMAG'][:,3]-nsa['ABSMAG'][:,5]


#NGFS
#OB18 only has 61 lines, excludes objects where no nuclear data was determined
ob=Table.read('ordenes-briceno18_tab1.tex')
#print("OB18 Columns: ", ob.colnames)
#75 nucleated galaxies, should use this for galaxy data
ei=Table.read('eigenthaler18_tab1.fits')
print("Ei18 Columns: ", ei.colnames)
#Uses Bell 2003 relations for stellar masses
ei['g-i']=ei['__g_-i__0']
ei['sloani']=ei['g_mag_lc']-ei['g-i']
ei['mli']=10.**(0.979*ei['g-i']-0.831)
ei['logmstar']=np.log10(ei['mli']*10**(-0.4*(ei['sloani']-4.53)))
ei['nucflag']=np.zeros(len(ei),dtype=int)
ei['nucflag'][(ei['Type'] == 'o')]=1
ei['source']='NGFS'


#NGVS
sjgal=Table.read('sanchez-janssen19_tab4.tex')
print("sjgal Columns: ", sjgal.colnames)
sjnuc=Table.read('sanchez-janssen19_tab5.tex')
#print("sjnuc Columns: ", sjnuc.colnames)
#Masses determined by Roediger relation fits to all available photometry
sjgal['g-i']=sjgal['gmag']-sjgal['imag']
sjgal['source']='NGVS'
#Georgiev2009 dwarfs
g09=Table.read('georgiev09_galaxies.dat',format='ascii')
print("g09 Columns: ", g09.colnames)
#Plan -- convert V-I and M_V to g-i and i and get Mstar using Roediger?
#Use Jordi+ 2006 relation to get V-I to g-i
g09['g-i']=1.481*g09['vi0']-0.536
#my own fit to I-i=-0.370282*(V-I)-0.161448 -- not a great fit
g09['sloani']=(np.array(g09['mv'],dtype='float')-g09['vi0'])-(-0.370282*g09['vi0']-0.161448)
g09['mli']=10.**(0.979*g09['g-i']-0.831)
g09['logmstar']=np.log10(g09['mli']*10**(-0.4*(g09['sloani']-4.53)))
g09['source']='G09'


#Georgiev2014 sample
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
g14['logmstar']=np.log10(g14['mli']*10**(-0.4*(g14['sloani']-g14['MOD']-4.53)))
g14good=np.isfinite(g14['Imag'])


#ACSVCS Sample
c06=Table.read('cote06_table1.dat',format='ascii')
print("c06 Columns: ", c06.colnames)
f06_tab12=Table.read('ferrarese06_tab12.fits')
#print("f06 Columns: ", f06_tab12.colnames)
f06_tab34=Table.read('ferrarese06_tab34.fits')
#print("f06 Columns: ", f06_tab34.colnames)
f06=join(f06_tab12,f06_tab34,keys='ACSVCS',join_type='left')
print("f06 Columns: ", f06.colnames)

#process f06
f06['g-i']=f06['g-z']*0.7851-0.006-(f06['E_B-V_']*3.1*(1.20585-0.49246))#from CMD website
f06['mlz']=10**(0.886*f06['g-i']-0.848) #from roediger
f06['logmstar']=np.log10(f06['mlz']*10**(-0.4*(f06['zmag_g']-f06['E_B-V_']*3.1*0.49246-31.087-4.50)))
f06['nucflag']=np.zeros(len(f06),dtype=int)
f06['nucflag'][(f06['N'] == 'Ia')]=1
f06['nucflag'][(f06['N'] == 'Ib')]=1
f06['source']='ACSVCS'


allgal=vstack([ei,sjgal,f06,g09,g14],join_type='inner')


logmstar_arr=np.arange(6,12,0.1)
gidivide=(logmstar_arr-6)*0.12+0.40 #rough dividing line blue cloud red sequence
top=np.zeros(len(gidivide))+1.4
bottom=np.zeros(len(gidivide))-1.4

plt.figure(figsize=(6,6))

indnuc=(allgal['nucflag'] == 1)
indnonuc=(allgal['nucflag'] == 0)
#plt.fill_between(logmstar_arr,gidivide,top,color='red',alpha=0.2)
#plt.fill_between(logmstar_arr,bottom,gidivide,color='blue',alpha=0.2)

h1=plt.scatter(allgal['logmstar'][indnonuc],allgal['g-i'][indnonuc],facecolors='none', edgecolors='gray',alpha=1.0)
h2=plt.scatter(allgal['logmstar'][indnuc],allgal['g-i'][indnuc],facecolors='black',alpha=0.5)

#plt.axes.Axes.tick_params(labelsize=16)
plt.xlim(6,12)
plt.ylim(-1.4,1.4)
plt.ylabel('($g-i$)$_0$',fontsize=ufontsize)
plt.xlabel('log($M_\star$)',fontsize=ufontsize)
#plt.legend()

legend=plt.legend(labels=['no NSC detected','w/ NSC'],handles=[h1,h2],fontsize=ufontsize/1.2,frameon=True)
frame = legend.get_frame()
#frame.set_facecolor('#AAAAAA')
frame.set_edgecolor('black')
plt.plot(logmstar_arr,gidivide,color='black',linestyle='--')
#plt.show()
plt.savefig('galaxy_sample_nscs_iau.eps',bbox_inches='tight')

# now work on occupation fraction 
allgal['gitest']=allgal['g-i']-((allgal['logmstar']-6)*0.12+0.40) #if this is positive it is red sequence, negative, blue cloud
redgal=allgal[(allgal['gitest']> 0)]
bluegal=allgal[(allgal['gitest']< 0)]
usebins=np.arange(5.5,12.3,0.7)
(totblue,outbins)=np.histogram(bluegal['logmstar'],bins=usebins)
(nucblue,outbins)=np.histogram(bluegal['logmstar'][(bluegal['nucflag'] == 1)],bins=usebins)
(totred,outbins)=np.histogram(redgal['logmstar'],bins=usebins)
(nucred,outbins)=np.histogram(redgal['logmstar'][(redgal['nucflag'] == 1)],bins=usebins)
occall=((nucblue+nucred)/(totblue+totred))
occblue=nucblue/totblue
occred=nucred/totred
errall=np.sqrt(occall*(1-occall)/(totblue+totred))
errblue=np.sqrt(occblue*(1-occblue)/totblue)
errred=np.sqrt(occred*(1-occred)/totred)

plt.figure(figsize=(6,6))

plotbins=usebins[:-1]+(usebins[1]-usebins[0])*0.5
plt.fill_between(plotbins,occall-errall,occall+errall,color='gray',alpha=0.2)
#plt.fill_between(plotbins,occred-errred,occred+errred,color='red',alpha=0.2)

plt.plot(plotbins,occall,color='black')

#plt.plot(plotbins,occblue,color='blue')
#plt.plot(plotbins,occred,color='red')
#plt.legend(['Late-Type','Early-Type'],fontsize=ufontsize/1.2,loc='lower center')
plt.ylim(0.0,1.0)

plt.ylabel('Fraction of Galaxies with NSC',fontsize=ufontsize)
plt.xlabel('log($M_\star$)',fontsize=ufontsize)
#plt.show()
plt.savefig('occupation_fraction_iau.eps',bbox_inches='tight')

check=allgal[(allgal['logmstar'] > 8.8)]
check2=check[(check['logmstar'] < 9.2)]
check3=check2[(check2['nucflag'] == 0)]

