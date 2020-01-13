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
from scipy.stats import norm


plt.close('all')

ufontsize=20
plt.rc('text', usetex=True)
plt.rc('font', family='sans-serif')
plt.rc('xtick', labelsize=ufontsize) 
plt.rc('ytick', labelsize=ufontsize) 
plt.rc('text.latex', preamble=r'\usepackage{cmbright}')


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

oe12=Table.read('erwin12_tab2.tex')
oe12['source']="E12"
ag=Table.read('additional_goodmass.dat',format='ascii')
ag['source']='N18'
e12=vstack([oe12,ag],join_type='inner')

l05=Table.read('lauer05_alltab.fits')
l05nuc=(l05['nucflag'] == 1)
l05=l05[l05nuc]

allnuc=vstack([sj,ngfs,ge,s17,e12,l05],join_type='inner')
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
plt.legend(['Late-Type','Early-Type','Dyn/Spec Masses'],fontsize=ufontsize/1.4)

fit=np.polyfit(allnuc['logmstargal']-9,allnuc['logmnuc'],1)
xarr=np.arange(5.5,11.2,0.1)
yarr=fit[0]*(xarr-9)+fit[1]
yfit=fit[0]*(allnuc['logmstargal']-9)+fit[1]
scatter=yfit-allnuc['logmnuc']
print('SCATTER',norm.fit(scatter))
plt.plot(xarr,yarr,linestyle='dashed',color='black',linewidth=2)

#highmass=(allnuc['logmstargal'] > 9)
#fithigh=np.polyfit(allnuc['logmstargal'][highmass],allnuc['logmnuc'][highmass],1)
#xarr=np.arange(9.0,11.2,0.1)
#yarr=fithigh[0]*xarr+fithigh[1]
#plt.plot(xarr,yarr,linestyle='dashed',color='cyan')

#highmass=(allnuc['logmstargal'] > 9)
fitgood=np.polyfit((e12['logmstargal']-9),e12['logmnuc'],1)
xarr=np.arange(9.0,11.2,0.1)
yarr=fitgood[0]*(xarr-9)+fitgood[1]
plt.plot(xarr,yarr,linestyle='dashdot',color='black',linewidth=2)

plt.savefig('mnuc_mgal.pdf',bbox_inches='tight')
plt.show()


fitarr=np.zeros((2,100))
fitgoodarr=np.zeros((2,100))
for i in range(0,100):
    resample=np.random.randint(0,high=len(allnuc),size=len(allnuc))
    fitarr[:,i]=np.polyfit(allnuc['logmstargal'][resample]-9,allnuc['logmnuc'][resample],1)
    regood=np.random.randint(0,high=len(e12),size=len(e12))
    fitgoodarr[:,i]=np.polyfit((e12['logmstargal'][regood]-9),e12['logmnuc'][regood],1)
   
print("Fit Errors",np.std(fitarr[0,:]),np.std(fitarr[1,:]))    

print("Fit FGood Errors",np.std(fitgoodarr[0,:]),np.std(fitgoodarr[1,:]))    



late=allnuc[lateind]
early=allnuc[earlyind]
usebins=np.arange(5.5,12.3,0.5)
plotbins=usebins[:-1]+0.25
latepercentiles=np.zeros((3,len(usebins)-1))
earlypercentiles=np.zeros((3,len(usebins)-1))
for i in range(len(usebins)-1):
    einbin=early[np.where((early['logmstargal'] > usebins[i]) & (early['logmstargal'] < usebins[i+1]))]
    linbin=late[np.where((late['logmstargal'] > usebins[i]) & (late['logmstargal'] < usebins[i+1]))]
    if len(einbin) > 4:
#        pdb.set_trace() 
        earlypercentiles[:,i]=np.percentile(einbin['dmass'],[25,50,75])
    if len(linbin) > 4:
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


xarr=np.arange(5.5,11.2,0.1)
yarr=fit[0]*(xarr-9)+fit[1]-xarr
plt.plot(xarr,yarr,linestyle='dashed',color='black',linewidth=2)
#xarr=np.arange(9.0,11.2,0.1)
#yarr=fitgood[0]*xarr+fitgood[1]-xarr
#plt.plot(xarr,yarr,linestyle='dashdot',color='black',linewidth=2)
plt.savefig('frac_mnsc_mgal.pdf',bbox_inches='tight')


plt.show()


#Table 2
earlylow=((allnuc['T'] <= 0) & (allnuc['logmstargal'] < 9.0))
earlyhigh=((allnuc['T'] <= 0) & (allnuc['logmstargal'] >= 9.0))
print("Early Low ",np.percentile(allnuc['logmnuc'][earlylow],(16,50,84)))
print("Early High ",np.percentile(allnuc['logmnuc'][earlyhigh],(16,50,84)))
print("Early Low ",10.**np.percentile(allnuc['dmass'][earlylow],(16,50,84)))
print("Early High ",10.**np.percentile(allnuc['dmass'][earlyhigh],(16,50,84)))

late=(allnuc['T'] > 0)

latelow=((allnuc['T'] > 0) & (allnuc['logmstargal'] < 9.0))
latehigh=((allnuc['T'] > 0) & (allnuc['logmstargal'] >= 9.0))
print("Late Low ",np.percentile(allnuc['logmnuc'][latelow],(16,50,84)))
print("Late High ",np.percentile(allnuc['logmnuc'][latehigh],(16,50,84)))
print("Late Low ",10.**np.percentile(allnuc['dmass'][latelow],(16,50,84)))
print("Late High ",10.**np.percentile(allnuc['dmass'][latehigh],(16,50,84)))
