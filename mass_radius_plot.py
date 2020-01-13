    #!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 13 12:03:39 2019

@author: seth
"""

import numpy as np
from astropy.table import Table, join, vstack
import matplotlib.pyplot as plt
import matplotlib.colors as colors
#from matplotlib.mlab import bivariate_normal


#nsc data first
s17=Table.read('spengler17_tab8.dat',format='ascii')
g16=Table.read('georgiev16.fits')
c06=Table.read('cote06_table1.dat',format='ascii')
acsvcs=join(s17,c06,keys='VCC',join_type='inner')
acsvcs['logmnuc']=acsvcs['logmstarnuc']
oe12=Table.read('erwin12_tab2.tex')
oe12['source']="E12"
ag=Table.read('additional_goodmass.dat',format='ascii')
ag['source']='N18'
e12=vstack([oe12,ag],join_type='inner')
#E12 now combines Nguyen+ 2018 and Erwin & Gadotti 2012

mylen = np.vectorize(len)

g16['nfilter_NSC']=mylen(g16['FNSC'])
g16['dist_pc']=10.0**(g16['m-M']/5+1)
g16good=np.where((g16['e_reffNSC'] < g16['reffNSC']/3) & (g16['nfilter_NSC'] > 1) & (g16['reffNSC'] > g16['dist_pc']*0.02/3600./180*3.1415)) #radius is more than 3x the lower error (NOT USING???)
g16['logsigmae']=np.log10(g16['MNSC']*1e4/2./(3.1415*g16['reffNSC']**2))
g16['logmnuc']=np.log10(g16['MNSC']*1e4)
g16['logmstargal']=np.log10(g16['Mgal_']*1.e7)
g16['reff']=g16['reffNSC']

#non NSC data
#nsa=Table.read('/Users/seth/astro/nsa/nsa_v0_1_2.fits')
#nsalogmass=np.log10(nsa['MASS'])
#h0=70. #km/s/Mpc
#c=2.99e5 # km/s
#nsasizepc=np.log10(nsa['SERSIC_TH50']/3600./180*3.1415*nsa['ZDIST']*c/h0*1.e6)

#harris=Table.read('harris10_mwgc.fits')
#harris['logmass']=np.log10(10**(-0.4*(harris['ABS_VMAG']-4.8))*2.)
#harris['logreff']=np.log10(harris['HELIO_DISTANCE']*1000.*harris['HALF_LIGHT_RADIUS']/60./180.*3.1415)
#harrisrh=np.array([i for i in harris['rh'] if '.' in i],dtype=float)
#harrisdist=np.array([i for i in harris['R_sun'] if '.' in i],dtype=float)
n14=Table.read('norris14_taba1.fits')
n14['logmass']=np.log10(n14['M_'])
n14['logre']=np.log10(n14['Re'])
n14['logsigmae']=np.log10(n14['M_']/2./(3.1415*n14['Re']**2))

n14gal=(n14['Type'] != 3)
n14co=np.where((n14['Type'] != 3))# & (n14['Type'] < 7))


e12['logsigmae']=np.log10(10**(e12['logmnuc'])/2./(3.1415*e12['reff']**2))
acsvcs['logsigmae']=np.log10(10**(acsvcs['logmstarnuc'])/2./(3.1415*(acsvcs['r_hg']*80.)**2))
acsvcs['reff']=acsvcs['r_hg']*80.
acsvcs['T']=-5
g16['T']=5
allnuc=vstack([e12,g16[g16good],acsvcs],join_type='inner')



plt.close('all')

ufontsize=24
plt.rc('text', usetex=True)
plt.rc('font', family='sans-serif')
plt.rc('xtick', labelsize=ufontsize) 
plt.rc('ytick', labelsize=ufontsize) 
plt.rc('text.latex', preamble=r'\usepackage{cmbright}')

#plt.hexbin(nsalogmass,nsasizepc,cmap=plt.cm.Greens,mincnt=50,gridsize=50)
#good=np.where((harris['logmass'] > 4.5) & (harris['logreff'] > 0))
#plt.hexbin(harris['logmass'][good],harris['logreff'][good],cmap=plt.cm.Purples,mincnt=1,gridsize=10)#,,mincnt=50)


#mass-radius small
plt.figure(figsize=(6,6))

h4=plt.hexbin(n14['logmass'][n14co],n14['logre'][n14co],cmap=plt.cm.Greys,mincnt=0,gridsize=(60,30),vmax=5, vmin=0.1,alpha=0.5)
#plt.hexbin(n14['logmass'][n14gal],n14['logre'][n14gal],cmap=plt.cm.Reds,mincnt=1,gridsize=(30,10))
plt.plot([7.5,9],[-0.3,0.2],'--',color='black')
plt.plot([7.5,9],[-0.1,0.65],'--',color='black')
plt.text(8.1,-0.05,'$r_{eff} \propto M_{NSC}^{1/3}$',fontsize=ufontsize/1.4,rotation=18)
plt.text(8.2,0.9,'$r_{eff} \propto M_{NSC}^{1/2}$',fontsize=ufontsize/1.4,rotation=30)
h1=plt.scatter(np.log10(g16['MNSC'][g16good]*1e4),np.log10(g16['reffNSC'][g16good]),color='b',alpha=0.5)
h2=plt.scatter(acsvcs['logmstarnuc'],np.log10(acsvcs['r_hg']*80.),color='r',alpha=0.5)
e12earlyind=(e12['T'] < 1)
e12lateind=(e12['T'] > 1)
h3=plt.scatter(e12['logmnuc'][e12lateind],np.log10(e12['reff'][e12lateind]),color='b',alpha=0.5,marker="*",s=200)
#plt.legend(labels=['Late-Type','Early-Type','Good Mass'],handles=[h1,h2,h3],fontsize=ufontsize/1.4,loc=2)
plt.scatter(e12['logmnuc'][e12earlyind],np.log10(e12['reff'][e12earlyind]),color='r',alpha=0.5,marker="*",s=200)
#plt.legend(labels=['Late-type','Early-type','Good Mass','GCs/UCDs/cEs'],handles=[h1,h2,h3,h4])
plt.xlim(5,9.5)
plt.ylim(-0.4,3.5)
plt.xlabel('log(Object Mass)',fontsize=ufontsize)
plt.ylabel('log($r_{eff}$ [pc])',fontsize=ufontsize)
plt.savefig('mass_radius.pdf',bbox_inches='tight')
plt.show()



#mass-sigmae small
plt.figure(figsize=(6,6))
h4=plt.hexbin(n14['logmass'][n14co],n14['logsigmae'][n14co],cmap=plt.cm.Greys,mincnt=0,gridsize=(60,30),vmax=5, vmin=0.1,alpha=0.5)
h1=plt.scatter(np.log10(g16['MNSC'][g16good]*1e4),g16['logsigmae'][g16good],color='b',alpha=0.5)
h2=plt.scatter(acsvcs['logmstarnuc'],acsvcs['logsigmae'],color='r',alpha=0.5)
e12earlyind=(e12['T'] < 1)
e12lateind=(e12['T'] > 1)
h3=plt.scatter(e12['logmnuc'][e12lateind],e12['logsigmae'][e12lateind],color='b',alpha=0.5,marker="*",s=200)
#plt.legend(labels=['Late-Type','Early-Type','Good Mass'],handles=[h1,h2,h3],fontsize=ufontsize/1.4,loc=4)
plt.scatter(e12['logmnuc'][e12earlyind],e12['logsigmae'][e12earlyind],color='r',alpha=0.5,marker="*",s=200)
plt.xlim(5,9.5)
#plt.ylim(-0.4,3.5)
plt.xlabel('log(Object Mass)',fontsize=ufontsize)
plt.ylabel('log($\Sigma_{eff}$ [M$_\odot$/pc$^2$])',fontsize=ufontsize)
plt.savefig('mass_density.pdf',bbox_inches='tight')
plt.show()


#mass-radius big
plt.figure(figsize=(6,6))
#plt.hexbin(nsalogmass,nsasizepc,cmap=plt.cm.Greens,mincnt=50,gridsize=50)
h4=plt.hexbin(n14['logmass'][n14co],n14['logre'][n14co],cmap=plt.cm.Greys,mincnt=0,gridsize=(60,30),vmax=5, vmin=0.1,alpha=0.5)
h1=plt.scatter(np.log10(g16['MNSC'][g16good]*1e4),np.log10(g16['reffNSC'][g16good]),color='b',alpha=0.5)
h2=plt.scatter(acsvcs['logmstarnuc'],np.log10(acsvcs['r_hg']*80.),color='r',alpha=0.5)
e12earlyind=(e12['T'] < 1)
e12lateind=(e12['T'] > 1)
h3=plt.scatter(e12['logmnuc'][e12lateind],np.log10(e12['reff'][e12lateind]),color='b',alpha=0.5,marker="*",s=200)
plt.legend(labels=['Late-Type NSCs','Early-Type NSCs','Dyn/Spec Masses'],handles=[h1,h2,h3],fontsize=ufontsize/1.2,loc=2)
plt.scatter(e12['logmnuc'][e12earlyind],np.log10(e12['reff'][e12earlyind]),color='r',alpha=0.5,marker="*",s=200)
#plt.legend(labels=['Late-type','Early-type','Good Mass','GCs/UCDs/cEs'],handles=[h1,h2,h3,h4])
plt.xlim(4.5,12)
plt.ylim(-0.4,5)
plt.xlabel('log(Object Mass)',fontsize=ufontsize)
plt.ylabel('log($r_{eff}$ [pc])',fontsize=ufontsize)
plt.savefig('mass_radius_wide.pdf',bbox_inches='tight')
plt.show()



#mass-sigmae big
plt.figure(figsize=(6,6))
h4=plt.hexbin(n14['logmass'][n14co],n14['logsigmae'][n14co],cmap=plt.cm.Greys,mincnt=0,gridsize=(60,30),vmax=5, vmin=0.1,alpha=0.7)
h1=plt.scatter(np.log10(g16['MNSC'][g16good]*1e4),g16['logsigmae'][g16good],color='b',alpha=0.5)
h2=plt.scatter(acsvcs['logmstarnuc'],acsvcs['logsigmae'],color='r',alpha=0.5)
e12earlyind=(e12['T'] < 1)
e12lateind=(e12['T'] > 1)
h3=plt.scatter(e12['logmnuc'][e12lateind],e12['logsigmae'][e12lateind],color='b',alpha=0.5,marker="*",s=200)
#plt.legend(labels=['Late-Type','Early-Type','Good Mass'],handles=[h1,h2,h3],fontsize=ufontsize/1.4,loc=4)
plt.scatter(e12['logmnuc'][e12earlyind],e12['logsigmae'][e12earlyind],color='r',alpha=0.5,marker="*",s=200)
plt.xlim(4.5,12)
#plt.ylim(-0.4,3.5)
plt.xlabel('log(Object Mass)',fontsize=ufontsize)
plt.ylabel('log($\Sigma_{eff}$ [M$_\odot$/pc$^2$])',fontsize=ufontsize)
plt.savefig('mass_density_wide.pdf',bbox_inches='tight')
plt.show()


#for Table 2
earlylow=((allnuc['T'] <= 0) & (allnuc['logmstargal'] < 9.0))
earlyhigh=((allnuc['T'] <= 0) & (allnuc['logmstargal'] >= 9.0))
print("Early Low ",np.percentile(allnuc['reff'][earlylow],(16,50,84)))
print("Early High ",np.percentile(allnuc['reff'][earlyhigh],(16,50,84)))
print("Early Low ",np.percentile(allnuc['logsigmae'][earlylow],(16,50,84)))
print("Early High ",np.percentile(allnuc['logsigmae'][earlyhigh],(16,50,84)))

late=(allnuc['T'] > 0)

latelow=((allnuc['T'] > 0) & (allnuc['logmstargal'] < 9.0))
latehigh=((allnuc['T'] > 0) & (allnuc['logmstargal'] >= 9.0))
print("Late Low ",np.percentile(allnuc['reff'][latelow],(16,50,84)))
print("Late High ",np.percentile(allnuc['reff'][latehigh],(16,50,84)))
print("Late Low ",np.percentile(allnuc['logsigmae'][latelow],(16,50,84)))
print("Late High ",np.percentile(allnuc['logsigmae'][latehigh],(16,50,84)))


#Plot Milky Way???