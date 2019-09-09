# Nuclear Star Cluster Review Data and Code to make Figures
## Neumayer, Seth, Boker 2019, Astronomy & Astrophysics Reviews
All code written/figures created by Anil Seth.

## Nucleated Galaxy Demographics Figure

galaxy_demographics.py -- primary code used to generate the galaxy color-mass diagram and the 

lauer_galaxies.py -- processes tables from Lauer et al. 2005 (lauer05_tab*.tex):
https://ui.adsabs.harvard.edu/abs/2005AJ....129.2138L/abstract
for (mostly) massive elliptical galaxies.  Some are nucleated, but the paper doesn't distinguish between possible AGN and stellar sources.  We consider all galaxies with absorption spectra and detected nuclei to be NSCs.  To put these in context for the demographics figures we need estimates of the galaxy colors.  We transform the V-I colors (reported at 1" and adding the gradient to get the estimated color for the galaxies at 10", typically 500 pc-1 kpc in these galaxies) to g-i colors, and then use the Roediger et al. 2015 relations to obtain M/L estimates.  In total 8/40 galaxies with color estimates have NSCs.  This significantly impacts the number of fraction of high mass ellipticals with reported NSCs.  