# Nuclear Star Cluster Review
## Github repository of data and code used to make figures
### Neumayer, Seth, Boker 2020, Astronomy & Astrophysics Reviews
All code written/figures created by Anil Seth except where otherwise noted.  This page will continue to be updated until final publication.  

### Fig. 1: NSCs in NGC300 and NGC205 -- images and surface brightness profiles.

review_figures.key -- keynote file with images and final SB profiles

plot_sbs.py -- code used to plot sb profiles.  Reads in NGC205_SBP_Valluri_2005.txt from Valluri et al. 2005, additional data using GraphClicks from that paper, ngc205_valluri_clicks.txt.  For NGC300, data from ngc0300_Bland-Hawthorn.txt with data from Bland-Hawthorn et al. 2005, Boeker et al. 2002, and Kim et al. 2004 respectively; the data from both Bland-Hawthorn et al. 2005 and Kim et al. 2004 was digitized.



### Fig. 2: NSCs are Special.  Data on NSCs and GCs from the Virgo cluster.

nscs_are_special.py -- only galaxies with galaxy stellar masses from spengler17_tab8.dat are plotted.  The GC data radii and magnitudes are taken from jordan09.fits for the same galaxies as the nuclei that are plotted.  The right panel contains the same information by organized galaxy-by-galaxy in a mass-ranked order.  Galaxy masses were labeled by hand (not in the code). Figures produced:
nscs_special_radial.pdf (left hand panel)
nscs_special_mass_nox.pdf (right hand panel)


### Fig. 3: Nucleated Galaxy Demographics Figure

galaxy_demographics.py -- primary code used to generate the galaxy color-mass diagram and the occupation fraction diagram.  For details see comments in code.  Produces: galaxy_sample_nscs.pdf, occupation_fraction.pdf

lauer_galaxies.py -- processes tables from Lauer et al. 2005 (lauer05_tab*.tex):
https://ui.adsabs.harvard.edu/abs/2005AJ....129.2138L/abstract
for (mostly) massive elliptical galaxies.  Some are nucleated, but the paper doesn't distinguish between possible AGN and stellar sources.  We consider all galaxies with absorption spectra and detected nuclei to be NSCs.  To put these in context for the demographics figures we need estimates of the galaxy colors.  We transform the V-I colors (reported at 1" and adding the gradient to get the estimated color for the galaxies at 10", typically 500 pc-1 kpc in these galaxies) to g-i colors, and then use the Roediger et al. 2015 relations to obtain M/L estimates.  In total 8/30 galaxies with color estimates have NSCs (after removing ACSCVCS duplicates).  This significantly impacts the number of fraction of high mass ellipticals with reported NSCs.  Additional information on duplicates available in lauer_notes.txt.  



### Fig. 4: NSC Effective Radii histogram

reff_hist.py -- reads in data from Georgiev et al. 2016 (georgiev16.fits) and Cote et al. 2006 (cote06_table1.dat). 
Creates: reff_histogram.pdf


### Fig. 5: Mass Ellipticity diagram

mass_ellipticity.py -- includes data from Geoergiev+ 2014, 2016, Spengler+ 2017 (obtained via private communication; see spengler17_privcomm.dat), and data from Seth+ 2006/2008.  

### Fig. 6: NSC Mass Histogram

mass_histogram.py -- Data from Georgiev+ 2016, Spengler+ 2017, Sanchez-Jannssen+ 2019, and globular cluster data from Harris Milky Way GC catalog, Harris 1996, updated 2010 (mw_gc_cat.fits).

### Fig. 7: NSC Mass/Radius/Density plots

mass_radius_plot.py -- 

### Fig. 8: NGC 247 spectrum

Plot made by Nikolay Kacharov.  See Kacharov+ 2018 for additional details.

### Fig. 9: Mass -- Metallicity Plots for Early Type Galaxies

metallicity_plots.py -- 

### Fig. 10: NSC Rotation Maps

Figures taken from Feldmeier+ 2014, Seth+ 2008b, Seth+ 2010, Lyubenova & Tsatsi 2019. Compiled and edited by the authors. 


### Fig. 11: NSC kinematics

vsig_eps.py -- Data compiled by Nadine Neumayer in v_sig_eps_table.txt.  Also plotted are globular clusater data from Kamann+ 2018 and Bianchini+ 2013 (GCs_Kamann.txt).  

### Fig. 12: NSC Mass -- Galaxy Stellar Mass Relation

mass_scaling.py -- plots as well as fits to data.  

### Fig. 13: BH Mass -- NSC Mass plots

plot_bh_nsc.py -- plots data from many sources compiled in bh_nsc_galmass.csv. 