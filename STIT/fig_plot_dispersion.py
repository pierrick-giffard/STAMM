#!/usr/bin/env python
#-*- coding:utf-8 -*-
# =============================================================================
# IMPORTS
# =============================================================================
import os
import sys
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
from matplotlib import gridspec

#Personal Librairies
sys.path.append('/data/rd_exchange2/tcandela/STAMM/LIB/') #localisation of STAMM librairies
import plot_lib as pl
import netCDF_lib as ncl
# =============================================================================
# USERS PARAMETERS
# =============================================================================
#Definition de la zone d'affichage
##NATL
#xmin = -110
#xmax = 30
#ymin = -10
#ymax = 70

#PAC
#xmin = -60
#xmax = 300
#ymin = -60
#ymax = 60

#PAC
xmin = 0
xmax = 120
ymin = -60
ymax = 20


lat_space = 20
lon_space = 40
# =============================================================================
# CLASS
# =============================================================================

# =============================================================================
# CODE
# =============================================================================
ifile = sys.argv[1]  
ofile = ifile.replace(".nc","_dispersion.png")
## -- Si l'affichage coupe les trace à 0° de longitude, bien vérifier que x n'est pas corrigé dans 
## -- PlotTrajectories2(display_trajectories_all) : #x = np.array([l+(((1-np.sign(l))/2)*360) for l in x])
## -- il faut que cette ligne soit bien commentée

# Plot figure ------
c = 0.89
f = plt.figure(figsize = (12*c/2.54,8*c/2.54))
gs = gridspec.GridSpec(2,1,height_ratios=[11,1],left=0.08, right=0.98, bottom=0.07, top=0.95)

dataFile = ncl.data(ifile)
#
zoom_out = 10
xmin = np.min(dataFile.lon) - zoom_out
xmax = np.max(dataFile.lon) + zoom_out
ymin = np.min(dataFile.lat) - zoom_out
ymax = np.max(dataFile.lat) + zoom_out
#
ax = plt.subplot(gs[0])
     
im,time = pl.display_trajectories(dataFile,f,ax,xmin)
pl.show_start_point(ax, dataFile.lat, dataFile.lon-360)
pl.plot_map(ax,ymin,ymax,xmin,xmax,lon_space=lon_space,lat_space=lat_space)

 
ax.spines['right'].set_linewidth(0.5)
ax.spines['left'].set_linewidth(0.5)
ax.spines['bottom'].set_linewidth(0.5)
ax.spines['top'].set_linewidth(0.5)
        

if np.max(time)>0*365: 
    ax_cb = plt.subplot(gs[1])
    label = u"Age (Years)"
    pl.display_colorbar(f,im, ax_cb, label)

else:
    ax_cb = plt.subplot(gs[2])
    label = u"Age (Days)"
    pl.display_colorbar(f,im, ax_cb, label)
            
        
plt.savefig(ofile,bbox_inches='tight',dpi=800)
print("wrote {}".format(ofile))
plt.close()
