#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 09:48:24 2020

@author: tcandela
"""
# =============================================================================
# IMPORTS
# =============================================================================
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import random as rd
from shapely.geometry import LineString
from scipy import interpolate, stats
#Personal Librairies
sys.path.append('/homelocal-px/px-179/tcandela/STAMM/LIB/') #localisation of STAMM librairies
import netCDF_lib as ncl
import plot_lib as pl
import turtle_lib as tul
# =============================================================================
# PATHS & FILES 
# =============================================================================
beach = 'NC_GouaroBay'
year = 2004

indir = '/data/rd_exchange2/tcandela/STAMM/run/PACIFIC/' + beach + '/'
filename_type = '_3y_leatherback_passives_'
filename = beach + filename_type + str(year)
outdir = indir
# =============================================================================
# USER PARAMETERS
# =============================================================================
#MAP PLOT PARAMETERS
#xmin = -100
#xmax = 40
#ymin = -10
#ymax = 60

#MAP PLOT PARAMETERS
xmin = 80
xmax = 220
ymin = -60
ymax = 10

#
lat_space = 10
lon_space = 20

#LETHARGY
coef_SMR = 5
lethargy = 30
To = 24

endday = 183
month = 6
# =============================================================================
# =============================================================================
# =============================================================================
#              /!\/!\/!\/!\/!\/!\ CODE /!\/!\/!\/!\/!\/!\        
# =============================================================================
# =============================================================================
# =============================================================================

# =============================================================================
# PRESENTATION
# =============================================================================

print('======================================================================')
print('=====================COMPUTATION OF TRAJECTORIES======================')
print('======================================================================')

# =============================================================================
# LOADING DATA
# =============================================================================

print('\nLoading data...')
data = ncl.read_nc(indir + filename +'.nc', ['traj_lat','traj_lon','date','traj_temp','init_t','traj_time'])
traj_lat = data['traj_lat']
traj_lon = data['traj_lon']
traj_temp = data['traj_temp']
init_t = data['init_t']
date = data['date']
traj_time = data['traj_time']

nb_turtles = traj_lat.shape[1]
nb_days = traj_lat.shape[0]

print('\nLoaded file : ' + indir + filename +'.nc')

#%% =============================================================================
# SEPARATION OF LIVING AND DEAD TURTLES
# =============================================================================

print('\nSeparation of living and dead turtles...')
days = traj_temp.shape[0]
age_year = np.arange(days)/365.
index_dead = []
index_alive = []

date_death = tul.find_date_death(nb_turtles, traj_temp, To, coef_SMR, lethargy, init_t, nb_days)

turtle = 0
for turtle in np.arange(nb_turtles):
    if date_death[turtle] == 1.0000e+34:
        index_alive.append(turtle)
    else:
        index_dead.append(turtle)
print('\n' + str(len(index_alive)) + ' alive turtles and ' + str(len(index_dead)) + ' dead turtles has been spotted !')

# =============================================================================
# PEAK TRAJECTORIES
# =============================================================================

print('\nComputation of trajectories for alive turtles...')
ifile =  indir + filename + '.nc'   
ofile = indir + filename + '_alive_dispersion_' + str(month) + 'm.png'

list_turtle = []
for turtle in np.arange(nb_turtles):
#    if turtle in index_alive:
    list_turtle.append(turtle)


lon = np.zeros((nb_days, len(list_turtle)))
lat = np.zeros((nb_days, len(list_turtle)))
time = np.zeros((nb_days, len(list_turtle)))
i = 0
for turtle in np.arange(nb_turtles):
    if turtle in list_turtle:            
        for day in np.arange(nb_days):
            lon[day, i] = traj_lon[day, turtle]
            lat[day, i] = traj_lat[day, turtle]            
            time[day, i] = traj_time[day, turtle]
        i += 1

# Plot figure ------
c = 0.89
f = plt.figure(figsize = (12*c/2.54,8*c/2.54))
gs = gridspec.GridSpec(2,1,height_ratios=[11,1],left=0.08, right=0.98, bottom=0.07, top=0.95)
ax = plt.subplot(gs[0])
im,time = pl.display_trajectories_particular(lon[:endday,:], lat[:endday,:], time[:endday,:], xmin, f,ax)
pl.show_start_point(ax, lat, lon)
pl.plot_map(ax,ymin,ymax,xmin,xmax,lon_space=lon_space,lat_space=lat_space)
ax.spines['right'].set_linewidth(0.5)
ax.spines['left'].set_linewidth(0.5)
ax.spines['bottom'].set_linewidth(0.5)
ax.spines['top'].set_linewidth(0.5)
    

if np.max(time)>0*365: 
    ax_cb = plt.subplot(gs[1])
    label = u"Age (Days)"
    pl.display_colorbar(f,im, ax_cb, label)

else:
    ax_cb = plt.subplot(gs[2])
    label = u"Age (Days)"
    pl.display_colorbar(f,im, ax_cb, label)
        
    
plt.savefig(ofile,bbox_inches='tight',dpi=800)
#plt.show()
print('\nPlot Saved at : ' + ofile)


















