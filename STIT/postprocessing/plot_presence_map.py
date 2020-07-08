#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
USE: python  plot_presence_map.py  output.nc
"""
# =============================================================================
# IMPORTS
# =============================================================================
import sys, os
import netCDF4 as nc

#Personal Librairies
sys.path.insert(1, os.path.join(sys.path[0], '../../LIB'))
import netCDF_lib as ncl
import turtle_lib as tul
import plot_lib as pl
# =============================================================================
# PATHS & FILES
# =============================================================================
file = sys.argv[1]
outdir = ncl.get_directory(file)
filename = ncl.get_name(file)
imageName = outdir+filename+'_alive_density_.png' 

# =============================================================================
# USERS PARAMETERS
# =============================================================================
dx = 1 #boxes size, in degrees
seuil_numpos = 25000
lat_space = 10
lon_space = 20
   
#MAP PLOT PARAMETERS
   
##### NATL
#xmin = -100
#xmax = 40
#ymin = -10
#ymax = 60

##### GOM
#xmin = -99.5
#xmax = -79.5
#ymin = 17.5
#ymax = 32.5

#PAC
#xmin = -60
#xmax = 300
#ymin = -60
#ymax = 60

#
xmin = 80
xmax = 220
ymin = -60
ymax = 10
   
coef_SMR= 5.
To = 24
lethargy = 30
alive = False
# =============================================================================
# CODE
# =============================================================================    
print('======================================================================')
print('=====================COMPUTATION OF PRESENCE MAP======================')
print('======================================================================')

mode = 'Density'
ncfile = nc.Dataset(file)
dico = ncl.read_nc(file, ['traj_lon','traj_lat'])
traj_lon = dico['traj_lon']
traj_lat = dico['traj_lat']

duration = traj_lon.shape[0]
turtles = traj_lon.shape[1]  
 
x_area,y_area, temp, t_init= ncl.extract_x_y(ncfile, duration, turtles, alive)
  
if alive == True : 
    print('ALIVE')
    x_alive, y_alive = tul.remove_dead_turtles(x_area,y_area,temp, t_init,To,lethargy,coef_SMR, duration)
    heatmap,extent,cmap = pl.compute_presence_map(dx,x_alive,y_alive,xmin,xmax,ymin,ymax,seuil_numpos,lat_space,lon_space)      
else : 
    print('ALL')
    heatmap,extent,cmap = pl.compute_presence_map(dx,x_area,y_area,xmin,xmax,ymin,ymax,seuil_numpos,lat_space,lon_space)
    
pl.plot_presence_map(heatmap,extent,cmap,xmin,ymin,xmax,ymax,lat_space,lon_space,imageName, mode, traj_lat, traj_lon,dx,seuil_numpos)
