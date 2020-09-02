#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Create initial position file from a .json file and a grid.
USE:   python  create_init_pos.py  myfile1.json  myfile2.json
"""
# =============================================================================
# IMPORTS
# =============================================================================
import numpy as np
import os
import sys
import datetime
import glob
import ntpath
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid.inset_locator import inset_axes

#personal libraries
sys.path.insert(1, os.path.join(sys.path[0], '../../LIB'))
import librumeau as brum
import netCDF_lib as ncl
import init_pos_lib as ipl
# =============================================================================
# PATHS & FILES
# =============================================================================
#gridfle parameters
gridfile = '/homelocal-px/px-171/pgiffard/DATA/PSY4_12_DAYLY/grid.nc'
lon_name = 'longitude'
lat_name = 'latitude'
mask_name = 'mask'
grid_type = 'regular' #orca or regular
#figure plot does not work with orca option, need to change gridtogeo

#infile & outfile
infile = sys.argv[1]
list_files = sys.argv[1:]
path = ncl.get_directory(infile)
outfile = infile.replace(".json","_initial_positions.txt")
append = True 

#figure
savefig_path = path + 'Release_Map.png'

'''
To annotate country names and change map scale, go to part named "PLOT"
'''
# =============================================================================
# USERS PARAMETERS
# =============================================================================
mode = 'lonlat' #lonlat or xy
disk_radius = 0.05 #for disk mode, in degrees
width = 1.852 #km, for rectangle mode
beach_orientation = 'E' #'W' or 'E', for beach mode and square mode.
release_zone_size = 0.25 #for square mode, release zone side size in degrees

#total of released turtles
nturtles = 1500

##zone to plot               #
xmin = -81              #
xmax = -78              #  
ymin = 26               #
ymax = 27.5              #
#                            #
lat_space = 1.5             #
lon_space = 3.0             # PACUARE
#                            #
##dot on the minimap         #
xdot = -80.0                #
ydot = 26.8                #

##zone to plot               #
#xmin = -89.0                #
#xmax = -79.0                #  
#ymin = 8.0                  #
#ymax = 14.5                 #
#                            #
#lat_space = 1.5             # PLAYA
#lon_space = 3.0             # GRANDE
#                            #
##dot on the minimap         #
#xdot = -85.8                #
#ydot = 10.3                 #

##zone to plot               #
#xmin = -83.0                #
#xmax = -73.0                #  
#ymin = 7.0                  #
#ymax = 13.5                 #
#                            # GULF
#lat_space = 1.5             # OF
#lon_space = 3.0             # URABA
#                            #
##dot on the minimap         #
#xdot = -77.0                #
#ydot = 8.50                 #

##zone to plot               #
#xmin = -57.0                #
#xmax = -47.0                #  
#ymin = 2.0                  #
#ymax = 8.5                  #
#                            #
#lat_space = 1.5             #
#lon_space = 3.0             # AWALA
#                            #
##dot on the minimap         #
#xdot = -57.0                #
#ydot = 5.0                  #

##zone to plot               #
#xmin = -86.0                #
#xmax = -76.0                #  
#ymin = 7.5                  #
#ymax = 14.0                 #
#                            #
#lat_space = 1.5             # CHIRIQUI
#lon_space = 3.0             # BEACH
#                            #
##dot on the minimap         #
#xdot = -81.6                #
#ydot = 8.95                 #

##zone to plot               #
#xmin = -70.0                #
#xmax = -60.0                #  
#ymin = 8.5                  #
#ymax = 15.0                 #
#                            #
#lat_space = 1.5             # GRANDE
#lon_space = 3.0             # RIVIERE
#                            #
##dot on the minimap         #
#xdot = -61.0                #
#ydot = 10.8                 #
#
##zone to plot               #
#xmin = -84.0                #
#xmax = -74.0                #  
#ymin = 23.0                 #
#ymax = 29.5                 #
#                            #
#lat_space = 1.5             # JUNO
#lon_space = 3.0             # BEACH
#                            #
##dot on the minimap         #
#xdot = -80.0                #
#ydot = 26.8                 #

##zone to plot               #
#xmin = -70.0                #
#xmax = -60.0                #  
#ymin = 15.0                 #
#ymax = 21.5                 #
#                            #
#lat_space = 1.5             # SAINT
#lon_space = 3.0             # CROIX
#                            #
##dot on the minimap         #
#xdot = -64.9                #
#ydot = 17.67                #

##zone to plot               #
#xmin = -90.0                #
#xmax = -80.0                #  
#ymin = 13                   #
#ymax = 19.5                 #
#                            #
#lat_space = 1.5             # CABO
#lon_space = 3.0             # CAMARON
#                            for square mode, #
##dot on the minimap         #
#xdot = -85.0                #
#ydot = 16.0                 #

#zone to plot               #
#xmin = 64.0                 #
#xmax = 74.0                 #  
#ymin = -8.0                 #
#ymax = -1.5                 #
                            #
#lat_space = 1.5             # PEROS
#lon_space = 3.0             # BANHOS
                            #
#dot on the minimap         #
#xdot = 71.7                 #
#ydot = -5.4                 #
# =============================================================================
# CODE
# =============================================================================
d_min = 40 - (111*release_zone_size)/2
d_max = 40 + (111*release_zone_size)/2

nb_year=1
nesting_year = 2020
date_ref = datetime.datetime(nesting_year,1,1) #

coord_mode = 'grid' #Anna:?
outfile = path + outfile

if os.path.exists(outfile) :
    if append == False:
        raise ValueError("The output file already exist, if you want to append, put append in True mode\n")
    else :
        print("\n The output file already exist and will be appended.\n")
    
file = open (outfile, 'a+') 
    
#lecture du fichier contenant les informations sur la grille pour la conversion (lon,lat) -> (x,y)
grid = ncl.read_nc(gridfile,[lat_name,lon_name,mask_name])

if grid_type == 'orca':
    lon_mat = np.squeeze(grid[lon_name])
    lat_mat = np.squeeze(grid[lat_name])

elif grid_type == 'regular':   
    try:
        lon_mat = np.squeeze(grid[lon_name])[0,:]
        lat_mat = np.squeeze(grid[lat_name])[:,0]
    except IndexError:
        lon_mat = grid[lon_name] 
        lat_mat = grid[lat_name]
else:
    raise ValueError('Please set grid_type to orca or regular')

if len(np.squeeze(grid[mask_name]).shape) == 3:
    griddata = np.squeeze(grid[mask_name])[0,:,:]
else:
    griddata =  np.squeeze(grid[mask_name])


            
for f in list_files :
    i=0
    beach_r = ipl.beach_json(release_zone_size, lon_name, lat_name, ipl.read_beach_json(f, grid, lon_name, lat_name), nesting_year, d_min, d_max, date_ref)
    print("Creating intial positions for beach: "+beach_r.beach_name)
    E = ipl.echantillon(beach_r)
    E.gen_positions(beach_r.nb_turtles, nb_year, lat_mat, lon_mat, griddata, coord_mode, grid_type, beach_orientation, disk_radius, width)
    if mode == 'lonlat':
        for i in np.arange(beach_r.nb_turtles):
            E.positions[i].x, E.positions[i].y = brum.grid_to_geo(E.positions[i].x, E.positions[i].y, lon_mat, lat_mat)
            if E.positions[i].x > 180:
                E.positions[i].x = E.positions[i].x-360
    E.write_all(file)
    filename = ntpath.basename(f).split('.')[0]+'.png'
#    E.histogram(filename)
#    E.release_map(filename,lon_mat,lat_mat)
#%% =============================================================================
# PLOT    
# =============================================================================
#Zone to plot in minimap
xmin_sl = -180
xmax_sl =  180
ymin_sl = -80
ymax_sl =  85

#yticks parameters
ymax_pos = ymax
dy = 0.1
import time
time.sleep(1)        
#Loading initial positions
initfile = open(outfile,'r')
print(outfile)
x_init, y_init, t_init = np.loadtxt(initfile,usecols=(0,1,3),unpack=True)
x_init, y_init, t_init = x_init[:nturtles], y_init[:nturtles], t_init[:nturtles]
print('\nInitial positions loaded')

if mode == 'lonlat':
    lon, lat = x_init, y_init


#Convert grid point into lon/lat
#if mode == 'xy':
#    lon, lat = brum.grid_to_geo(x_init, y_init, lon_mat, lat_mat)
#    for i in np.arange(len(lon)):
#        if lon[i] > 180:
#            lon[i] = lon[i]-360
    
if mode == 'xy':
    lon = np.zeros(len(x_init))
    lat = np.zeros(len(x_init))
    for i in np.arange(len(x_init)):
        lon[i], lat[i] = brum.grid_to_geo(x_init[i], y_init[i], lon_mat, lat_mat)
        if lon[i] > 180:
            lon[i] = lon[i]-360
            
# PLOT CODE  
#Figure parameters
fig = plt.figure(figsize=(7/2.54,5/2.54))
ax_pos = fig.add_subplot(111)
ax_pos_small = inset_axes(ax_pos, width="50%", height="25%", loc='upper right')

ax_pos.set_yticks(np.arange(ymin,ymax_pos,dy))
ax_pos.spines['right'].set_linewidth(0.5)
ax_pos.spines['left'].set_linewidth(0.5)
ax_pos.spines['bottom'].set_linewidth(0.5)
ax_pos.spines['top'].set_linewidth(0.5)
    
#Plot map    
m = Basemap(ax=ax_pos, projection='merc',lon_0=5,lat_0=90.0,llcrnrlon=xmin,urcrnrlon=xmax,llcrnrlat=ymin,urcrnrlat=ymax,resolution='h')
m.fillcontinents(color='0.65',alpha=1, lake_color='w')
m.drawcoastlines(color='0.3',linewidth=0.2)
m.drawcountries(color='w',linewidth=0.1)
m.drawparallels(np.arange(ymin,ymax_pos,lat_space), labels=[1,0,0,0], fontsize=5, linewidth=0.1)
m.drawmeridians(np.arange(xmin,xmax,lon_space), labels=[0,0,0,1], fontsize=5, linewidth=0.1)
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''    
# =============================================================================
# BEGINNING OF PART TO CHANGE
# =============================================================================
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
#Plot annotations
##PACUARE#
#ax_pos.annotate(u'Nicaragua', xy=(0.20, 0.83), fontsize=5,fontweight='normal',xycoords='axes fraction' ,horizontalalignment='right', verticalalignment='center',color='k')
#ax_pos.annotate(u'Costa Rica', xy=(0.295, 0.35), fontsize=5,fontweight='normal',xycoords='axes fraction' ,horizontalalignment='right', verticalalignment='center',color='k')
#ax_pos.annotate(u'Panama', xy=(0.50, 0.16), fontsize=5,fontweight='normal',xycoords='axes fraction' ,horizontalalignment='right', verticalalignment='center',color='k')

##PLAYA GRANDE#
#ax_pos.annotate(u'Nicaragua', xy=(0.50, 0.80), fontsize=5,fontweight='normal',xycoords='axes fraction' ,horizontalalignment='right', verticalalignment='center',color='k')
#ax_pos.annotate(u'Costa Rica', xy=(0.55, 0.35), fontsize=5,fontweight='normal',xycoords='axes fraction' ,horizontalalignment='right', verticalalignment='center',color='k')
#ax_pos.annotate(u'Panama', xy=(0.83, 0.08), fontsize=5,fontweight='normal',xycoords='axes fraction' ,horizontalalignment='right', verticalalignment='center',color='k')
#ax_pos.annotate(u'Honduras', xy=(0.25, 0.95), fontsize=5,fontweight='normal',xycoords='axes fraction' ,horizontalalignment='right', verticalalignment='center',color='k')
#ax_pos.annotate(u'Guatemala', xy=(0.11, 0.85), fontsize=3,fontweight='normal',xycoords='axes fraction' ,horizontalalignment='right', verticalalignment='center',color='k')

##GULF OF URABA#
#ax_pos.annotate(u'Costa Rica', xy=(0.25, 0.25), fontsize=5,fontweight='normal',xycoords='axes fraction' ,horizontalalignment='right', verticalalignment='center',color='k')
#ax_pos.annotate(u'Panama', xy=(0.475, 0.33), fontsize=5,fontweight='normal',xycoords='axes fraction' ,horizontalalignment='right', verticalalignment='center',color='k')
#ax_pos.annotate(u'Colombia', xy=(0.85, 0.16), fontsize=5,fontweight='normal',xycoords='axes fraction' ,horizontalalignment='right', verticalalignment='center',color='k')

##AWALA-YALIMAPO / CAYENNE#
#ax_pos.annotate(u'Suriname', xy=(0.22, 0.30), fontsize=5,fontweight='normal',xycoords='axes fraction' ,horizontalalignment='right', verticalalignment='center',color='k')
#ax_pos.annotate(u'French Guiana', xy=(0.50, 0.33), fontsize=5,fontweight='normal',xycoords='axes fraction' ,horizontalalignment='right', verticalalignment='center',color='k')
#ax_pos.annotate(u'Brazil', xy=(0.575, 0.12), fontsize=5,fontweight='normal',xycoords='axes fraction' ,horizontalalignment='right', verticalalignment='center',color='k')

##CHIRIQUI BEACH#
#ax_pos.annotate(u'Nicaragua', xy=(0.20, 0.83), fontsize=5,fontweight='normal',xycoords='axes fraction' ,horizontalalignment='right', verticalalignment='center',color='k')
#ax_pos.annotate(u'Costa Rica', xy=(0.295, 0.35), fontsize=5,fontweight='normal',xycoords='axes fraction' ,horizontalalignment='right', verticalalignment='center',color='k')
#ax_pos.annotate(u'Panama', xy=(0.50, 0.16), fontsize=5,fontweight='normal',xycoords='axes fraction' ,horizontalalignment='right', verticalalignment='center',color='k')

##GRANDE RIVIERE#
#ax_pos.annotate(u'Trinidad \n& \nTobago', xy=(0.875, 0.3125), fontsize=2,fontweight='normal',xycoords='axes fraction' ,horizontalalignment='center', verticalalignment='center',color='k')
#ax_pos.annotate(u'Venezuela', xy=(0.50, 0.10), fontsize=5,fontweight='normal',xycoords='axes fraction' ,horizontalalignment='right', verticalalignment='center',color='k')

##JUNO BEACH#
#ax_pos.annotate(u'USA', xy=(0.30, 0.5), fontsize=5,fontweight='normal',xycoords='axes fraction' ,horizontalalignment='center', verticalalignment='center',color='k')
#ax_pos.annotate(u'Bahamas', xy=(0.75, 0.335), fontsize=5,fontweight='normal',xycoords='axes fraction' ,horizontalalignment='right', verticalalignment='center',color='k')

##SAINT CROIX#
#ax_pos.annotate(u'Dominican \nRepublic', xy=(0.0745, 0.57), fontsize=4,fontweight='normal',xycoords='axes fraction' ,horizontalalignment='center', verticalalignment='center',color='k')
#ax_pos.annotate(u'Porto Rico', xy=(0.427, 0.5), fontsize=4,fontweight='normal',xycoords='axes fraction' ,horizontalalignment='right', verticalalignment='center',color='k')
#ax_pos.annotate(u'U.S. \nVirgin \nIslands', xy=(0.56, 0.37), fontsize=3,fontweight='normal',xycoords='axes fraction' ,horizontalalignment='center', verticalalignment='center',color='k')

##CABO CAMARON#
#ax_pos.annotate(u'Honduras', xy=(0.55, 0.30), fontsize=5,fontweight='normal',xycoords='axes fraction' ,horizontalalignment='right', verticalalignment='center',color='k')
#ax_pos.annotate(u'El Salvador', xy=(0.205, 0.10), fontsize=4,fontweight='normal',xycoords='axes fraction' ,horizontalalignment='right', verticalalignment='center',color='k')
#ax_pos.annotate(u'Belize', xy=(0.15, 0.60), fontsize=3,fontweight='normal',xycoords='axes fraction' ,horizontalalignment='right', verticalalignment='center',color='k')
#ax_pos.annotate(u'Mexico', xy=(0.15, 0.875), fontsize=5,fontweight='normal',xycoords='axes fraction' ,horizontalalignment='right', verticalalignment='center',color='k')
#ax_pos.annotate(u'Guatemala', xy=(0.105, 0.34), fontsize=3,fontweight='normal',xycoords='axes fraction' ,horizontalalignment='right', verticalalignment='center',color='k')

#PEROS BANHOS#
ax_pos.annotate(u'Chagos \nArchipelago', xy=(0.79, 0.28), fontsize=4,fontweight='normal',xycoords='axes fraction' ,horizontalalignment='center', verticalalignment='center',color='k')
ax_pos.annotate(u'Diego \nGarcia', xy=(0.81, 0.11), fontsize=2,fontweight='normal',xycoords='axes fraction' ,horizontalalignment='center', verticalalignment='center',color='k')

''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
# =============================================================================
# END OF PART TO CHANGE
# =============================================================================
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
#PLot turtles
lon_m,lat_m = m(lon,lat)
m.scatter(lon_m,lat_m,color='b', marker="o",edgecolors='None',s = 0.1)#,alpha=0.8,marker=".")
#m.plot([lon0,lon1],[lat0,lat1])
#m.scatter([lon],[lat],color='b', marker="o",edgecolors='None',s = 10)#,alpha=0.8,marker=".")
m_pt = Basemap(ax=ax_pos_small, projection='cyl',lon_0=-180,lat_0=90.0,llcrnrlon=xmin_sl,urcrnrlon=xmax_sl,llcrnrlat=ymin_sl,urcrnrlat=ymax_sl,resolution='l')
m_pt.fillcontinents(color='0.65',alpha=1, lake_color='w')
m_pt.scatter([xdot],[ydot], color='b', marker="o", edgecolors='None', s=4,zorder=12)
    
#Save figure
plt.subplots_adjust(left=0.1, right=0.95, bottom=0.02, top=0.98)
plt.savefig(savefig_path,dpi = 800)
   
print('\nMission accomplie ! (plot saved at ' + savefig_path +')\n')
