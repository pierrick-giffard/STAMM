#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Animate STAMM output files.
Choose parameters.
Execute as:
python plot_animation.py output.nc namelist
"""
# =============================================================================
# IMPORTS
# =============================================================================
import sys, os
from pathlib import Path

#Personal librairies
sys.path.insert(1, os.path.join(sys.path[0], '../../LIB'))
import netCDF_lib as ncl
import plot_lib as pl
import IOlib as IO


# =============================================================================
# USERS PARAMETERS
# =============================================================================
#background
hab_mode = 'current' # 'food', 'temp', 'tot', 'current', 'void'
#option 'void' for no background,
mortality = True #set to False not to calculate dead turtles
jeanette = True
nb_cat = 10
colors = ['red', 'darkgreen', 'blue', 'yellow', 'magenta', 'cyan', 'darkorange', 'black', 'lawngreen', 'darkviolet']
#zone: 'Atlantic, 'Caribbean', 'Pacific', 'Indian', 'Gulf_stream', 'Acores', 'Med', 'COR', 'tmp'
zone = 'test'

# time delta between 2 frames (in days)
h = 1

#Plot first last_turtle. -1 to plot all turtles
last_turtle = 1000

#Dates
start_day = 0
end_day = 365

#gridfile for NPP lon/lat
gridfile = '/data/rd_exchange2/tcandela/STAMM/ressources/VGPM/VGPM_083_mesh.nc'

# Video
fps = 8 #images/sec
dpi = 350 #images resolution
#
tracer = "PP" #PP or mnk
# =============================================================================
# PATHS & FILES
# =============================================================================
file_path = sys.argv[1]
namelist = sys.argv[2]
#Defaults save paths
directory = ncl.get_directory(file_path)
filename = ncl.get_name(file_path)
save_path = directory + 'animation_' + filename + '/'
if not Path(save_path).exists():
    Path(save_path).mkdir(parents=True)
#
videofile = save_path + '/' + filename + '_animation.avi'
#   
param = IO.read_namelist(namelist, display=False)

print('\n')
print('********************************************************************************')
print("Writing frames in ", save_path)
print('********************************************************************************')
print('\n')



# =============================================================================
# HABITAT & MORTALITY PARAMETERS
# =============================================================================
#Température optimale pour le calcul d'habitat
To = 22.

# Temps de léthargie (nombre de jours maximum dans Tw<Tmin)
lethargy=10.

# Variables pour le calcul de l'habitat thermique et alimentaire
coef_SMR=5.




# =============================================================================
#ZONES
# =============================================================================
if zone == 'Atlantic':
    lonmin = -110
    lonmax = 36
    latmin = 0
    latmax = 62

elif zone == 'Caribbean':
    lonmin = 260.
    lonmax = 305.
    latmin = 7
    latmax = 31

elif zone == 'Pacific':
    lonmin = 175
    lonmax = 285
    latmin = -0
    latmax = 40

elif zone == 'Indian':
    lonmin = 0
    lonmax = 120
    latmin = -60
    latmax = 20
    
elif zone == 'Gulf_stream':
    lonmin = -85
    lonmax = -35
    latmin = 25
    latmax = 50
    
elif zone == 'Acores':
    lonmin = -30
    lonmax = -24
    latmin = 37
    latmax = 40
    
elif zone == 'Med':
    lonmin = -10
    lonmax = 20
    latmin = 30
    latmax = 45

elif zone == 'COR': #coral sea
    lonmin = 140
    lonmax = 180
    latmin = -40
    latmax = -10
    
elif zone == 'tmp':
    lonmin = -35
    lonmax = -32
    latmin = 38.4
    latmax = 41

elif zone == 'tmp2':
    lonmin = -42
    lonmax = -25
    latmin = 32
    latmax = 42

elif zone == 'test':
    lonmin = -85
    lonmax = -60
    latmin = 25
    latmax = 45


# =============================================================================
# CODE
# =============================================================================
# Lecture du fichier d'entrée
nc_dico=ncl.read_nc(file_path,['traj_lat'])
variables = ['traj_lat','traj_lon','init_t', 'traj_time']
if hab_mode != 'void' and mortality:
    variables.append('SCL')
    variables.append('traj_temp')

       
# Read nc file
dico = ncl.read_nc(file_path, variables)

if jeanette:
    group = ncl.classify_lon_init(dico, nb_cat)
    print(group)
    if len(colors) != nb_cat:
        print('nb_cat has to be equal to len(colors)')
else:
    group = []

data_lists = ncl.data_lists(param, end_day, dico['init_t'])
pl.plot_animation_frames(gridfile, dico, hab_mode, To, lethargy, coef_SMR, start_day, end_day, h, [latmin, latmax],
                          [lonmin, lonmax], tracer, save_path, param, data_lists, last_turtle, mortality, group, nb_cat, colors, dpi)  
pl.convert_frames_to_video(save_path, videofile, fps)
