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
import numpy as np

#Personal librairies
sys.path.insert(1, os.path.join(sys.path[0], '../../LIB'))
import netCDF_lib as ncl
import plot_lib as pl
import IOlib as IO


# =============================================================================
# AUTO PARAMETERS
# =============================================================================
file_path = sys.argv[1]
namelist = sys.argv[2]
param = IO.read_namelist(namelist, display=False)

last_turtle = param['nturtles'] # Plot first n turtles
start_day = 0
end_day = param['ndays_simu']

# =============================================================================
# USER PARAMETERS
# =============================================================================

#background
hab_mode = 'current' # 'food', 'temp', 'tot', 'current', 'void'
#option 'void' for no background,

mortality = False #set to False not to calculate dead turtles

# Jeanette project: different colors depending on release lon
jeanette = True
nb_cat = 4
colors = ['darkviolet','blue','green','orange','red']

# plot zone
zone = 'jeanette_1month'

# time delta between 2 frames (in days)
h = 1

## Overwrite Auto parameters
#last_turtle = 100
#start_day = 50
#end_day = 365

#gridfile for NPP lon/lat
gridfile = '/data/rd_exchange2/tcandela/STAMM/ressources/VGPM/VGPM_083_mesh.nc'

# Video
fps = 2 #images/sec
dpi = 350 #images resolution

# =============================================================================
# PATHS & FILES
# =============================================================================

#Defaults save paths
directory = ncl.get_directory(file_path)
filename = ncl.get_name(file_path)
save_path = directory + 'animation_' + filename + '/'
if not Path(save_path).exists():
    Path(save_path).mkdir(parents=True)
#
videofile = save_path + '/' + filename + '_animation.avi'
#   

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

elif zone == 'jeanette_1year':
    lonmin = -85
    lonmax = -10
    latmin = 20
    latmax = 50

elif zone == 'jeanette_1month':
    lonmin = -81
    lonmax = -78#-59
    latmin = 26#24
    latmax = 27#42


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
    if len(colors) != nb_cat:
        print('nb_cat has to be equal to len(colors)')
else:
    group = []

data_lists = ncl.data_lists(param, end_day, np.float64(dico['init_t']))
pl.plot_animation_frames(gridfile, dico, hab_mode, To, lethargy, coef_SMR, start_day, end_day, h, [latmin, latmax],
                          [lonmin, lonmax], save_path, param, data_lists, last_turtle, mortality, group, nb_cat, colors, dpi)  
pl.convert_frames_to_video(save_path, videofile, fps)
