#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Animate STAMM output files.
Execute as:
python plot_animation.py output.nc namelist
"""
# =============================================================================
# IMPORTS
# =============================================================================
import sys, os
import numpy as np
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
hab_mode = 'food' # 'food', 'temp', 'current', 'all', 'void'
#option 'void' for no background
mortality = False #set to False not to calculate dead turtles

#zone: 'Atlantic, 'Caribbean', 'Pacific', 'Indian', 'Gulf_stream', 'Acores', 'Med'
zone = 'Acores'

# time delta between 2 frames (in days)
h = 100


#Dates
start_day = 0
end_day = 5000


# =============================================================================
# PATHS & FILES
# =============================================================================
file_path = sys.argv[1]
namelist = sys.argv[2]
param = IO.read_namelist(namelist, display=False)
directory = str(Path(file_path).parents[0])
#
gridfile = '/data/rd_exchange2/tcandela/STAMM/ressources/VGPM/VGPM_083_mesh.nc' #gridfile for NPP lon/lat
food_path = param['food_dir'] + '/'
#Defaults save paths
save_path = directory + '/animation/'
filename = Path(namelist).stem.replace('namelist','')
videofile = save_path + '/' + filename + '_animation.avi'
if not Path(save_path).exists():
    Path(save_path).mkdir(parents=True)

print('\n')
print('********************************************************************************')
print("Writing frames in ", save_path)
print('********************************************************************************')
print('\n')



# =============================================================================
# AUTO PARAMETERS
# =============================================================================
#Température optimale pour le calcul d'habitat
To = 24.

# Temps de léthargie (nombre de jours maximum dans Tw<Tmin)
lethargy=30.

# Variables pour le calcul de l'habitat thermique et alimentaire
coef_SMR=5.
Fa = param['P0']

# Nombre d'images par seconde pour la vidéo
fps = 8 #8
dpi = 150


# =============================================================================
#ZONES
# =============================================================================
if zone == 'Atlantic':
    lonmin = -110.
    lonmax = 36.
    latmin = -7
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
    lonmin = -80
    lonmax = -50
    latmin = 30
    latmax = 45
    
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


#
tracer = "PP" #PP or mnk
species = param['species']
# =============================================================================
# CODE
# =============================================================================
# Lecture du fichier d'entrée
nc_dico=ncl.read_nc(file_path,['traj_lat'])
nsteps,nturtles=np.shape(nc_dico['traj_lat'])
variables = ['traj_lat','traj_lon','init_t', 'traj_time']
if hab_mode != 'void' and mortality:
    variables.append('traj_temp')

       
# Read nc file
dico = ncl.read_nc(file_path, variables)
data_lists = ncl.data_lists(param, end_day, dico['init_t'])
pl.plot_animation_frames(gridfile, food_path, dico, hab_mode, To, lethargy,
                          coef_SMR, Fa, start_day, end_day, nturtles, h, [latmin, latmax],
                          [lonmin, lonmax], tracer, species, save_path, param, data_lists, mortality, dpi)  
pl.convert_frames_to_video(save_path, videofile, fps)
