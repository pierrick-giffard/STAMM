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
from pathlib import Path, PurePath

#Personal librairies
sys.path.insert(1, os.path.join(sys.path[0], '../../LIB'))
import netCDF_lib as ncl
import plot_lib as pl
import IOlib as IO


# =============================================================================
# USERS PARAMETERS
# =============================================================================
#background
hab_mode = 'void' # 'food', 'temp', 'current', 'all', 'void'
#option 'void' for no background

#zone: 'Atlantic, 'Caribbean', 'Pacific', 'Indian'
zone = 'Atlantic'

# time delta between 2 frames (in days)
h = 10

#Starting age
start_age = 0

#Dates
start_day = 0
end_day = 50


# =============================================================================
# PATHS & FILES
# =============================================================================
file_path = sys.argv[1]
namelist = sys.argv[2]
param = IO.read_namelist(namelist, display=False)
directory = str(Path(file_path).parents[0])
#
gridfile = param['mesh_phy']
U_path = param['U_dir'] + '/'
V_path = param['V_dir'] + '/'
T_path = param['T_dir'] + '/'
food_path = param['food_dir'] + '/'
#Defaults save paths
save_path = directory + '/animation/'
filename = Path(namelist).stem.replace('namelist','')
videofile = save_path + '/' + filename + '_animation.avi'
#Rewrite them manually
# save_path =
# videofile =
if not Path(save_path).exists():
    save_path.mkdir(parents=True)

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
dpi = 350
##Zone de plot
#Pacifique

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

# Espacement des méridiens et des parallèles sur les figures
lon_space = 60
lat_space = 20
#
tracer = "PP"
mode = param['mode']
species = param['species']
# =============================================================================
# CODE
# =============================================================================
# Lecture du fichier d'entrée
nc_dico=ncl.read_nc(file_path,['traj_lat'])
nsteps,turtles=np.shape(nc_dico['traj_lat'])
variables = ['traj_lat','traj_lon','init_t', 'traj_time']
if hab_mode != 'void':
    variables.append('traj_temp')    

        
# Read nc file
dico = ncl.read_nc(file_path, variables)
pl.plot_animation_frames(gridfile, food_path, T_path, U_path, V_path, dico, hab_mode, To, lethargy,
                         coef_SMR, Fa, start_day, end_day, turtles, h, [latmin, latmax],
                         [lonmin, lonmax], lat_space, lon_space, tracer, mode, species,
                         start_age, save_path, filename, variables, dpi)  
pl.convert_frames_to_video(save_path, videofile, fps)
