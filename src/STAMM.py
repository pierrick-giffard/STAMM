#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Main code of STAMM.
Need parcels version 2.0.0.
The function AdvectionRK4 has to be modified.

Authors: Pierrick Giffard, Philippe Gaspar at Mercator Ocean.
February 2020.
"""

# =============================================================================
# IMPORTS
# =============================================================================
#Python libraries
from parcels import FieldSet, ParticleSet, JITParticle, ParticleFile, plotTrajectoriesFile, Variable,ErrorCode,Field
import numpy as np
from datetime import timedelta as delta
from glob import glob
import time
import sys
import netCDF4 as nc
import matplotlib.pyplot as plt

#Personal libraries
import IOlib as IO
import TurtleClass as tc
import Advection_kernel as adv
import Additional_kernels as add
import Turtle_kernels as tk
import Functions as fc

# =============================================================================
# INITIALIZATION
# =============================================================================
#Initial time
t0=time.time()

#Read users arguments
namelist = sys.argv[1]
OutputFile = sys.argv[2]

#Read namelist
param = IO.read_namelist(namelist)
IO.check_param(param,OutputFile)

# =============================================================================
# PARAMETERS
# =============================================================================
#Time step in seconds
tstep = param['tstep']
#Number of steps
nsteps_simu = param['nsteps_simu']
#Number of turtles
nturtles = param['nturtles']
#Mode can be 'passive' (no intended movement, the IBM works exactly as Ariane in this mode), 'diffusion' (passive + random difffusion) or 'active' if directed movement is to be modelled
mode = param['mode']
#alpha is a parameter representing the easiness that turtles have to follow the habitat gradient, the lower alpha is, the higher the diffusion will be
alpha = param['alpha']
#Boolean, True if all tracers are to be used
key_alltracers = param['key_alltracers']
#Species
species = param['species']

# =============================================================================
# DATA tmp !!!!!!!!!!!!!
# =============================================================================
data_path = '/homelocal/pgiffard/test_PARCELS/PSY_025degORCA/test_anna/'
mesh_mask = data_path+'mesh_hgr_PSY4V3_deg.nc'
filenames = {'U': {'lon': mesh_mask, 'lat': mesh_mask, 'data': data_path+'ORCA025deg_PSY4V3R1_1dAV_gridU_y2018m12d29.nc'},
            'V': {'lon': mesh_mask, 'lat': mesh_mask, 'data': data_path+'ORCA025deg_PSY4V3R1_1dAV_gridV_y2018m12d29.nc'},
            'T': {'lon': mesh_mask, 'lat': mesh_mask, 'data': data_path+'ORCA025deg_PSY4V3R1_1dAV_gridT_y2018m12d29.nc'},
            'NPP': {'lon': mesh_mask, 'lat': mesh_mask, 'data': data_path+'mercatorpsy4v3_deg_GLO_20181225_2D.nc'}}
#add NPP
variables = {'U': 'vozocrtx','V': 'vomecrty','T': 'votemper','NPP': 'npp'}
dimensions = {'lon': 'glamf', 'lat': 'gphif', 'time': 'time_counter'}   #need f nodes 


# =============================================================================
# FIELDSET
# =============================================================================
fieldset = FieldSet.from_netcdf(filenames, variables, dimensions, allow_time_extrapolation=True, deferred_load=True)  #time_periodic=delta(days=365)

 
# =============================================================================
# Read initial positions and time
# =============================================================================
lon_init, lat_init, t_init = IO.read_positions(param)


# =============================================================================
# CLASS AND PARTICLESET
# =============================================================================
turtle = tc.define_Turtle_Class(fieldset)
pset = ParticleSet(fieldset, pclass=turtle, lon=lon_init, lat=lat_init,time=t_init)

for p in pset:  
    p.vscale = param['vscale']
    p.P0 = param['P0']
    p.dx = 10000 #tmp, dx for gradient calculation
    p.alpha = alpha
    if mode == 'active':
        p.mode = 1
    elif mode == 'passive':
        p.mode = 0
    p.tstep = tstep

# =============================================================================
# KERNELS
# =============================================================================
adv_scheme = 'RK4' #tmp
grid_type = 'standard'#tmp
#
k_adv = adv.define_advection_kernel(pset, mode, adv_scheme)
k_turtle = tk.define_turtle_kernels(pset, mode, species)
k_add = add.define_additional_kernels(fieldset, pset, key_alltracers, param['key_periodic'], grid_type) 
#
kernels = add.sum_kernels(k_adv, k_turtle, k_add)
   
    
##################### OUTPUT FILE #######################################
t_output=24 #tmp




# =============================================================================
# COMPUTATION
# =============================================================================
output_file = pset.ParticleFile(name=OutputFile, outputdt=delta(hours=t_output))
pset.execute(kernels, runtime=delta(seconds=nsteps_simu*tstep), dt=delta(seconds=tstep),output_file=output_file)



################## PLOT #####################################
plotTrajectoriesFile(OutputFile)
# test = nc.Dataset('/homelocal/pgiffard/SRC/STAMM/EXP/PARCELS/namelist.nc')
# traj_lat = np.squeeze(test.variables['lat'])
# print(traj_lat.shape)
# plt.figure()
# for k in range(nturtles):
#     plt.plot(traj_lat)
# plt.show()
