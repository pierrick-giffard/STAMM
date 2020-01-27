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
from parcels import AdvectionRK4
import numpy as np
from datetime import timedelta as delta
from glob import glob
import time
import sys

#Personal libraries
import IOlib as IO
import Kernels as ke
import TurtleClass as tc

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

# =============================================================================
# DATA tmp !!!!!!!!!!!!!
# =============================================================================
data_path = '/homelocal/pgiffard/test_PARCELS/PSY_025degORCA/test_anna/'
mesh_mask = data_path+'mesh_hgr_PSY4V3_deg.nc'
filenames = {'U': {'lon': mesh_mask, 'lat': mesh_mask, 'data': data_path+'ORCA025deg_PSY4V3R1_1dAV_gridU_y2018m12d29.nc'},
            'V': {'lon': mesh_mask, 'lat': mesh_mask, 'data': data_path+'ORCA025deg_PSY4V3R1_1dAV_gridV_y2018m12d29.nc'},
            'T': {'lon': mesh_mask, 'lat': mesh_mask, 'data': data_path+'ORCA025deg_PSY4V3R1_1dAV_gridT_y2018m12d29.nc'}}

variables = {'U': 'vozocrtx','V': 'vomecrty','T': 'votemper'}
dimensions = {'lon': 'glamf', 'lat': 'gphif', 'time': 'time_counter'}   #need f nodes 


##################### FIELDSET #######################################
fieldset = FieldSet.from_netcdf(filenames, variables, dimensions, allow_time_extrapolation=True, deferred_load=True)  #time_periodic=delta(days=365)

 
# =============================================================================
# Read initial positions and time
# =============================================================================
lon_init, lat_init, t_init = IO.read_positions(param)


##################### PARTICLE SET #######################################
turtle = tc.define_Turtle_Class(fieldset)
pset = ParticleSet(fieldset, pclass=turtle, lon=lon_init, lat=lat_init,time=0)


# =============================================================================
# 
# =============================================================================
def SampleP(particle, fieldset, time): 
    particle.T = fieldset.T[time, particle.depth, particle.lat, particle.lon]

k_sample = pset.Kernel(SampleP)
  

def add_swimming_velocity(particle,fieldset,time):
    particle.u_swim, particle.v_swim = 0., 0.
    #print(particle.age)
    
    
k_swim = pset.Kernel(add_swimming_velocity)

kernels = pset.Kernel(AdvectionRK4) + k_swim + k_sample

##################### OUTPUT FILE #######################################
t_output=24 #tmp
output_file = pset.ParticleFile(name=OutputFile, outputdt=delta(hours=t_output))


################## COMPUTATION #####################################
pset.execute(kernels, runtime=delta(seconds=nsteps_simu*tstep), dt=delta(seconds=tstep),output_file=output_file)

################## PLOT #####################################
plotTrajectoriesFile(OutputFile)