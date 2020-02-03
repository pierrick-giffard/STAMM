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
from parcels import ParticleSet, plotTrajectoriesFile
from datetime import timedelta as delta
import time
import sys

#Personal libraries
import IOlib as IO
import TurtleClass as tc
import Advection_kernel as adv
import Additional_kernels as add
import Turtle_kernels as tk
import Functions as fc

# =============================================================================
# PARAMETERS
# =============================================================================
#Initial time
t0=time.time()

#Users arguments
namelist = sys.argv[1]
OutputFile = sys.argv[2]

#Read namelist
param = IO.read_namelist(namelist)
IO.check_param(param,OutputFile)

#Time step in seconds
tstep = param['tstep']
#Number of time steps
nsteps_simu = param['nsteps_simu']
#Output frequency
t_output = param['t_output']


 
# =============================================================================
# Read initial positions and time
# =============================================================================
lon_init, lat_init, t_init = IO.read_positions(param)


# =============================================================================
# FIELDSET, CLASS AND PARTICLESET
# =============================================================================
fieldset = fc.build_fieldset(param)
turtle = tc.define_Turtle_Class(fieldset,param)
pset = ParticleSet(fieldset, pclass=turtle, lon=lon_init, lat=lat_init,time=t_init)
#
fc.initialization(pset, param)

# =============================================================================
# KERNELS
# =============================================================================
k_adv = adv.define_advection_kernel(pset, param)
k_turtle = tk.define_turtle_kernels(pset, param)
k_add = add.define_additional_kernels(fieldset, pset, param) 
#
kernels = add.sum_kernels(k_adv, k_turtle, k_add)
   
    
# =============================================================================
# COMPUTATION
# =============================================================================
output_file = pset.ParticleFile(name=OutputFile, outputdt=delta(seconds=t_output))
pset.execute(kernels, runtime=delta(seconds=nsteps_simu*tstep), dt=delta(seconds=tstep),output_file=output_file)


tt=time.time()-t0
print('\n')
print('Total execution time: '+ str(delta(seconds=int(tt))))
print('\n')

################## PLOT #####################################
plotTrajectoriesFile(OutputFile)
