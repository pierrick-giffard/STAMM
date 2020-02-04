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
import Passive_kernels as pk
import Active_kernels as ak
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
ndays_simu = param['ndays_simu']
#Output frequency
t_output = param['t_output']

 
# =============================================================================
# Read initial positions and time
# =============================================================================
lon_init, lat_init, t_init = IO.read_positions(param)
t_init *= 86400 #from days to seconds

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
k_active = ak.define_turtle_kernels(pset, param)
k_passive = pk.define_passive_kernels(fieldset, pset, param) 
#
kernels = pk.sum_kernels(k_adv, k_active, k_passive)

    
# =============================================================================
# COMPUTATION
# =============================================================================
output_file = pset.ParticleFile(name=OutputFile, outputdt=delta(seconds=t_output))
pset.execute(kernels, runtime=delta(days=ndays_simu), dt=delta(seconds=tstep),output_file=output_file)


tt=time.time()-t0
print('\n')
print('Total execution time: '+ str(delta(seconds=int(tt))))
print('\n')

################## PLOT #####################################
plotTrajectoriesFile(OutputFile)
