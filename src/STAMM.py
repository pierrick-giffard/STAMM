#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Main code of STAMM, the Sea Turtle Active Movement Model.
Built with parcels version 2.2.0

Author: Pierrick Giffard working in Philippe Gaspar's team at Mercator Ocean.
February - May 2020.
"""

# =============================================================================
# IMPORTS
# =============================================================================
#Python libraries
from parcels import plotTrajectoriesFile,ErrorCode
from datetime import timedelta as delta
import sys, os
import time

#Personal libraries
import TurtleClass as tc
import Passive_kernels as pk
import Functions as fc
sys.path.insert(1, os.path.join(sys.path[0], '../LIB'))
import IOlib as IO

#Initial time
t0=time.time()
# =============================================================================
# PARAMETERS
# =============================================================================
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
#Temporal loop over data for time_periodic days
time_periodic = param['time_periodic'] 

# =============================================================================
# Read initial positions and time
# =============================================================================
lon_init, lat_init, t_init = IO.read_positions(param)
ndays_simu += max(t_init) - min(t_init)


# =============================================================================
# FIELDSET, CLASS AND PARTICLESET
# ========================================================================
fieldset = fc.create_fieldset(param, ndays_simu, t_init)
print((fieldset.U.__dict__['lat'][-1] - fieldset.U.__dict__['lat'][0])/len(fieldset.U.__dict__['lat']))
fc.initialization(fieldset, ndays_simu, param)
t_release = fc.compute_t_release(t_init, fieldset, param)

turtle = tc.define_Turtle_Class(fieldset, param)
pset = fc.create_particleset(fieldset, turtle, lon_init, lat_init, t_release, param) 

# =============================================================================
# KERNELS
# =============================================================================
k_adv = fc.define_advection_kernel(pset, param)
k_active = fc.define_active_kernels(pset, param)
k_passive = fc.define_passive_kernels(fieldset, pset, param) 
#
kernels = fc.sum_kernels(k_adv, k_active, k_passive)

    
# =============================================================================
# COMPUTATION
# =============================================================================
output_file = pset.ParticleFile(name=OutputFile, outputdt=delta(seconds=t_output))
pset.execute(kernels, runtime=delta(days=ndays_simu), dt=delta(seconds=tstep),\
              output_file=output_file,\
              recovery={ErrorCode.ErrorOutOfBounds: pk.DisableParticle})


# =============================================================================
# OUTPUT
# =============================================================================
output_file.export()
#plotTrajectoriesFile(OutputFile)
fc.modify_output(OutputFile, t_init, param)


tt=time.time()-t0
print('\n')
print('Total execution time: '+ str(delta(seconds=int(tt))))
print('\n')
