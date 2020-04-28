#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Main code of STAMM, the Sea Turtle Active Movement Model.
Built with parcels version 2.1.4.

Author: Pierrick Giffard working in Philippe Gaspar's team at Mercator Ocean.
February 2020.
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
fc.PSY_patch(fieldset,param)
fc.initialization(fieldset, ndays_simu, param)
turtle = tc.define_Turtle_Class(fieldset,param)
pset = fc.create_particleset(fieldset, turtle, lon_init, lat_init, t_init, param) 

"""
print("\n")
print("\n")
print(fieldset.__dict__)
print("\n")
print("\n")
print(fieldset.U.__dict__)
print("\n")
print("\n")
print(fieldset.npp.__dict__)
print("\n")
print("\n")
print(fieldset.U.grid.__dict__)
print("\n")
print("\n")
print(fieldset.npp.grid.__dict__)
"""



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
#plotTrajectoriesFile('pb_interp.nc',tracerfile='/data/rd_exchange2/pgiffard/DATA/GLORYS12/GLORYS12_PGS_2D_20020606_UVT.nc',tracerlon='longitude',tracerlat='latitude',tracerfield='uo');
#plotTrajectoriesFile(OutputFile)
fc.modify_output(OutputFile, t_init, param)


tt=time.time()-t0
print('\n')
print('Total execution time: '+ str(delta(seconds=int(tt))))
print('\n')
