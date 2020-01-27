#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 15:58:58 2020

@author: pgiffard
"""

# =============================================================================
# IMPORTS
# =============================================================================
from parcels import FieldSet, ParticleSet, JITParticle, ParticleFile, plotTrajectoriesFile, Variable,ErrorCode,Field
from parcels import AdvectionRK4
import numpy as np
from datetime import timedelta as delta
from glob import glob
import matplotlib.pyplot as plt
from netCDF4 import Dataset


########################### PARAMETERS ##############################
nb_turtles = 4
t_simu = 200 #simulation duration (days)
DT = 24 #time step (hours)
t_output = 24 #writing time step (hours)
output_name = 'orca_passive.nc'


##################### DATA #######################################
data_path = '/homelocal/pgiffard/test_PARCELS/PSY_025degORCA/test_anna/'
mesh_mask = data_path+'mesh_hgr_PSY4V3_deg.nc'
filenames = {'U': {'lon': mesh_mask, 'lat': mesh_mask, 'data': data_path+'ORCA025deg_PSY4V3R1_1dAV_gridU_y2018m12d29.nc'},
            'V': {'lon': mesh_mask, 'lat': mesh_mask, 'data': data_path+'ORCA025deg_PSY4V3R1_1dAV_gridV_y2018m12d29.nc'},
            'T': {'lon': mesh_mask, 'lat': mesh_mask, 'data': data_path+'ORCA025deg_PSY4V3R1_1dAV_gridT_y2018m12d29.nc'}}

variables = {'U': 'vozocrtx','V': 'vomecrty','T': 'votemper'}
dimensions = {'lon': 'glamf', 'lat': 'gphif', 'time': 'time_counter'}   #need f nodes 


##################### FIELDSET #######################################
fieldset = FieldSet.from_netcdf(filenames, variables, dimensions, allow_time_extrapolation=True, deferred_load=True)  #time_periodic=delta(days=365)


################### RELEASE DATA #################################
init = open(data_path+'../YalimapoCayenne_summer_sample.txt','r')
x_init, y_init, t_init = np.loadtxt(init,usecols=(0,1,3),unpack=True)
x_init = (36-(1440-x_init)/4)[0:nb_turtles]
y_init = (y_init/4-60)[0:nb_turtles]


##################### PARTICLE SET #######################################
class turtle(JITParticle): #Leatherback?
    T = Variable('T', initial=fieldset.T)
    u_swim = Variable('u_swim', to_write=True, dtype=np.float32)
    v_swim = Variable('v_swim', to_write=True, dtype=np.float32)
    
pset = ParticleSet(fieldset, pclass=turtle, lon=x_init, lat=y_init,time=0)


################## ADDITIONAL KERNELS ###############################
def SampleP(particle, fieldset, time): 
    particle.T = fieldset.T[time, particle.depth, particle.lat, particle.lon]

k_sample = pset.Kernel(SampleP)
  

def add_swimming_velocity(particle,fieldset,time):
    particle.u_swim, particle.v_swim = compute_swimming_velocity(particle,fieldset,time)
    
    
k_swim = pset.Kernel(add_swimming_velocity)

kernels = pset.Kernel(AdvectionRK4) + k_swim + k_sample

##################### OUTPUT FILE #######################################
output_file = pset.ParticleFile(name=output_name, outputdt=delta(hours=t_output))


################## COMPUTATION #####################################
pset.execute(kernels, runtime=delta(days=t_simu), dt=delta(hours=DT),output_file=output_file)

################## PLOT #####################################
plotTrajectoriesFile(output_name)