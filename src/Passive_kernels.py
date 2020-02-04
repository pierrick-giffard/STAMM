#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 08:53:40 2020

@author: pgiffard
"""
import math

def IncrementAge(particle, fieldset, time):
   "Age in days. Needs to be the last kernel to actualize the age at the end of the timestep"
   particle.age += particle.dt/86400.


def SampleTracers(particle, fieldset, time): 
    particle.T = fieldset.T[time, particle.depth, particle.lat, particle.lon]
    particle.NPP = fieldset.NPP[time, particle.depth, particle.lat, particle.lon]
            


def Periodic(particle, fieldset, time):
    if particle.lon < fieldset.halo_west:
        particle.lon += fieldset.halo_east - fieldset.halo_west
    elif particle.lon > fieldset.halo_east:
        particle.lon -= fieldset.halo_east - fieldset.halo_west



def Distance(particle, fieldset, time):
    # Calculate the distance in latitudinal direction (using 1.11e2 kilometer per degree latitude)
    particle.lat_dist = (particle.lat - particle.prev_lat) * particle.deg
    # Calculate the distance in longitudinal direction, using cosine(latitude) - spherical earth
    particle.lon_dist = (particle.lon - particle.prev_lon) * particle.deg * math.cos(particle.lat * math.pi / 180)
    # Calculate the total Euclidean distance travelled by the particle
    particle.distance = math.sqrt(math.pow(particle.lon_dist, 2) + math.pow(particle.lat_dist, 2))

    particle.prev_lon = particle.lon
    particle.prev_lat = particle.lat



def CurrentVelocity(particle, fieldset, time):
    """
    Compute current mean velocity during a tstep.
    In active mode, correct if V swim is constant over the time step.
    Compute first TotalDistance, then CurrentVelocity.
    """
    if particle.mode == 1: #active
        particle.u_current = particle.lon_dist / particle.tstep - particle.u_swim
        particle.v_current = particle.lat_dist / particle.tstep - particle.v_swim
    
    elif particle.mode == 0: #passive
        particle.u_current = particle.lon_dist / particle.tstep
        particle.v_current = particle.lat_dist / particle.tstep


def define_passive_kernels(fieldset, pset, param):
    """
    Parameters:
        -fieldset
        -pset
        -param: needs key_alltracers (if True, T and NPP are sampled),
        periodicBC (True for east/west periodicity) and grid_type ('orca' or 'standard')
    """
    key_alltracers = param['key_alltracers']
    periodicBC = param['periodicBC']
    grid_type = param['grid_type']
    #
    kernels_list = [Distance, CurrentVelocity]
    if key_alltracers:
        kernels_list.append(SampleTracers)
    
    if periodicBC:
        kernels_list.append(Periodic)
        if grid_type == 'orca':
            #give halos east and west: how to do it?????
            a=1#give the first and the last column?
        elif grid_type == 'standard':
            try:
                fieldset.add_constant('halo_west', fieldset.U.grid.lon[0,0])
                fieldset.add_constant('halo_east', fieldset.U.grid.lon[0,-1])
            except:
                fieldset.add_constant('halo_west', fieldset.U.grid.lon[0])
                fieldset.add_constant('halo_east', fieldset.U.grid.lon[-1])
            #
            fieldset.add_periodic_halo(zonal=True)
    #
    kernels_list.append(IncrementAge)
    #       
    for k in range(len(kernels_list)):
        kernels_list[k]=pset.Kernel(kernels_list[k])  
    return kernels_list    



def sum_kernels(k_adv, k_active, k_passive):
    """
    Sums all the kernels and returns the summed kernel.
    The position is important.
    """
    if k_active != []:
        kernels = k_active[0]
        print_kernels = [k_active[0].funcname]
        #
        for k in k_active[1:]:
            kernels = kernels + k
            print_kernels.append(k.funcname)
        kernels += k_adv
        print_kernels.append(k_adv.funcname)
    #    
    else:   
        kernels = k_adv
        print_kernels = [k_adv.funcname]
    #
    for k in k_passive:
        kernels = kernels + k
        print_kernels.append(k.funcname)
    #
    print('****************************************************')
    print("These kernels will be used for computation: \n")
    for k in print_kernels:
        print(k, "\n")
    print('****************************************************')
    
    return kernels
        
    
    
    
    
    
    
    
    
    
    
    
    
    