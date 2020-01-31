#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 08:53:40 2020

@author: pgiffard
"""
import math

def SampleTracers(particle, fieldset, time): 
    particle.T = fieldset.T[time, particle.depth, particle.lat, particle.lon]
    particle.NPP = fieldset.NPP[time, particle.depth, particle.lat, particle.lon]
            


def periodicBC(particle, fieldset, time):
    if particle.lon < fieldset.halo_west:
        particle.lon += fieldset.halo_east - fieldset.halo_west
    elif particle.lon > fieldset.halo_east:
        particle.lon -= fieldset.halo_east - fieldset.halo_west



def TotalDistance(particle, fieldset, time):
    # Calculate the distance in latitudinal direction (using 1.11e2 kilometer per degree latitude)
    particle.lat_dist = (particle.lat - particle.prev_lat) * particle.deg
    # Calculate the distance in longitudinal direction, using cosine(latitude) - spherical earth
    particle.lon_dist = (particle.lon - particle.prev_lon) * particle.deg * math.cos(particle.lat * math.pi / 180)
    # Calculate the total Euclidean distance travelled by the particle
    particle.distance += math.sqrt(math.pow(particle.lon_dist, 2) + math.pow(particle.lat_dist, 2))

    particle.prev_lon = particle.lon  # Set the stored values for next iteration.
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


def define_additional_kernels(fieldset, pset, key_alltracers, key_periodic, grid_type):
    """
    Parameters:
        -fieldset
        -pset
        -key_alltracers: if True, T and NPP are sampled
        -key_periodic: True for east/west periodicity
        -grid_type: 'orca' or 'standard'
    """
    kernels_list = [TotalDistance, CurrentVelocity]
    if key_alltracers:
        kernels_list.append(SampleTracers)
    
    if key_periodic:
        kernels_list.append(periodicBC)
        if grid_type == 'orca':
            #give halos east and west: how to do it?????
            a=1
        elif grid_type == 'standard':
            try:
                fieldset.add_constant('halo_west', fieldset.U.grid.lon[0,0])
                fieldset.add_constant('halo_east', fieldset.U.grid.lon[0,-1])
            except:
                fieldset.add_constant('halo_west', fieldset.U.grid.lon[0])
                fieldset.add_constant('halo_east', fieldset.U.grid.lon[-1])
            
            fieldset.add_periodic_halo(zonal=True)
    for k in kernels_list:
        k=pset.Kernel(k)
    return kernels_list
        



def sum_kernels(k_adv, k_turtle, k_add):
    """
    Sums all the kernels and returns the summed kernel.
    """
    kernels = k_adv
    for k in k_add:
        kernels = kernels + k
    for k in k_turtle: #TurtleAge has to be the last kernel
        kernels = kernels + k
    
    return kernels
        
    
    
    
    
    
    
    
    
    
    
    
    
    