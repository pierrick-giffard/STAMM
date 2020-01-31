#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 08:53:40 2020

@author: pgiffard
"""

def SampleP(particle, fieldset, time): 
    particle.T = fieldset.T[time, particle.depth, particle.lat, particle.lon]
    particle.NPP = fieldset.NPP[time, particle.depth, particle.lat, particle.lon]
            


def periodicBC(particle, fieldset, time):
    if particle.lon < fieldset.halo_west:
        particle.lon += fieldset.halo_east - fieldset.halo_west
    elif particle.lon > fieldset.halo_east:
        particle.lon -= fieldset.halo_east - fieldset.halo_west




def define_additional_kernels(fieldset, pset, key_alltracers, key_periodic, grid_type):
    """
    Parameters:
        -fieldset
        -pset
        -key_alltracers: if True, T and NPP are sampled
        -key_periodic: True for east/west periodicity
        -grid_type: 'orca' or 'standard'
    """
    kernels_list = []
    if key_alltracers:
        kernels_list.append(SampleP)
    
    if key_periodic:
        kernels_list.append(periodicBC)
        if grid_type == 'orca':
            #give halos east and west
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
        
    
    
    
    
    
    
    
    
    
    
    
    
    