#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Definition of kernels used for all particles (both in active and passive modes).
"""

import math


def IncrementAge(particle, fieldset, time):
   """"
   Increment turtles age (in days).
   """
   particle.age += particle.dt/86400.
   #
   if particle.age > fieldset.ndays_simu + 1:
       particle.active = 0



def SampleTracers(particle, fieldset, time):
    """
    Sample tracers at particle location in passive mode.
    In active mode sampling is integrated to function compute_habitat.
    """
    if particle.active == 1:
        particle.T = fieldset.T[time, particle.depth, particle.lat, particle.lon]
        particle.NPP = fieldset.NPP[time, particle.depth, particle.lat, particle.lon]
            


def Periodic(particle, fieldset, time):
    """
    So that particles move from east boundary to west boundary.
    """
    if particle.active == 1:
        if particle.lon < fieldset.halo_west:
            particle.lon += fieldset.halo_east - fieldset.halo_west
        elif particle.lon > fieldset.halo_east:
            particle.lon -= fieldset.halo_east - fieldset.halo_west



    
    
def DeleteParticle(particle, fieldset, time):
    """ 
    Delete out of bounds particles.
    """
    print("\n")
    print("Turtle [%d] deleted at lon,lat = %f,%f and time = %f"%(particle.id,particle.lon,particle.lat,time))
    particle.delete()
    
    
def DisableParticle(particle, fieldset, time):
    """ 
    Delete out of bounds particles.
    """
    print("\n")
    print("Turtle [%d] disabled at lon,lat = %f,%f and time = %f"%(particle.id,particle.lon,particle.lat,time))
    particle.active = 0


def BeachTesting(particle, fieldset, time):
    """
    If particle is on land, particle.beached is set to 1.
    """
    if particle.active == 1:
        (u, v) = fieldset.UV[time, particle.depth, particle.lat, particle.lon]
        if math.fabs(u) < 1e-14 and math.fabs(v) < 1e-14:
            particle.beached = 1




           
def UndoMove(particle, fieldset, time):
    """
    Send particle back to last position in case it is on land.
    If it is on land more than onland_max times in a row, it is deleted.
    """
    onland_max = 50  
    if particle.beached == 1:
        particle.beached = 0
        particle.lon = particle.prev_lon
        particle.lat = particle.prev_lat
        particle.onland += 1
        #
        if particle.onland > onland_max:
            print("Particle [%d] was disabled after beaching 50 times in a row at lon,lat = %f,%f"%(particle.id,particle.lon,particle.lat))
            particle.active = 0
    else:
        particle.onland = 0



def CheckOnLand(particle,fieldset,time):
    """ 
    Check if particles are released on land. 
    Executed once with dt=0 just after ParticleSet is created.
    """
    (u, v) = fieldset.UV[0, particle.depth, particle.lat, particle.lon]
    if math.fabs(u) < 1e-14 and math.fabs(v) < 1e-14:
        print("Particle [%d] is released on land at lon,lat = %f,%f. Execution stops."%(particle.id,particle.lon,particle.lat))
        exit(0)
        
    
    
    
    
    
    
    
    
    
    
    
    
    