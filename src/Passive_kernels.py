#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
All passive kernels (also used in active mode).
"""

import math

def IncrementAge(particle, fieldset, time):
   "Increment turtles age (in days)."
   particle.age += particle.dt/86400.


def SampleTracers(particle, fieldset, time):
    """
    Sample tracers at particle location.
    """
    particle.T = fieldset.T[time, particle.depth, particle.lat, particle.lon]
    particle.NPP = fieldset.NPP[time, particle.depth, particle.lat, particle.lon]
            


def Periodic(particle, fieldset, time):
    """
    Particles pass from east boundary to west boundary.
    """
    if particle.lon < fieldset.halo_west:
        particle.lon += fieldset.halo_east - fieldset.halo_west
    elif particle.lon > fieldset.halo_east:
        particle.lon -= fieldset.halo_east - fieldset.halo_west



def Distance(particle, fieldset, time):
    """
    Calculate the distance travelled at each time step.
    """
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
    In active mode, this calculation is correct if V swim is constant over the whole time step.
    """
    if particle.mode == 1: #active
        particle.u_current = particle.lon_dist / particle.tstep - particle.u_swim
        particle.v_current = particle.lat_dist / particle.tstep - particle.v_swim
    
    elif particle.mode == 0: #passive
        particle.u_current = particle.lon_dist / particle.tstep
        particle.v_current = particle.lat_dist / particle.tstep


def DeleteParticle(particle, fieldset, time):
    """ 
    Delete out of Bounds particles
    """
    print("\n")
    print("Turtle [%d] deleted at lon,lat = %f,%f and time = %f"%(particle.id,particle.lon,particle.lat,time))
    particle.delete()



def UndoMove(particle, fieldset, time):
    """
    Send turtle back to last position in case it is on land.
    """
    (u, v) = fieldset.UV[time, particle.depth, particle.lat, particle.lon]
    if fabs(u) < 1e-14 and fabs(v) < 1e-14: #does not work with abs
        print("Particle [%d] is on land at lon,lat = %f,%f"%(particle.id,particle.lon,particle.lat))
        print(" => It is sent back to its last position")
        particle.lon = particle.prev_lon
        particle.lat = particle.prev_lat

        
    

        
    
    
    
    
    
    
    
    
    
    
    
    
    