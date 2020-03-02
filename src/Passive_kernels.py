#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Definition of kernels used for all particles (both in active and passive modes).
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
    So that particles move from east boundary to west boundary.
    """
    if particle.lon < fieldset.halo_west:
        particle.lon += fieldset.halo_east - fieldset.halo_west
    elif particle.lon > fieldset.halo_east:
        particle.lon -= fieldset.halo_east - fieldset.halo_west



def Distance(particle, fieldset, time):
    """
    Calculate the distance travelled at each time step.
    Flat earth assumption during a timestep.
    """
    # Calculate the distance in latitudinal direction (using 1.11e2 kilometer per degree latitude)
    particle.lat_dist = (particle.lat - particle.prev_lat) * fieldset.deg
    # Calculate the distance in longitudinal direction, using cosine(latitude) - spherical earth
    particle.lon_dist = (particle.lon - particle.prev_lon) * fieldset.deg * math.cos(particle.lat * math.pi / 180)
    # Calculate the total Euclidean distance travelled by the particle
    particle.distance = math.sqrt(math.pow(particle.lon_dist, 2) + math.pow(particle.lat_dist, 2))

    particle.prev_lon = particle.lon
    particle.prev_lat = particle.lat



def CurrentVelocity(particle, fieldset, time):
    """
    Compute current mean velocity during a tstep.
    This calculation is correct if the swimming velocity is constant over the whole time step.
    """
    if fieldset.active == 1:
        particle.u_current = particle.lon_dist / fieldset.tstep - particle.u_swim
        particle.v_current = particle.lat_dist / fieldset.tstep - particle.v_swim
    else:
        particle.u_current = particle.lon_dist / fieldset.tstep
        particle.v_current = particle.lat_dist / fieldset.tstep

    
    
def DeleteParticle(particle, fieldset, time):
    """ 
    Delete out of bounds particles.
    """
    print("\n")
    print("Turtle [%d] deleted at lon,lat = %f,%f and time = %f"%(particle.id,particle.lon,particle.lat,time))
    particle.delete()



def BeachTesting(particle, fieldset, time):
    """
    If particle is on land, particle.beached is set to 1.
    """
    (u, v) = fieldset.UV[time, particle.depth, particle.lat, particle.lon]
    if fieldset.active == 1:
        npp = fieldset.NPP[time, particle.depth, particle.lat, particle.lon]
        if (math.fabs(u) < 1e-14 and math.fabs(v) < 1e-14) or math.fabs(npp) < 1e-14:
            particle.beached = 1
    else: 
        if math.fabs(u) < 1e-14 and math.fabs(v) < 1e-14:
            particle.beached = 1


           
def UndoMove(particle, fieldset, time):
    """
    Send particle back to last position in case it is on land.
    If it is on land more than onland_max times in a row, it is deleted.
    """
    onland_max = 10  
    if particle.beached == 1:
        particle.beached = 0
        particle.lon = particle.prev_lon
        particle.lat = particle.prev_lat
        particle.onland += 1
        #
        if particle.onland > onland_max:
            print("Particle [%d] was deleted after beaching %f times in a row."%(particle.id,particle.onland))
            particle.delete()
    else:
        particle.onland = 0


        

def CheckOnLand(particle,fieldset,time):
    """ 
    Check if particles are released on land. 
    Executed once with dt=0 just after ParticleSet is created.
    """
    (u, v) = fieldset.UV[0, particle.depth, particle.lat, particle.lon]
    if fieldset.active == 1:
        npp = fieldset.NPP[time, particle.depth, particle.lat, particle.lon]
        if (math.fabs(u) < 1e-14 and math.fabs(v) < 1e-14) or math.fabs(npp) < 1e-14:
            print("Particle [%d] is released on land at lon,lat = %f,%f. Execution stops."%(particle.id,particle.lon,particle.lat))
            exit(0)
    else:
        if math.fabs(u) < 1e-14 and math.fabs(v) < 1e-14:
            print("Particle [%d] is released on land at lon,lat = %f,%f. Execution stops."%(particle.id,particle.lon,particle.lat))
            exit(0)
    
    
    
    
    
    
    
    
    
    
    
    
    