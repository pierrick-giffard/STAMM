#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Define advection kernels with added swimming velocity.
"""

import math


def RK4_swim(particle, fieldset, time):
    """
    Advection of particles using fourth-order Runge-Kutta integration.
    Swimming velocity is added to the current velocity.
    Swimming velocity is supposed constant during the whole timestep.
    """
    if particle.active == 1:
        #From m/s to °/s
        u_swim = particle.u_swim * math.cos(particle.lat * math.pi / 180) / fieldset.deg
        v_swim = particle.v_swim / fieldset.deg
        #   
        (u1, v1) = fieldset.UV[time, particle.depth, particle.lat, particle.lon]
        u1 += u_swim
        v1 += v_swim
        lon1, lat1 = (particle.lon + u1*.5*particle.dt, particle.lat + v1*.5*particle.dt)
    
        (u2, v2) = fieldset.UV[time + .5 * particle.dt, particle.depth, lat1, lon1]
        u2 += u_swim
        v2 += v_swim
        lon2, lat2 = (particle.lon + u2*.5*particle.dt, particle.lat + v2*.5*particle.dt)
    
        (u3, v3) = fieldset.UV[time + .5 * particle.dt, particle.depth, lat2, lon2]
        u3 += u_swim
        v3 += v_swim
        lon3, lat3 = (particle.lon + u3*particle.dt, particle.lat + v3*particle.dt)
    
        (u4, v4) = fieldset.UV[time + particle.dt, particle.depth, lat3, lon3]
        u4 += u_swim
        v4 += v_swim
    
        particle.lon += (u1 + 2*u2 + 2*u3 + u4) / 6. * particle.dt
        particle.lat += (v1 + 2*v2 + 2*v3 + v4) / 6. * particle.dt



def Euler_swim(particle, fieldset, time):
    """
    Advection of particles using Explicit Euler (aka Euler Forward) integration.
    Swimming velocity is added to the current velocity.
    """
    if particle.active == 1:
        #From m/s to °/s
        u_swim = particle.u_swim * math.cos(particle.lat * math.pi / 180) / fieldset.deg
        v_swim = particle.v_swim / fieldset.deg
        #
        (u1, v1) = fieldset.UV[time, particle.depth, particle.lat, particle.lon]
        u1 += u_swim
        v1 += v_swim
        particle.lon += u1 * particle.dt
        particle.lat += v1 * particle.dt
            
    


