#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
All kernels needed.
"""
from parcels import FieldSet, ParticleSet, JITParticle, ParticleFile, plotTrajectoriesFile, Variable,ErrorCode,Field
import numpy as np
from parcels import AdvectionRK4
from parcels import AdvectionEE as Euler

def define_advection_kernel(pset, mode, adv_scheme):
    """
    Function that defines the kernel that will be used for advection.
    Parameters:
        -mode: active or passive
        -adv_scheme: RK4 or Euler
    Return: the advection kernel (kernel object)
    """
    """
        PASSIVE MODE
    """
    if mode == 'passive':
        if adv_scheme == 'RK4':
            adv_kernel = AdvectionRK4
        elif adv_scheme == 'Euler':
            adv_kernel = Euler
            
    """
        ACTIVE MODE
    """        
    if mode == 'active':
        if adv_scheme == 'RK4':
            def RK4(particle, fieldset, time):
                """
                    Advection of particles using fourth-order Runge-Kutta integration.
                    Swimming velocity is added to the current velocity.
                    Swimming velocity is supposed constant during the whole timestep.
                """
                (u1, v1) = fieldset.UV[time, particle.depth, particle.lat, particle.lon]
                u1 += particle.u_swim
                v1 += particle.v_swim
                lon1, lat1 = (particle.lon + u1*.5*particle.dt, particle.lat + v1*.5*particle.dt)
            
                (u2, v2) = fieldset.UV[time + .5 * particle.dt, particle.depth, lat1, lon1]
                u2 += particle.u_swim
                v2 += particle.v_swim
                lon2, lat2 = (particle.lon + u2*.5*particle.dt, particle.lat + v2*.5*particle.dt)
            
                (u3, v3) = fieldset.UV[time + .5 * particle.dt, particle.depth, lat2, lon2]
                u3 += particle.u_swim
                v3 += particle.v_swim
                lon3, lat3 = (particle.lon + u3*particle.dt, particle.lat + v3*particle.dt)
            
                (u4, v4) = fieldset.UV[time + particle.dt, particle.depth, lat3, lon3]
                u4 += particle.u_swim
                v4 += particle.v_swim
            
                particle.lon += (u1 + 2*u2 + 2*u3 + u4) / 6. * particle.dt
                particle.lat += (v1 + 2*v2 + 2*v3 + v4) / 6. * particle.dt
            adv_kernel = RK4
            
        elif adv_scheme == 'Euler':
            def AdvectionEE(particle, fieldset, time):
                """
                    Advection of particles using Explicit Euler (aka Euler Forward) integration.
                    Swimming velocity is added to the current velocity.
                """
                (u1, v1) = fieldset.UV[time, particle.depth, particle.lat, particle.lon]
                u1 += particle.u_swim
                v1 += particle.v_swim
                particle.lon += u1 * particle.dt
                particle.lat += v1 * particle.dt
            adv_kernel = AdvectionEE
    
    return pset.Kernel(adv_kernel)


def define_additional_kernels(pset, mode, species, adv_scheme):
    """
        Function that defines additional kernel that will be used for computation.
        mode: active or passive
        species: leatherback or loggerhead
        adv_scheme: RK4 or Euler
    """
    if mode == 'active':  
       
        def AgeParticle(particle, fieldset, time, dt):
            particle.age += dt
        
        
        def SampleP(particle, fieldset, time): 
            particle.T = fieldset.T[time, particle.depth, particle.lat, particle.lon]
        
        """
            LEATHERBACK
        """ 
        if species == 'leatherback':    
            def compute_SCL(particle, fieldset, time):
                if mode == 'active':
                    if species == 'leatherback':
                        particle.SCL = 1.43*(1-np.exp(-0.226*(particle.age/365.+0.17)))
                    elif species == 'green':
                        particle.SCL = 0.983*(1-np.exp(-0.949*(particle.age/365.+0.074)))
                    elif species == 'loggerhead':
                        particle.SCL = 0
                
            def compute_M(particle, fieldset, time):
                if mode == 'active':
                    if species == 'leatherback':
                        particle.M = 112.31*(particle.SCL)**2.86
                    else:
                        particle.M = 0


def sum_kernels(adv_kernel, additional_kernels):
    return adv_kernel