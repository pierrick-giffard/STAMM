#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Definition of kernels used for all particles (both in active and passive modes).
"""

import math

def store_variables(particle, fieldset, time):
    """
    Store lon/lat before the are re-written.
    Has to be the first kernel.
    """
    if particle.active:
        particle.prev_lon = particle.lon
        particle.prev_lat = particle.lat
 
       
def compute_SCL_VGBF(particle, fieldset, time):
    """
    Compute Straight Carapace Length (meters) at time t based on SCL at time t-1.
    Uses a Von Bertalanffy function (VGBF).
    """
    particle.SCL = particle.SCL + fieldset.k * (fieldset.SCLmax - particle.SCL) * particle.dt / 31536000 #dt has to be in years --> 86400*365



def compute_SCL_Gompertz(particle, fieldset, time):
    """
    Compute Straight Carapace Length (meters). Age is in days.
    Uses a modified Gompertz equation in which growth depends on habitat.
    This model needs SCL in cm.
    """
    if particle.active:
        prev_SCL = particle.SCL * 100 #this model needs centimeters
        prev_K = particle.K
        #
        SCL = prev_SCL + fieldset.alpha_gomp * particle.hab * log(prev_K / prev_SCL) * prev_K *  particle.dt / 86400
        particle.K = prev_K + fieldset.beta * particle.hab * 1 / (1 + exp(-(fieldset.M0 - prev_SCL) / fieldset.S)) *  particle.dt / 86400
        #
        particle.SCL = SCL / 100 #back to meters


def IncrementAge(particle, fieldset, time):
   """"
   Increment turtles age (in days).
   Deactivate turtles that have already lived 'ndays_simu' days.
   """
   if particle.active:
       particle.age += particle.dt / 86400.
       #
       if particle.age > fieldset.ndays_simu + 1:
           particle.active = 0



def SampleTracers(particle, fieldset, time):
    """
    Sample tracers at particle location in passive mode.
    In active mode sampling is integrated to function compute_habitat.
    """
    if particle.active:
        particle.T = fieldset.T[time, particle.depth, particle.lat, particle.lon]
        particle.NPP = fieldset.NPP[time, particle.depth, particle.lat, particle.lon]
        
        
def SampleCurrent(particle, fieldset, time):
    """
    Sample u_current and v_current in passive mode.
    In active mode, sampling is integrated to advection kernel.
    """
    if particle.active:
        uc, vc = fieldset.UV[time, particle.depth, particle.lat, particle.lon]
        particle.u_current = uc * fieldset.deg * math.cos(particle.lat * math.pi / 180) 
        particle.v_current = vc * fieldset.deg


def Periodic(particle, fieldset, time):
    """
    So that particles move from east boundary to west boundary.
    """
    if particle.active:
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
    if particle.active:
        (u, v) = fieldset.UV[time, particle.depth, particle.lat, particle.lon]
        if math.fabs(u) < 1e-14 and math.fabs(v) < 1e-14:
            particle.beached = 1
        else:
            particle.beached = 0



           
def UndoMove(particle, fieldset, time):
    """
    Send particle back to last position in case it is on land.
    If it is on land more than onland_max times in a row, it is deleted.
    The tactic factor is set to 1 during 1 time step when a particle is on land.
    """
    onland_max = 50
    if particle.active:        
        if particle.beached == 1:
            #print('Particle [%d] beached at lon,lat = %f,%f and time = %f'%(particle.id,particle.lon,particle.lat,particle.time))                       
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
    Time is taken at 41400 sec from origin (11.5 h) because of an issue with hourly forcings (time = 0 is OK with daily forcings).
    This value should be the minimum time of all particles.
    """
    (u, v) = fieldset.UV[41400, particle.depth, particle.lat, particle.lon]
    if math.fabs(u) < 1e-14 and math.fabs(v) < 1e-14:
        print("Particle [%d] is released on land at lon,lat = %f,%f. Execution stops."%(particle.id,particle.lon,particle.lat))
        exit(0)
        
    

    
    
    
