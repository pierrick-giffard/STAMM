#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Kernels used to compute swimming velocity and cold induced mortality.
"""


from math import sqrt, cos, sin, atan2, exp, log
import math
from parcels import random


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
    if particle.active == 1:
        prev_SCL = particle.SCL * 100 #this model needs centimeters
        prev_K = particle.K
        #
        SCL = prev_SCL + fieldset.alpha_gomp * particle.hab * log(prev_K / prev_SCL) * prev_K *  particle.dt / 86400
        particle.K = prev_K + fieldset.beta * particle.hab * 1 / (1 + exp(-(fieldset.M0 - prev_SCL) / fieldset.S)) *  particle.dt / 86400
        #
        particle.SCL = SCL / 100 #back to meters



def compute_Mass(particle, fieldset, time):
    "Compute mass (kg). SCL is in meters."
    if particle.active == 1:
        particle.M = fieldset.a*(particle.SCL)**fieldset.b


def compute_Tmin_Topt(particle, fieldset, time): 
    if particle.active == 1:
        particle.Topt = fieldset.T0 - fieldset.to * sqrt(particle.M)
        particle.Tmin = fieldset.T0 - fieldset.tm * sqrt(particle.M)



def compute_PPmax_VGBF(particle, fieldset, time):
    """
    Compute food threshold.
    """
    if particle.active == 1:
        x = particle.SCL / fieldset.SCLmax
        PPnorm = fieldset.b * fieldset.beta_jones * (1 - x) * x**(fieldset.b-1) / (1 - x**(fieldset.b * fieldset.beta_jones))
        particle.PPmax = PPnorm * fieldset.P0
    

def compute_PPmax_Gompertz(particle, fieldset, time):
    """
    Assume F = c*M
    """
    if particle.active == 1:
        PPnorm = fieldset.c * particle.M
        particle.PPmax = PPnorm * fieldset.P0


def compute_vmax(particle, fieldset, time):
    """
    Compute maximum speed at current age.
    """
    if particle.active == 1:
        particle.vmax = fieldset.vscale*(particle.SCL**fieldset.d)



def compute_habitat(particle, fieldset, time):
    """
    Computes habitat at position, left, right, bottom and top.
    Computes habitat gradients.
    Save habT, habPP and hab at the particle location.
    Save xgradh and ygradh.
    Save T and NPP at particle location.
    """
    if particle.active == 1:
        #Convert dx to lon and lat
        dx_lon =  fieldset.grad_dx * cos(particle.lat * math.pi / 180) / fieldset.deg
        dx_lat =  fieldset.grad_dx / fieldset.deg
        #
        #Get 5 T and 5 NPP
        #
        T0 = [fieldset.T[time, particle.depth, particle.lat, particle.lon],#position
               fieldset.T[time, particle.depth, particle.lat, particle.lon - dx_lon],#left
               fieldset.T[time, particle.depth, particle.lat, particle.lon + dx_lon],#right
               fieldset.T[time, particle.depth, particle.lat - dx_lat, particle.lon],#bottom
               fieldset.T[time, particle.depth, particle.lat + dx_lat, particle.lon]]#top
        
        NPP0 = [fieldset.NPP[time, particle.depth, particle.lat, particle.lon],#position
              fieldset.NPP[time, particle.depth, particle.lat, particle.lon - dx_lon],#left
              fieldset.NPP[time, particle.depth, particle.lat, particle.lon + dx_lon],#right
              fieldset.NPP[time, particle.depth, particle.lat - dx_lat, particle.lon],#bottom
              fieldset.NPP[time, particle.depth, particle.lat + dx_lat, particle.lon]]#top       
        #Save T and NPP at particle location
        particle.T = T0[0]
        particle.NPP = NPP0[0]
        #
        #Temperature habitat
        #
        Tmin = particle.Tmin
        Topt = particle.Topt
        T_hab = [0, 0, 0, 0, 0] #position, left, right, bottom and top
        #
        if T0[0] >= Topt:
            T_hab[0] = 1.0
        else:
            T_hab[0] = exp(-2*((T0[0]-Topt)/(Topt-Tmin))**2)
        #
        if T0[1] >= Topt:
            T_hab[1] = 1.0
        else:
            T_hab[1] = exp(-2*((T0[1]-Topt)/(Topt-Tmin))**2)
        #
        if T0[2] >= Topt:
            T_hab[2] = 1.0
        else:
            T_hab[2] = exp(-2*((T0[2]-Topt)/(Topt-Tmin))**2)
        #
        if T0[3] >= Topt:
            T_hab[3] = 1.0
        else:
            T_hab[3] = exp(-2*((T0[3]-Topt)/(Topt-Tmin))**2)
        #
        if T0[4] >= Topt:
            T_hab[4] = 1.0
        else:
            T_hab[4] = exp(-2*((T0[4]-Topt)/(Topt-Tmin))**2)
        #
        #Food habitat
        #
        food_hab = [0, 0, 0, 0, 0] #position, left, right, bottom and top
        #
        if NPP0[0] < 0:
            print('WARNING: negative NPP at lon,lat = %f,%f and time = %f: set to 0'%(particle.lon,particle.lat,time))
            food_hab[0] = 0
        else:
            food_hab[0] = min(NPP0[0]/particle.PPmax,1)
        #
        if NPP0[1] < 0:
            print('WARNING: negative NPP at lon,lat = %f,%f and time = %f: set to 0'%(particle.lon,particle.lat,time))
            food_hab[1] = 0
        else:
            food_hab[1] = min(NPP0[1]/particle.PPmax,1)
        #
        if NPP0[2] < 0:
            print('WARNING: negative NPP at lon,lat = %f,%f and time = %f: set to 0'%(particle.lon,particle.lat,time))
            food_hab[2] = 0
        else:
            food_hab[2] = min(NPP0[2]/particle.PPmax,1)
        #
        if NPP0[3] < 0:
            print('WARNING: negative NPP at lon,lat = %f,%f and time = %f: set to 0'%(particle.lon,particle.lat,time))
            food_hab[3] = 0
        else:
            food_hab[3] = min(NPP0[3]/particle.PPmax,1)
        #
        if NPP0[4] < 0:
            print('WARNING: negative NPP at lon,lat = %f,%f and time = %f: set to 0'%(particle.lon,particle.lat,time))
            food_hab[4] = 0
        else:
            food_hab[4] = min(NPP0[4]/particle.PPmax,1)
        #
        #Total habitat
        #
        particle.habT = T_hab[0]
        particle.habPP = food_hab[0]
        particle.hab = particle.habT * particle.habPP
        h_left = T_hab[1] * food_hab[1]
        h_right = T_hab[2] * food_hab[2]
        h_bot = T_hab[3] * food_hab[3]
        h_top = T_hab[4] * food_hab[4]
        #
        #Habitat gradient
        #
        particle.xgradh = (h_right - h_left)/(2 * fieldset.grad_dx)
        particle.ygradh = (h_top - h_bot)/(2 * fieldset.grad_dx)
        #
        #Safety check
        #
        if particle.hab < 0 or particle.hab > 1:
            print("Habitat is %f at lon,lat = %f,%f. Execution stops."%(particle.hab,particle.lon,particle.lat))
            exit(0)



def compute_swimming_direction(particle, fieldset, time):
    """
    Compute particule.theta
    Theta has to be between 0 and 2*pi for random.vommises, 0 corresponding to east.
    A tactic factor is used (t = [0,1]) as a memory effect (Benhamou 1991, Elementary orientation mechanisms).
    """
    if particle.active == 1:
        #Compute theta0
        theta0 = atan2(particle.ygradh,particle.xgradh) #return 0 in case xgradh=ygradh=0 (east !)
                                                        #but ok because vonmises becomes uniform without gradient
        if theta0 < 0:
            theta0 += 2 * math.pi #theta0 has to be between 0 and 2*pi
        
        grad = sqrt(math.pow(particle.xgradh, 2) + math.pow(particle.ygradh, 2))
        
        #Compute theta
        prev_theta = particle.theta
        current_theta = random.vonmisesvariate(theta0,fieldset.alpha*grad)
        particle.theta = atan2((fieldset.t * sin(current_theta) + (1 - fieldset.t) * sin(prev_theta)) , (fieldset.t * cos(current_theta) + (1 - fieldset.t) * cos(prev_theta)))  
        #
        if particle.theta < 0:
            particle.theta += 2 * math.pi #theta0 has to be between 0 and 2*pi



def compute_swimming_velocity(particle, fieldset, time):
    """
    Compute particule.u_swim and particle.v_swim
    """
    if particle.active == 1:
        particle.u_swim = particle.vmax * (1-particle.hab) * cos(particle.theta)
        particle.v_swim = particle.vmax * (1-particle.hab) * sin(particle.theta)
    


def cold_induced_mortality(particle, fieldset, time):
    """
    Increment particle.lethargy_time if T < Tmin.
    If particle.lethargy_time > cold_resistance, then delete particle.
    PB: how to keep in memory dead turtles ?
    """
    if particle.active == 1:
        if fieldset.T[time, particle.depth, particle.lat, particle.lon] < particle.Tmin:
            particle.lethargy_time += particle.dt
            if particle.lethargy_time > fieldset.cold_resistance:
                particle.cold_death = 1
                particle.active = 0
        else:
            particle.lethargy_time = 0

 
               








