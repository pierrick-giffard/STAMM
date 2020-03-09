#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Kernels relative to Leatherback turtles.
Note: parcels kernels don't accept homemade functions, for loopsor numpy function.    
"""


from math import exp, sqrt, cos, log
import math


  
def compute_habitat(particle, fieldset, time):
    """
    Computes habitat at position, left, right, bottom and top.
    Computes habitat gradients.
    Save habT, habPP and hab at the particle location.
    Save xgradh and ygradh.
    """
    #Convert dx to lon and lat
    dx_lon =  fieldset.grad_dx * cos(particle.lat * math.pi / 180) / fieldset.deg
    dx_lat =  fieldset.grad_dx / fieldset.deg
    #
    """
    Get 5 T and 5 NPP
    """
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
    """
    Check temperature values
    """
    if T0[0] < -10. or T0[0] > 100:
        print("WARNING: Incorrect temperature at lon,lat = %f,%f and time = %f: set to 0"%(particle.lon,particle.lat,time))
        T0[0] = 0
    if T0[1] < -10. or T0[1] > 100:
        print("WARNING: Incorrect temperature at lon,lat = %f,%f and time = %f: set to 0"%(particle.lon,particle.lat,time))
        T0[1] = 0
    if T0[2] < -10. or T0[2] > 100:
        print("WARNING: Incorrect temperature at lon,lat = %f,%f and time = %f: set to 0"%(particle.lon,particle.lat,time))
        T0[2] = 0
    if T0[3] < -10. or T0[3] > 100:
        print("WARNING: Incorrect temperature at lon,lat = %f,%f and time = %f: set to 0"%(particle.lon,particle.lat,time))
        T0[3] = 0
    if T0[4] < -10. or T0[4] > 100:
        print("WARNING: Incorrect temperature at lon,lat = %f,%f and time = %f: set to 0"%(particle.lon,particle.lat,time))
        T0[4] = 0
    
    """
    Temperature habitat
    """
    Topt = 24. - 0.21*sqrt(particle.M)
    Tmin = 24. - 1.05*sqrt(particle.M)
    particle.Tmin = Tmin
    #
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
    """
    Food habitat
    """
    food_hab = [0, 0, 0, 0, 0] #position, left, right, bottom and top
    #
    if NPP0[0] < 0 or NPP0[0] > 100000:
        print('WARNING: NPP out of range at lon,lat = %f,%f and time = %f: set to 0'%(particle.lon,particle.lat,time))
        food_hab[0] = 0
    else:
        food_hab[0] = min(NPP0[0]/particle.PPmax,1)
    #
    if NPP0[1] < 0 or NPP0[1] > 100000:
        print('WARNING: NPP out of range at lon,lat = %f,%f and time = %f: set to 0'%(particle.lon,particle.lat,time))
        food_hab[1] = 0
    else:
        food_hab[1] = min(NPP0[1]/particle.PPmax,1)
    #
    if NPP0[2] < 0 or NPP0[2] > 100000:
        print('WARNING: NPP out of range at lon,lat = %f,%f and time = %f: set to 0'%(particle.lon,particle.lat,time))
        food_hab[2] = 0
    else:
        food_hab[2] = min(NPP0[2]/particle.PPmax,1)
    #
    if NPP0[3] < 0 or NPP0[3] > 100000:
        print('WARNING: NPP out of range at lon,lat = %f,%f and time = %f: set to 0'%(particle.lon,particle.lat,time))
        food_hab[3] = 0
    else:
        food_hab[3] = min(NPP0[3]/particle.PPmax,1)
    #
    if NPP0[4] < 0 or NPP0[4] > 100000:
        print('WARNING: NPP out of range at lon,lat = %f,%f and time = %f: set to 0'%(particle.lon,particle.lat,time))
        food_hab[4] = 0
    else:
        food_hab[4] = min(NPP0[4]/particle.PPmax,1)
    #
    """
    Total habitat
    """
    particle.habT = T_hab[0]
    particle.habPP = food_hab[0]
    particle.hab = particle.habT * particle.habPP
    h_left = T_hab[1] * food_hab[1]
    h_right = T_hab[2] * food_hab[2]
    h_bot = T_hab[3] * food_hab[3]
    h_top = T_hab[4] * food_hab[4]
    #
    """
    Habitat gradient
    """ 
    particle.xgradh = (h_right - h_left)/(2 * fieldset.grad_dx)
    particle.ygradh = (h_top - h_bot)/(2 * fieldset.grad_dx)
    #
    """
    Safety check
    """
    if particle.hab < 0 or particle.hab > 1:
        print("Habitat is %f at lon,lat = %f,%f. Execution stops."%(particle.hab,particle.lon,particle.lat))
        exit(0)


            
def compute_SCL_VGBF(particle, fieldset, time):
    """
    Compute Straight Carapace Length (meters). Age is in days.
    Uses a Von Bertalanffy function (VGBF).
    Ref : Jones, T.T., Hastings, M.D., Bostrom, B.L., Pauly, D., Jones, D.R., 2011. Growth of captive leatherback turtles, Dermochelys coriacea, with inferences on growth in the wild: Implications for population decline and recovery. Journal of Experimental Marine Biology and Ecology 399, 84–92.
    """
    k = 0.226
    SCLmax = 1.43
    t0 = -0.17    
    particle.SCL = SCLmax * (1 - exp(-k * (particle.age/365. - t0)))



def compute_SCL_Gompertz(particle, fieldset, time):
    """
    Compute Straight Carapace Length (meters). Age is in days.
    Uses a modified Gompertz equation in which growth depends on habitat.
    This model needs SCL in cm.
    Ref: D. Chevallier, B. Mourrain, M. Girondot, 2020. Modelling leatherback biphasic indeterminate growth using a modified Gompertz equation.
    """
    alpha0 = 0.009839
    beta0 = 0.084164
    M0 = 105.2563
    S0 = 15.5433
    #
    prev_SCL = particle.SCL * 100 #this model needs centimeters
    prev_K = particle.K
    #
    SCL = prev_SCL + alpha0 * particle.hab * log(prev_K / prev_SCL) * prev_K *  particle.dt / 86400
    particle.K = prev_K + beta0 * particle.hab * 1 / (1 + exp(-(M0 - prev_SCL) / S0)) *  particle.dt / 86400
    #
    particle.SCL = SCL / 100 #back to meters
    
    

            
def compute_Mass(particle, fieldset, time):
    "Compute mass (kg). SCL is in meters."
    a = 112.31
    b = 2.86
    particle.M = a*(particle.SCL)**b

    
    
def compute_PPmax_VGBF(particle, fieldset, time):
    """
    Compute food threshold
    Ref : Jones, T.T., Bostrom, B.L.,  Hastings, M.D., Van Houtan K.S., Pauly, D., 2012. Resource requirements of the Pacific Leatherback Turtle Population
    """
    b = 2.86
    beta = 0.0328
    SCLmax = 1.43
    #
    x = particle.SCL / SCLmax
    PPnorm = b * beta * (1 - x) * x**(b-1) / (1 - x**(b * beta))
    particle.PPmax = PPnorm * fieldset.P0
    

def compute_PPmax_Gompertz(particle, fieldset, time):
    """
    Assume F = c*M
    """
    c = 0.0031#1/325
    PPnorm = c * particle.M
    particle.PPmax = PPnorm * fieldset.P0

  
             
def compute_vmax(particle, fieldset, time):
    """
    Compute maximum speed at current age.
    Ref : Gaspar, P., Benson, S., Dutton, P., Réveillère, A., Jacob, G., Meetoo, C., Dehecq, A., Fossette, S., 2012. Oceanic dispersal of juvenile leatherback turtles: going beyond passive drift modeling. Marine Ecology Progress Series.
    """
    particle.vmax = fieldset.vscale*(particle.SCL**0.126)

    


        
        

        
        
        