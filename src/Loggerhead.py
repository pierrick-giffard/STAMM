#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Kernels relative to Loggerhead turtles.
Note: parcels kernels don't accept homemade functions, for loops or numpy function.    
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
        #
        #Check temperature values
        #
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
        
        #
        #Temperature habitat
        #
        Topt = 10.
        Tmin = 18.
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
        #Food habitat
        #
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


            
def compute_SCL_VGBF_age(particle, fieldset, time):
    """
    Compute Straight Carapace Length (meters) as a function of age. Age is in days.
    Uses a Von Bertalanffy function (VGBF).
    Ref : Jones, T.T., Hastings, M.D., Bostrom, B.L., Pauly, D., Jones, D.R., 2011. Growth of captive leatherback turtles, Dermochelys coriacea, with inferences on growth in the wild: Implications for population decline and recovery. Journal of Experimental Marine Biology and Ecology 399, 84–92.
    """
    if particle.active == 1:
        k = 0.0981
        SCLmax = 1.09
        t0 = -0.36 
        particle.SCL = SCLmax * (1 - exp(-k * (particle.age/365. - t0)))





def compute_SCL_VGBF(particle, fieldset, time):
    """
    Compute Straight Carapace Length (meters) at time t based on SCL at time t-1.
    Uses a Von Bertalanffy function (VGBF).
    Ref : Jones, T.T., Hastings, M.D., Bostrom, B.L., Pauly, D., Jones, D.R., 2011. Growth of captive leatherback turtles, Dermochelys coriacea, with inferences on growth in the wild: Implications for population decline and recovery. Journal of Experimental Marine Biology and Ecology 399, 84–92.
    """
    k = 0.0981
    SCLmax = 1.09 
    #
    L = particle.SCL + k * (SCLmax - particle.SCL) * particle.dt / 31536000 #dt has to be in years --> 86400*365
    particle.SCL = L





def compute_Mass(particle, fieldset, time):
    "Compute mass (kg). SCL is in meters."
    if particle.active == 1:
        a = 104
        b = 2.5
        particle.M = a*(particle.SCL)**b

    
    
def compute_PPmax_VGBF(particle, fieldset, time):
    """
    Compute food threshold
    Ref:
    """
    if particle.active == 1:
        b = 2.5
        beta = 0.078
        SCLmax = 1.09
        #
        x = particle.SCL / SCLmax
        PPnorm = b * beta * (1 - x) * x**(b-1) / (1 - x**(b * beta))
        particle.PPmax = PPnorm * fieldset.P0
    

  
             
def compute_vmax(particle, fieldset, time):
    """
    Compute maximum speed at current age.
    Ref : Gaspar, P., Benson, S., Dutton, P., Réveillère, A., Jacob, G., Meetoo, C., Dehecq, A., Fossette, S., 2012. Oceanic dispersal of juvenile leatherback turtles: going beyond passive drift modeling. Marine Ecology Progress Series.
    """
    if particle.active == 1:
        particle.vmax = fieldset.vscale*(particle.SCL**0.025)

    


        
        

        
        
        
