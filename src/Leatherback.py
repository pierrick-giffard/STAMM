#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Kernels relative to Leatherback turtles.    
"""


from math import exp, sqrt, cos
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
    """
    Habitat gradient
    """ 
    particle.xgradh = (h_right - h_left)/(2 * fieldset.grad_dx)
    particle.ygradh = (h_top - h_bot)/(2 * fieldset.grad_dx)
    """
    Safety checks
    """ 
    #Ensure habitat is set to 0 on land
    # if fieldset.LandMask == 0:
    #         T_hab = 0
    #         food_hab = 0
    if particle.hab < 0:
        a=particle.habT
        b=particle.habPP
        c=particle.hab
        print(a,b,c)



            
def compute_SCL(particle, fieldset, time):
    "Compute Straight Carapace Length"
    particle.SCL = 1.43*(1-exp(-0.226*(particle.age/365.+0.17)))

            
            
def compute_Mass(particle, fieldset, time):
    "Compute mass"
    particle.M = 112.31*(particle.SCL)**2.86
    
    
def compute_PPmax(particle, fieldset, time):
    """
    Compute food threshold
    Ref : Jones, T.T., Bostrom, B.L.,  Hastings, M.D., Van Houtan K.S., Pauly, D., 2012. Resource requirements of the Pacific Leatherback Turtle Population
    """
    #Besoin annuels en valeur absolue

    PPmax = 266.80368*(((1-exp(-0.299*(particle.age/365.+0.17)))**(2.86-1))*(exp(-0.299*(particle.age/365.+0.17))))/(1-(1-exp(-0.299*(particle.age/365.+0.17)))**(2.86*0.0328))
    particle.PPmax = PPmax/(2835.24/fieldset.P0)
  
             
def compute_vmax(particle, fieldset, time):
    """
    Compute maximum speed at current age.
    Ref : Gaspar, P., Benson, S., Dutton, P., Réveillère, A., Jacob, G., Meetoo, C., Dehecq, A., Fossette, S., 2012. Oceanic dispersal of juvenile leatherback turtles: going beyond passive drift modeling. Marine Ecology Progress Series.
    """
    particle.vmax = fieldset.vscale*(particle.SCL**0.126)

    
def cold_induced_mortality(particle, fieldset, time):
    """
    Increment particle.lethargy_time if T < Tmin.
    If particle.lethargy_time > cold_resistance, then delete particle.
    """
    if fieldset.T[time, particle.depth, particle.lat, particle.lon] < particle.Tmin:
        particle.lethargy_time += fieldset.tstep
        if particle.lethargy_time > fieldset.cold_resistance:
            particle.cold_death = 1
            particle.delete()
    else:
        particle.lethargy_time = 0

        
        

        
        
        