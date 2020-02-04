#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Define kernels needed for computation.
3 functions:
    -define_advection_kernel
    -define_additional_kernels
    -sum_kernels
    
"""


from math import exp, sqrt, cos, sin, atan2
import math
from parcels import random




"""
SHARED KERNELS
"""
  
def compute_habitat(particle, fieldset, time):
    """
    Computes habitat at position, left, right, bottom and top.
    Computes habitat gradients.
    Save habT, habPP and hab at the particle location.
    Save xgradh and ygradh.
    """
    #Convert dx to lon and lat
    dx_lon =  particle.grad_dx * cos(particle.lat * math.pi / 180) / particle.deg
    dx_lat =  particle.grad_dx / particle.deg
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
    if T0[0] < 0. or T0[0] > 100:
        print("Incorrect temperature value at lon,lat = %f,%f: set to 0"%(particle.lon,particle.lat))
        T0[0] = 0
    if T0[1] < 0. or T0[1] > 100:
        print("Incorrect temperature value at lon,lat = %f,%f: set to 0"%(particle.lon-dx_lon,particle.lat))
        T0[1] = 0
    if T0[2] < 0. or T0[2] > 100:
        print("Incorrect temperature value at lon,lat = %f,%f: set to 0"%(particle.lon+dx_lon,particle.lat))
        T0[2] = 0
    if T0[3] < 0. or T0[3] > 100:
        print("Incorrect temperature value at lon,lat = %f,%f: set to 0"%(particle.lon,particle.lat-dx_lat))
        T0[3] = 0
    if T0[4] < 0. or T0[4] > 100:
        print("Incorrect temperature value at lon,lat = %f,%f: set to 0"%(particle.lon,particle.lat+dx_lat))
        T0[4] = 0
    
    """
    Temperature habitat
    """
    Topt = 24. - 0.21*sqrt(particle.M)
    Tmin = 24. - 1.05*sqrt(particle.M)
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
        print('WARNING: NPP out of range at lon,lat = %f,%f: set to 0'%(particle.lon,particle.lat))
        food_hab[0] = 0
    else:
        food_hab[0] = min(NPP0[0]/particle.PPmax,1)
    #
    if NPP0[1] < 0 or NPP0[1] > 100000:
        print('WARNING: NPP out of range at lon,lat = %f,%f: set to 0'%(particle.lon-dx_lon,particle.lat))
        food_hab[1] = 0
    else:
        food_hab[1] = min(NPP0[1]/particle.PPmax,1)
    #
    if NPP0[2] < 0 or NPP0[2] > 100000:
        print('WARNING: NPP out of range at lon,lat = %f,%f: set to 0'%(particle.lon+dx_lon,particle.lat))
        food_hab[2] = 0
    else:
        food_hab[2] = min(NPP0[2]/particle.PPmax,1)
    #
    if NPP0[3] < 0 or NPP0[3] > 100000:
        print('WARNING: NPP out of range at lon,lat = %f,%f: set to 0'%(particle.lon,particle.lat-dx_lat))
        food_hab[3] = 0
    else:
        food_hab[3] = min(NPP0[3]/particle.PPmax,1)
    #
    if NPP0[4] < 0 or NPP0[4] > 100000:
        print('WARNING: NPP out of range at lon,lat = %f,%f: set to 0'%(particle.lon,particle.lat+dx_lat))
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
    particle.xgradh = (h_right - h_left)/(2 * particle.grad_dx)
    particle.ygradh = (h_top - h_bot)/(2 * particle.grad_dx)
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



def compute_swimming_direction(particle, fieldset, time):
    """
    Compute particule.theta
    Theta has to be between 0 and 2*pi for random.vommises, 0 corresponding to east.
    """
    #Compute theta0
    theta0 = atan2(particle.ygradh,particle.xgradh) #return 0 in case xgradh=ygradh=0 (east !)
                                                    #but ok because vonmises becomes uniform without gradient
    if theta0 < 0:
        theta0 = 2 * math.pi + theta0 #theta0 has to be between 0 and 2*pi
    
    grad = sqrt(math.pow(particle.xgradh, 2) + math.pow(particle.ygradh, 2))
    
    #Compute theta
    particle.theta = random.vonmisesvariate(theta0,particle.alpha*grad)


  


def compute_swimming_velocity(particle, fieldset, time):
    """
    Compute particule.u_swim and particle.v_swim
    """
    particle.u_swim = particle.vmax * (1-particle.hab) * cos(particle.theta)
    particle.v_swim = particle.vmax * (1-particle.hab) * sin(particle.theta)
    
      
"""
LEATHERBACK
"""             
def SCL_leath(particle, fieldset, time):
    "Compute Straight Carapace Length"
    particle.SCL = 1.43*(1-exp(-0.226*(particle.age/365.+0.17)))  
            
            
def Mass_leath(particle, fieldset, time):
    "Compute mass"
    particle.M = 112.31*(particle.SCL)**2.86
    
    
def PPmax_leath(particle, fieldset, time):
    """
    Compute food threshold
    Ref : Jones, T.T., Bostrom, B.L.,  Hastings, M.D., Van Houtan K.S., Pauly, D., 2012. Resource requirements of the Pacific Leatherback Turtle Population
    """
    #Besoin annuels en valeur absolue

    PPmax = 266.80368*(((1-exp(-0.299*(particle.age/365.+0.17)))**(2.86-1))*(exp(-0.299*(particle.age/365.+0.17))))/(1-(1-exp(-0.299*(particle.age/365.+0.17)))**(2.86*0.0328))
    particle.PPmax = PPmax/(2835.24/particle.P0)
  
             
def vmax_leath(particle, fieldset, time):
    """
    Compute maximum speed at current age.
    Ref : Gaspar, P., Benson, S., Dutton, P., Réveillère, A., Jacob, G., Meetoo, C., Dehecq, A., Fossette, S., 2012. Oceanic dispersal of juvenile leatherback turtles: going beyond passive drift modeling. Marine Ecology Progress Series.
    """
    particle.vmax = particle.vscale*(particle.SCL**0.126)    
    

"""
GREEN
"""                
def SCL_green(particle, fieldset, time):
                particle.SCL = 0.983*(1-exp(-0.949*(particle.age/365.+0.074)))


                
def define_turtle_kernels(pset, param):
    """
    Function that defines additional kernel that will be used for computation.
    Parameters:
        -pset: ParticleSet
        -param: needs mode (active or passive) and species (leatherback, loggerhead or green)
    """
    #
    mode = param['mode']
    species = param['species']
    #
    kernels_list = [] 
    if mode == 'active':       
        if species == 'leatherback':
            kernels_list.append(SCL_leath)      
            kernels_list.append(Mass_leath)
            kernels_list.append(PPmax_leath)
            kernels_list.append(vmax_leath)
        
        elif species == 'green':
            kernels_list.append(SCL_green)      
            kernels_list.append(Mass_green)
            kernels_list.append(PPmax_green)
            kernels_list.append(vmax_green) 
            

        
        kernels_list.append(compute_habitat)
        kernels_list.append(compute_swimming_direction)
        kernels_list.append(compute_swimming_velocity)

        
    for k in range(len(kernels_list)):
        kernels_list[k]=pset.Kernel(kernels_list[k]) 
    return kernels_list







