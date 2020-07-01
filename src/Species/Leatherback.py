#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Parameters for Leatherback sea turtles.  
"""


# =============================================================================
# SCL
# =============================================================================
### Von Bertalanffy (VGBF) ###
#Ref : Jones, T.T., Hastings, M.D., Bostrom, B.L., Pauly, D., Jones, D.R., 2011. Growth of captive leatherback turtles, Dermochelys coriacea, with inferences on growth in the wild: Implications for population decline and recovery. Journal of Experimental Marine Biology and Ecology 399, 84–92.
k = 0.226
SCLmax = 1.43 

### Gompertz modified (habitat dependant growth) ###
#Ref: D. Chevallier, B. Mourrain, M. Girondot, 2020. Modelling leatherback biphasic indeterminate growth using a modified Gompertz equation.  
alpha_gomp = 0.009839
beta = 0.084164
M0 = 105.2563
S = 15.5433
K0 = 20.13


# =============================================================================
# Mass
# =============================================================================
#M = a * SCL**b          
a = 112.31
b = 2.86
    
    
# =============================================================================
# Tmin and Topt (T1 and T2)
# ============================================================================= 
#Tmin_Topt: 'variable' or 'constant'. 
#If 'constant', needs Topt and Tmin. If 'variable', needs T0, to and tm and Topt, Tmin are calculated as:
#Topt = T0. - to * sqrt(M)
#Tmin = T0. - tm * sqrt(M)  
Tmin_Topt = 'variable' 
T0 = 22.
to = 0.21
tm = 1.05       

  
# =============================================================================
# PPmax
# =============================================================================              

### VGBF ###
#Ref : Jones, T.T., Bostrom, B.L.,  Hastings, M.D., Van Houtan K.S., Pauly, D., 2012. Resource requirements of the Pacific Leatherback Turtle Population
beta_jones = 0.0328

### Gompertz ###
#Assume PPmax = c * Mass
c = 0.0032#1/312        
        

# =============================================================================
# Vmax
# =============================================================================             
# vmax = vscale * SCL**d
#Ref : Gaspar, P., Benson, S., Dutton, P., Réveillère, A., Jacob, G., Meetoo, C., Dehecq, A., Fossette, S., 2012. Oceanic dispersal of juvenile leatherback turtles: going beyond passive drift modeling. Marine Ecology Progress Series.     
vscale = 1.2
d = 0.126


# =============================================================================
# Frenzy swimming
# =============================================================================
frenzy_theta = 0 # angle in rad defined with respect to east
frenzy_duration = 7 # duration of frenzy swimming in days
frenzy_speed = 0.25 # speed in m/s
