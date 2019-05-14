#!/usr/bin/python2.7
#-*- coding:utf-8 -*-

"""
	Author : Amaury Dehecq (intern at CLS 2011-2012 working with Phillipe Gaspar)
	Updated : Hirohiti Raapoto (intern at CLS 2013 working with Phillipe Gaspar)
	Module description : Individual Based Model for sea turtle populations

        To launch the program : > python IBM2D.py namelist outfile
        where :
        - namelist : str, name of the namelist file containing all parameters of the simulation. By default, the file must be situated in the 'input' directory.
        - outfile : str, name of the file where to save output. By default, the file will be saved in the 'output' directory.

        For more informations, read 'IBM-tutorial' file
   

"""


#Python libraries
import numpy as np
import time
import os, sys
from datetime import date as jour

#Personal libraries
import IOlib as IO
import Transports_vonMises_anticrash as ot
import TurtleClass as tc


#Read arguments
input_file = sys.argv[1]
output_file = sys.argv[2]
#input_dir = '../../input/'
#output_dir = '../../output/'

print(output_file)
print(input_file)
#Read input parameters in namelist file
param = IO.read_namelist(input_file)


"""
Input parameters for the simulation
"""
#Time step in seconds (imposed by the time step of input files U,V,T)
#For now, the program does not allow a sub-sampling of the position
tstep = param['tstep']
#Max number of steps for the simulation (max number of data file available)
nsteps_max = param['nsteps_max']
#Number of steps simulated
nsteps_simu = param['nsteps_simu']
#Number of trajectories to be computed (less or equal to number of lines in intital_positions file)
nturtles = param['nturtles']
#Boolean, True if all tracers are to be used
key_alltracers = param['key_alltracers']
#Mode can be 'passive' (no intended movement, the IBM works exactly as Ariane in this mode), 'diffusion' (passive + random difffusion) or 'active' if directed movement is to be modelled
mode = param['mode']


"""
Input grid parameters
"""
#Object of class Grid containing all parameters of the grid
Cgrid =IO.Grid(param)


"""
Read initial positions and time
"""
x_init, y_init, t_init = IO.read_positions(param,Cgrid.glamu,Cgrid.gphiv,Cgrid.mask)
t_init = np.floor(t_init)+0.5
#Starting and ending index of calendar
nstart = int(np.min(t_init))
nstop = int(param['nsteps_simu']+np.max(t_init))+1


"""
Initialisation of object of class TurtleClass, containing informations on all turtles at current time (position, ambient temperature, SCL...) 
"""
print("SPECIES", param['species'])
if param['species']=='leatherback':
    print("ok species")
    pop = tc.Leatherback(x_init, y_init, t_init, nsteps_simu, key_alltracers, mode)
elif param['species']=='green':
    pop = tc.Green(x_init, y_init, t_init, nsteps_simu, key_alltracers, mode) 

elif param['species']=='loggerhead':
    
    pop = tc.Loggerhead(x_init,y_init,t_init,nsteps_simu,key_alltracers,mode)
else:
    sys.exit("ERROR : wrong species")

pop.UpDatePosition(Cgrid)
pop.compute_SCL()


"""
Time variables
"""
#Current time for each turtle, in number of time steps
t=np.zeros(nturtles,np.float64)
t[:]=t_init[:]
t_init_min=np.min(t)


"""
Objects which gather every extern parameters
"""
#Zonal current velocity
U = IO.Variable(param['c_dir_zo'],param['c_prefix_zo'],param['c_suffix_zo'],param['maxsize_zo'],param['nc_var_zo'],param['nc_lon_zo'],param['nc_lat_zo'],param['nc_att_mask_zo'],Cgrid.umask,param['ind0_zo'])

#Meridional current velocity
V = IO.Variable(param['c_dir_me'],param['c_prefix_me'],param['c_suffix_me'],param['maxsize_me'],param['nc_var_me'],param['nc_lon_me'],param['nc_lat_me'],param['nc_att_mask_me'],Cgrid.vmask,param['ind0_me'])


#Optionnal tracers
if key_alltracers == True:
    #Sea surface temperature
    Temp = IO.Variable(param['c_dir_te'],param['c_prefix_te'],param['c_suffix_te'],param['maxsize_te'],param['nc_var_te'],param['nc_lon_te'],param['nc_lat_te'],param['nc_att_mask_te'],Cgrid.mask,param['ind0_te'])
    pop.UpDateTemperature(Temp)
    
    #Primary production (or micronekton)
    PP = IO.Variable(param['c_dir_pp'],param['c_prefix_pp'],param['c_suffix_pp'],param['maxsize_pp'],param['nc_var_pp'],param['nc_lon_pp'],param['nc_lat_pp'],param['nc_att_mask_pp'],Cgrid.mask,param['ind0_pp'])
    TotalTransport = ot.Transport(pop,U.values,V.values,Cgrid,mode,Temp.values,PP.values)
    pop.UpDatePP(PP)

 
    
else:
    Temp = 0
    PP = 0
    TotalTransport = ot.Transport(pop,U.values,V.values,Cgrid,mode)


"""
Loop on time
"""
#At the beginning, all turtles are inactive (state = 0)
pop.state = np.zeros(pop.nturtles,dtype='int32') 
physio_turtle = np.zeros((nsteps_simu+1,nturtles))
active_turtles = []

date = np.zeros((nsteps_simu+1,nturtles))	#

i=0			#

########### Compute fonction age Jones

for step in range(nstart-1,nstop-1):
    print("Step ", step)

    # Loop on time if not enough input files
    step = step - (2557*(int(step/2557.)))+1 

    """
    Loop on active turtles
    """
    
    for j in range(nturtles):
        #Activate turtle if current time > initial time (day start at 12am, turtles are released at 12pm -> +0.5). Turtle out of domain stay where they are
        if (pop.t_init[j]<step+0.5) & (pop.state[j]!=2) & (pop.state[j]!=3):
            pop.state[j]=1
        #De-activate turtle if number of simulated steps > max number of steps to be simulated 
        if ((t[j]-pop.t_init[j])>=nsteps_simu) & (pop.state[j]!=3):
            pop.state[j]=0


    #Update active turtles and incremente time
    active_turtles = np.where(pop.state==1)[0]
    t[active_turtles]=t[active_turtles]+1.

    #Update turtles 'out of domain' and incremente time
    turtles_out = np.where(pop.state==2)[0]					
    t[turtles_out]=t[turtles_out]+1.						
    pop.age_now[turtles_out] = pop.age_now[turtles_out]+tstep/86400		
    pop.age[pop.tind[turtles_out],turtles_out] = pop.age_now[turtles_out]	

    if (len(active_turtles)>0):
        # 1998-01-01
        origin=729390
        # 2002-01-01
        #origin=730851
        print("Day ", jour.fromordinal(origin-1 + step))
        t1=time.clock()

        """
        Update external parameters
        """
        U.UpdateDay(param['ind0_zo']+step-1)
        V.UpdateDay(param['ind0_me']+step-1)
        if key_alltracers == True:
            Temp.UpdateDay(param['ind0_te']+step-1)
            PP.UpdateDay(param['ind0_pp']+step-1)

        #Update turtle's caracteristics
        pop.compute_SCL()
        pop.compute_M()
        if param['mode'] == 'active':
            pop.compute_PPmax(param['Fa'])
            pop.compute_vmax(param['vscale'])
        #Update transport
        if key_alltracers == True:
            TotalTransport.UpDate(pop,U.values,V.values,Temp.values,PP.values)
        else:
            TotalTransport.UpDate(pop,U.values,V.values)


        TotalTransport.compute_TotalTransport(param['alpha'],param['c_prefix_pp'])
        TotalTransport.Advection2D(tstep)   
        """
        Update position and caracteristics of active turtles
        """
        #Save initial temperature for first loop
        if key_alltracers == True:
            for j in range(nturtles):
                if (pop.t_init[j]==step+0.5):
                    pop.Tenv[0,j] = Temp.values[pop.j0[j],pop.i0[j]]

        #Update position and tracers
        pop.UpDatePosition(Cgrid)
        if key_alltracers == True:
            pop.UpDateTemperature(Temp)
            pop.UpDatePP(PP)

        t2=time.clock()
        print('Done in '+str(t2-t1)+' sec. \n') 
    """
    Compute date instead of turtles age
    """
    if i<(nsteps_simu+1):
        print(i)
        for a in range(nturtles):
            date[i,a] = t_init[a]+i+0.5
        i=i+1
	

"""
Save output :
- initial positions (grid indices) and times : x_init, y_init, t_init
- final positions (grid indices) and times
- positions (longitude, latitude) and age at all time
- positions (grid indices) at all time
- tracers along trajectories at all time (temperature, primary production...)
- if in 'active' mode, habitat along trajectories at all time
"""

print('*******************')
print('==> Saving results ')
print('*******************')

title=output_file
IO.FillNetcdfFile(param,title,x_init,y_init,t_init,pop.x,pop.y,t,pop.trajlon,pop.trajlat,pop.age,date,pop.traji0,pop.trajj0,pop.u_current,pop.v_current,pop.Tenv,pop.PP,pop.habT,pop.habPP,pop.hab,pop.xgrad,pop.ygrad,pop.u_swim,pop.v_swim,pop.u_tot,pop.v_tot, pop.SCL)

print('\n')
print('****************************************************')
print('==> Results saved in : \n'+str(title))
print('****************************************************')


