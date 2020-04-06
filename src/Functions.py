#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Functions needed to:
    -create FieldSet
    -create ParticleSet
    -define kernels to be used in execute
    -modify output file
"""

# =============================================================================
# IMPORTS
# =============================================================================
#Python libraries
from parcels import FieldSet, ParticleSet, AdvectionRK4, AdvectionEE
from datetime import timedelta
import numpy as np
import time
import netCDF4
import math
import subprocess
import sys, os

#Personal libraries
import Advection_kernel as adv
import Passive_kernels as pk
import Active_kernels as ak
sys.path.insert(1, os.path.join(sys.path[0], 'Species'))
import Leatherback as leath
import Loggerhead as log
sys.path.insert(1, os.path.join(sys.path[0], '../LIB'))
import IOlib as IO


# =============================================================================
# FIELDSET
# =============================================================================
def create_fieldset(param, ndays_simu, t_init):
    """
    Build fieldset from data files, dimensions and variables.
    """
    print('****************************************************')
    print("Building FieldSet...")
    print('****************************************************')
    t0 = time.time()
    key_alltracers = param['key_alltracers']
    mesh_phy = param['mesh_phy']
    mesh_food = param['mesh_food']
    #Forcings
    last_date = IO.find_last_date(param)
    date_start, date_end, time_periodic = IO.define_start_end(ndays_simu, param, t_init, last_date)
    ufiles = IO.forcing_list(param['U_dir'], param['U_suffix'], date_start, date_end, print_date=True)
    vfiles = IO.forcing_list(param['V_dir'], param['V_suffix'], date_start, date_end)
    #Filenames
    filenames = {'U': {'lon': mesh_phy, 'lat': mesh_phy, 'data': ufiles},
                 'V': {'lon': mesh_phy, 'lat': mesh_phy, 'data': vfiles}}
    #Variables
    variables = {'U': param['U_var'],
                 'V': param['V_var']}
    #Dimensions: Caution, C-grids need f nodes
    dimensions = {'U': {'lon': param['lon_phy'], 'lat': param['lat_phy'], 'time': param['time_var_phy']},
                  'V': {'lon': param['lon_phy'], 'lat': param['lat_phy'], 'time': param['time_var_phy']}} 
    if key_alltracers:
        tfiles = IO.forcing_list(param['T_dir'], param['T_suffix'], date_start, date_end)
        ffiles = IO.forcing_list(param['food_dir'], param['food_suffix'], date_start, date_end, vgpm=param['vgpm'])
        #
        filenames['T'] = {'lon': mesh_phy, 'lat': mesh_phy, 'data': tfiles}
        filenames['NPP'] = {'lon': mesh_food, 'lat': mesh_food, 'data': ffiles}
        #
        variables['T'] = param['T_var']
        variables['NPP'] = param['food_var']
        #
        dimensions['T'] = {'lon': param['lon_phy'], 'lat': param['lat_phy'], 'time': param['time_var_phy']}
        dimensions['NPP'] = {'lon': param['lon_food'], 'lat': param['lat_food'], 'time': param['time_var_food']}
 
    #Fieldset creation
    if time_periodic:
        fieldset = FieldSet.from_netcdf(filenames, variables, dimensions, time_periodic=timedelta(days=time_periodic))
    else:
        fieldset = FieldSet.from_netcdf(filenames, variables, dimensions)#, field_chunksize=(800, 800))
    #East/West periodicity
    if param['periodicBC'] and not param['halo']:
        add_halo(fieldset)

    #Time
    tt=time.time()-t0
    print('\n')
    print(' => FieldSet created in: '+ str(timedelta(seconds=int(tt))))
    print('\n')
    return fieldset




def PSY_patch(fieldset,param):
    """
    Our files PSY4 interpolated on A-grid coarsened to 1/4° have a 
    problem at the equator and at Greenwich meridian: lon and lat are NaN.
    This function passes them to 0.
    """
    try:
        fieldset.U.grid.lon[:,720] = 0
        fieldset.U.grid.lat[320,:] = 0
        fieldset.V.grid.lon[:,720] = 0
        fieldset.V.grid.lat[320,:] = 0
        print('Using PSY patch')
        if param['key_alltracers']:
            fieldset.T.grid.lon[:,720] = 0
            fieldset.T.grid.lat[320,:] = 0
            fieldset.NPP.grid.lon[:,720] = 0
            fieldset.NPP.grid.lat[320,:] = 0
    except:
        print('Not Using PSY patch')


def add_halo(fieldset):
    try:
        fieldset.add_constant('halo_west', fieldset.U.grid.lon[0,0])
        fieldset.add_constant('halo_east', fieldset.U.grid.lon[0,-1])
    except:
        fieldset.add_constant('halo_west', fieldset.U.grid.lon[0])
        fieldset.add_constant('halo_east', fieldset.U.grid.lon[-1])
    #          
    fieldset.add_periodic_halo(zonal=True) 

# =============================================================================
# PARTICLESET
# =============================================================================
def create_particleset(fieldset, pclass, lon, lat, t_init, param):
    print('\n')
    print('****************************************************')
    print("Building ParticleSet...")
    print('****************************************************')
    #
    t0 = time.time()
    #
    t_release = (t_init - int(np.min(t_init))) * 86400
    pset = ParticleSet(fieldset, pclass=pclass, lon=lon, lat=lat, time = t_release)
    #
    pset.execute(pk.CheckOnLand, dt=0)
    #Time
    tt=time.time()-t0
    print('\n')
    print(' => ParticleSet created in: '+ str(timedelta(seconds=int(tt))))
    print('\n')
    
    return pset


# =============================================================================
# CONSTANT PARAMETERS
# =============================================================================
def initialization(fieldset, ndays_simu, param):
    """
    Links constant parameters to fieldset in order to use them within kernels.
    """
    fieldset.deg = 111195 #1degree = 111,195 km approx
    fieldset.cold_resistance = param['cold_resistance'] * 86400 #from days to seconds
    fieldset.ndays_simu = ndays_simu
    if param['mode'] == 'active':
        ### NAMELIST PARAMETERS ###
        fieldset.active = 1
        fieldset.vscale = param['vscale']
        fieldset.P0 = param['P0']
        fieldset.grad_dx = param['grad_dx']
        fieldset.alpha = param['alpha']
        ### SPECIES PARAMETERS ###
        if param['species'] == 'leatherback':
            file = leath
        elif param['species'] == 'loggerhead':
            file = log
        #
        fieldset.a = file.a
        fieldset.b = file.b
        fieldset.d = file.d
        fieldset.SCL0 = file.SCL0
        if param['growth'] == 'VGBF':
            fieldset.k = file.k
            fieldset.SCLmax = file.SCLmax
            fieldset.beta_jones = file.beta_jones
        elif param['growth'] == 'Gompertz':
            fieldset.alpha = file.alpha
            fieldset.beta = file.beta
            fieldset.M0 = file.M0
            fieldset.S = file.S
            fieldset.K0 = file.K0
            fieldset.c = file.c
        if file.Tmin_Topt == 'constant':
            fieldset.Tmin = file.Tmin
            fieldset.Topt = file.Topt
        elif file.Tmin_Topt == 'variable':
            fieldset.T0 = file.T0
            fieldset.to = file.to
            fieldset.tm = file.tm
            fieldset.Tmin = 0.
            fieldset.Topt = 0.
        else:
            raise ValueError('Please set Tmin_Topt to constant or variable')
        param['Tmin_Topt'] = file.Tmin_Topt
    else:
        fieldset.active = 0
    if param['key_alltracers']:
        fieldset.key_alltracers = 1
    else:
        fieldset.key_alltracers = 0



# =============================================================================
# DEFINE KERNELS
# =============================================================================
def define_advection_kernel(pset, param):
    """
    Function that defines the kernel that will be used for advection.
    Parameters:
        -pset: ParticleSet
        -param: needs mode (active or passive) and adv_scheme (RK4 or Euler)
    Return: the advection kernel (kernel object)
    """
    mode = param['mode']
    adv_scheme = param['adv_scheme']
    #passive
    if mode == 'passive':
        if adv_scheme == 'RK4':
            adv_kernel = AdvectionRK4
        elif adv_scheme == 'Euler':
            adv_kernel = AdvectionEE
    #active
    elif mode == 'active':
        if adv_scheme == 'RK4':
            adv_kernel = adv.RK4_swim 
        elif adv_scheme == 'Euler':
            adv_kernel = adv.Euler_swim
    
    return pset.Kernel(adv_kernel)



def define_passive_kernels(fieldset, pset, param):
    """
    Parameters:
        -fieldset
        -pset
        -param: needs key_alltracers (if True, T and NPP are sampled) and
        periodicBC (True for east/west periodicity)
    """
    key_alltracers = param['key_alltracers']
    periodicBC = param['periodicBC']
    #
    kernels_list = [pk.IncrementAge, 
                    pk.BeachTesting, 
                    pk.UndoMove, 
                    pk.Distance, 
                    pk.CurrentVelocity]
    #
    if key_alltracers:
        kernels_list.append(pk.SampleTracers)
    #
    if periodicBC:
        kernels_list.append(pk.Periodic)
    #
    for k in range(len(kernels_list)):
        kernels_list[k]=pset.Kernel(kernels_list[k])  
    return kernels_list



def define_active_kernels(pset, param):
    """
    Function that defines additional kernel that will be used for computation.
    Parameters:
        -pset: ParticleSet
        -param: needs mode (active or passive) and species (leatherback, loggerhead or green)
    """
    #
    mode = param['mode']
    species = param['species']
    growth = param['growth']
    #
    kernels_list = []
    if mode == 'active':      
        if growth == 'VGBF':
            compute_SCL = ak.compute_SCL_VGBF
            compute_PPmax = ak.compute_PPmax_VGBF   
        elif growth == 'Gompertz':
            compute_SCL = ak.compute_SCL_Gompertz
            compute_PPmax = ak.compute_PPmax_Gompertz
        #
        kernels_list = [compute_SCL,
                        ak.compute_Mass,
                        compute_PPmax,
                        ak.compute_vmax,
                        ak.compute_habitat,
                        ak.compute_swimming_direction,
                        ak.compute_swimming_velocity]
        if param['Tmin_Topt'] == 'variable':
            kernels_list.insert(2, ak.compute_Tmin_Topt) #Needs to be after compute_Mass
    if param['cold_death']:
        kernels_list.append(ak.cold_induced_mortality)
        
    for k in range(len(kernels_list)):
        kernels_list[k]=pset.Kernel(kernels_list[k]) 
    return kernels_list



def sum_kernels(k_adv, k_active, k_passive):
    """
    Sums all the kernels and returns the summed kernel.
    WARNING: The position is important.
    """
    if k_active != []:
        kernels = k_active[0]
        print_kernels = [k_active[0].funcname]
        #
        for k in k_active[1:]:
            kernels = kernels + k
            print_kernels.append(k.funcname)
        kernels += k_adv
        print_kernels.append(k_adv.funcname)
    #    
    else:   
        kernels = k_adv
        print_kernels = [k_adv.funcname]
    #
    for k in k_passive:
        kernels = kernels + k
        print_kernels.append(k.funcname)
    #
    print('****************************************************')
    print("These kernels will be used for computation: \n")
    for k in print_kernels:
        print(k, "\n")
    print('****************************************************')
    
    return kernels


# =============================================================================
# OUTPUT
# =============================================================================
def modify_output(OutputFile, t_init, param):
    """
    Modify output file so that variables names are the same as in STAMM 2.0 and
    turtles live exactly ndays_simu days.
    """
    dt = math.ceil(max(t_init) - min(t_init))
    #
    nc_i = netCDF4.Dataset(OutputFile, 'r')
    name_out = OutputFile.replace('.nc', '0.nc')
    nc_o = netCDF4.Dataset(name_out, 'w')
    #
    nsteps = nc_i.dimensions['obs'].size
    nturtles = nc_i.dimensions['traj'].size
    #
    nc_o.createDimension('nsteps', size = nsteps - dt)
    nc_o.createDimension('nturtles', size = nturtles)
    #
    for var_name in nc_i.variables:
        if var_name not in ['time','trajectory','z']:
            var = nc_i.variables[var_name]
            if dt > 0:
                values = np.transpose(np.squeeze(var))[:-dt, :]
            else:
                values = np.transpose(np.squeeze(var))[:]
            tmp = nc_o.createVariable(var_name, 'f', ('nsteps','nturtles'))
            tmp[:] = values
    #
    nc_o.renameVariable('lat','traj_lat')
    nc_o.renameVariable('lon','traj_lon')
    nc_o.renameVariable('age','traj_time')
    if param['key_alltracers']:
        nc_o.renameVariable('T','traj_temp')
        nc_o.renameVariable('NPP','traj_pp')
    if param['mode'] == 'active':
        nc_o.renameVariable('xgradh','xgrad')
        nc_o.renameVariable('ygradh','ygrad')
    init_t = nc_o.createVariable('init_t', 'f', ('nturtles'))
    init_t[:] = t_init
    nc_o.close()
    nc_i.close()
    #delete initial OutputFile 
    subprocess.run(["mv", "-f", name_out, OutputFile])
    print('\n')
    print('********************************************************************************')
    print("Wrote", OutputFile)
    print('********************************************************************************')
    print('\n')
    
  
        
        
       
