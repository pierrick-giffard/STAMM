#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Functions needed to:
    -create FieldSet
    -create ParticleSet
    -define kernels to be used    
"""

# =============================================================================
# IMPORTS
# =============================================================================
#Python libraries
from glob import glob
from parcels import FieldSet, ParticleSet, AdvectionRK4, AdvectionEE
from datetime import timedelta, date
import numpy as np
import time

#Personal libraries
import Advection_kernel as adv
import Passive_kernels as pk
import Active_kernels as ak
import Leatherback as leath
import Loggerhead as log



# =============================================================================
# FIELDSET
# =============================================================================
def create_fieldset(param, t_init):
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
    ufiles = forcing_list(param['U_dir'], param['U_suffix'], param['ystart'], param['ndays_simu'], t_init, param['time_periodic'], print_date=True)
    vfiles = forcing_list(param['V_dir'], param['V_suffix'], param['ystart'], param['ndays_simu'], t_init, param['time_periodic'])
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
        tfiles = forcing_list(param['T_dir'], param['T_suffix'], param['ystart'], param['ndays_simu'], t_init, param['time_periodic'])
        ffiles = forcing_list(param['food_dir'], param['food_suffix'], param['ystart'], param['ndays_simu'], t_init, param['time_periodic'], vgpm=param['vgpm'])
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
    if param['time_periodic']:
        fieldset = FieldSet.from_netcdf(filenames, variables, dimensions, time_periodic=timedelta(days=param['time_periodic']))
    else:
        fieldset = FieldSet.from_netcdf(filenames, variables, dimensions)#, field_chunksize=(800, 800))
    #East/West periodicity
    if param['periodicBC']:
        if param['grid_type'] == 'standard':
            #halo should be automatically included in orca grids 
            try:
                fieldset.add_constant('halo_west', fieldset.U.grid.lon[0,0])
                fieldset.add_constant('halo_east', fieldset.U.grid.lon[0,-1])
            except:
                fieldset.add_constant('halo_west', fieldset.U.grid.lon[0])
                fieldset.add_constant('halo_east', fieldset.U.grid.lon[-1])
            #          
            fieldset.add_periodic_halo(zonal=True) 
    #Time
    tt=time.time()-t0
    print('\n')
    print(' => FieldSet created in: '+ str(timedelta(seconds=int(tt))))
    print('\n')
    return fieldset


def forcing_list(f_dir, f_suffix, ystart, ndays_simu, t_init, time_periodic, print_date = False, vgpm = False):
    """
    Return a list with data needed for simulation.
    This function highly depends on files names, it might not work for particular names.
    It works for names format:
        -   *_YYYY*suffix
        -   *YYYY*suffix
        -   vgpm files with vgpm=True in namelist
    If none of this format is found, consider first file is 01/01/ystart.
    """
    tmin = timedelta(days=int(np.min(t_init)))
    date_start = date(ystart, 1, 1) + tmin
    #
    if time_periodic == False or time_periodic > (ndays_simu + np.max(t_init)):
        date_end = date_start + timedelta(days=ndays_simu+1)
    else:
        date_end = date(ystart, 1, 1) + timedelta(days=time_periodic+1)
    #
    list_years = np.arange(ystart, date_end.year + 1)
    files = []
    for yr in list_years:
        files += sorted(glob(f_dir + '/*' + '_' + str(yr) + '*' + f_suffix))
    if files == []:
        for yr in list_years:
            files += sorted(glob(f_dir + '/*' + str(yr) + '*' + f_suffix))
    #remove useless files
    del(files[:tmin.days])
    del(files[(date_end-date_start).days+1:])
    #
    if vgpm:
        files = []
        for yr in list_years:
            files += sorted(glob(f_dir + '/*' + '.' + str(yr) + '*' + f_suffix))
        #remove useless files
        del(files[:tmin.days//8])
        del(files[((date_end-date_start).days)//8+3:]) 
    #
    if files == []:
        print('   Years do not appear in file names of '+f_dir+', considering first file is 01/01/%d. \n'%ystart)
        files = sorted(glob(f_dir + '/*' + f_suffix))
        del(files[:tmin.days])
        del(files[(date_end-date_start).days+1:]) 
    #
    if print_date:
        print('   Date first file: ', date_start)
        print('   Date last file:  ', date_end)
        if date_end < date_start:
            raise ValueError("Date of last file is lower than date of first file: \
                             please note that time_periodic is defined respect to 01/01/ystart")
        print('\n')
    return files



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
    
        if param['key_alltracers']:
            fieldset.T.grid.lon[:,720] = 0
            fieldset.T.grid.lat[320,:] = 0
            fieldset.NPP.grid.lon[:,720] = 0
            fieldset.NPP.grid.lat[320,:] = 0
    except:
        print('Not Using PSY patch')




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


def initialization(fieldset, param):
    """
    Links constant parameters to fieldset in order to use them within kernels.
    """
    fieldset.tstep = param['tstep']
    fieldset.deg = 111195 #1degree = 111,195 km approx
    fieldset.cold_resistance = param['cold_resistance'] * 86400 #from days to seconds
    if param['mode'] == 'active':
        fieldset.active = 1
        fieldset.vscale = param['vscale']
        fieldset.P0 = param['P0']
        fieldset.grad_dx = param['grad_dx']
        fieldset.alpha = param['alpha']
    else:
        fieldset.active = 0



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
    kernels_list = [pk.IncrementAge, pk.BeachTesting, pk.UndoMove, pk.Distance, pk.CurrentVelocity]
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
    #
    kernels_list = [] 
    if mode == 'active':      
        if species == 'leatherback':
            file = leath
        elif species == 'loggerhead':
            file = log

        kernels_list.append(file.compute_SCL)      
        kernels_list.append(file.compute_Mass)
        kernels_list.append(file.compute_PPmax)
        kernels_list.append(file.compute_vmax)      
        kernels_list.append(file.compute_habitat)
        kernels_list.append(ak.compute_swimming_direction)
        kernels_list.append(ak.compute_swimming_velocity)
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