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
from glob import glob
from parcels import FieldSet, ParticleSet, AdvectionRK4, AdvectionEE
from datetime import timedelta, date, datetime
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
    last_date = find_last_date(param)
    date_start, date_end, time_periodic = define_start_end(ndays_simu, param, t_init, last_date)
    ufiles = forcing_list(param['U_dir'], param['U_suffix'], date_start, date_end, print_date=True)
    vfiles = forcing_list(param['V_dir'], param['V_suffix'], date_start, date_end)
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
        tfiles = forcing_list(param['T_dir'], param['T_suffix'], date_start, date_end)
        ffiles = forcing_list(param['food_dir'], param['food_suffix'], date_start, date_end, vgpm=param['vgpm'])
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


def find_last_date(param):
    """
    Returns the last date it is possible to use for computation.
    It is the smaller last date of all data files. 
    """
    U_dir = param['U_dir']
    U_suffix = param['U_suffix']
    V_dir = param['V_dir']
    V_suffix = param['V_suffix']
    time_var_phy = param['time_var_phy']
    if param['key_alltracers']:
        T_dir = param['T_dir']
        T_suffix = param['T_suffix']
        food_dir = param['food_dir']
        food_suffix = param['food_suffix']
        time_var_food = param['time_var_food']
    #
    last_U = sorted(glob(U_dir + '/*' + U_suffix))[-1]
    file_U = netCDF4.Dataset(last_U)
    t_unit = file_U.variables[time_var_phy].units
    t_value = int(file_U.variables[time_var_phy][:].data)
    time_U = netCDF4.num2date(t_value, t_unit)
    file_U.close()
    #
    last_V = sorted(glob(V_dir + '/*' + V_suffix))[-1]
    file_V = netCDF4.Dataset(last_V)
    t_unit = file_V.variables[time_var_phy].units
    t_value = int(file_V.variables[time_var_phy][:].data)
    time_V = netCDF4.num2date(t_value, t_unit)
    file_V.close()
    #
    if param['key_alltracers']:
        last_T = sorted(glob(T_dir + '/*' + T_suffix))[-1]
        file_T = netCDF4.Dataset(last_T)
        t_unit = file_T.variables[time_var_phy].units
        t_value = int(file_T.variables[time_var_phy][:].data)
        time_T = netCDF4.num2date(t_value, t_unit)
        file_T.close()
        #
        last_food = sorted(glob(food_dir + '/*' + food_suffix))[-1]
        file_food = netCDF4.Dataset(last_food)
        t_unit = file_food.variables[time_var_food].units
        t_value = int(file_food.variables[time_var_food][:].data)
        time_food = netCDF4.num2date(t_value, t_unit)
        file_food.close()
        #
        last_file = min(time_U, time_V, time_T, time_food)
    else:
        last_file = min(time_U, time_V)
    return last_file


def define_start_end(ndays_simu, param, t_init, last_date):
    """
    Returns date of first simulation day, date of last data file to use and update time_periodic.
    """
    #
    ystart = param['ystart']
    time_periodic = param['time_periodic']
    #
    tmin = timedelta(days=int(np.min(t_init)))
    date_start = datetime(ystart, 1, 1) + tmin
    #
    if isinstance(time_periodic, int) and time_periodic > ndays_simu:
        time_periodic = False
    #
    if time_periodic == False:
        date_end = date_start + timedelta(days=ndays_simu)
    #
    elif time_periodic == 'auto':
        if date_start + timedelta(days=ndays_simu) > last_date:
            if last_date.month == 12 and last_date.day == 31:
                last_year = last_date.year
            else:
                last_year = last_date.year - 1
            date_end = datetime(last_year, 12, 31)
            time_periodic = (date_end - date_start).days
            print('time_periodic is set to %d'%time_periodic)
            if date_end < date_start:
                raise ValueError('Not enough data files. You can try to set time_periodic manually')
        else:
           time_periodic = False
           date_end = date_start + timedelta(days=ndays_simu)

    #if time_periodic is integer
    else:
        date_end = date_start + timedelta(days=time_periodic)
    #
    if date_end > last_date:
        raise ValueError("Simulation ends after the date of last data file available. Please check parameter time_periodic or set it to auto")
    #
    print('   Date of first file: ', date_start)
    print('   Date of last file:  ', date_end)    
    print('\n')
    return date_start, date_end, time_periodic

def forcing_list(f_dir, f_suffix, date_start, date_end, print_date = False, vgpm = False):
    """
    Return a list with data needed for simulation.
    It is important that the first file is the first day of release.
    This function highly depends on files names, it might not work for particular names.
    It works for names format:
        -   *_YYYY*suffix
        -   *YYYY*suffix
        -   vgpm files with vgpm=True in namelist
    If none of this format is found, consider first file is 01/01/ystart.
    """
    tmin = date_start - datetime(date_start.year, 1, 1)
    #
    list_years = np.arange(date_start.year, date_end.year + 1)
    files = []
    for yr in list_years:
        files += sorted(glob(f_dir + '/*' + '_' + str(yr) + '*' + f_suffix))
    if files == []:
        for yr in list_years:
            files += sorted(glob(f_dir + '/*' + str(yr) + '*' + f_suffix))
    #remove useless files
    del(files[:tmin.days])
    del(files[(date_end-date_start).days+2:])
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
        print('   Years do not appear in file names of '+f_dir+', considering first file is 01/01/%d. \n'%date_start.year)
        files = sorted(glob(f_dir + '/*' + f_suffix))
        del(files[:tmin.days])
        del(files[(date_end-date_start).days+2:])
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
    
  
        
        
       
