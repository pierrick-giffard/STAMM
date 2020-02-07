#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 31 15:56:33 2020

@author: pgiffard
"""
from glob import glob
from parcels import FieldSet, ParticleSet
from datetime import timedelta, date
import numpy as np
import time

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
    ufiles = forcing_list(param['U_dir'], param['U_suffix'], param['ystart'], param['ndays_simu'], t_init, print_date=True)
    vfiles = forcing_list(param['V_dir'], param['V_suffix'], param['ystart'], param['ndays_simu'], t_init)
    #Filenames
    filenames = {'U': {'lon': mesh_phy, 'lat': mesh_phy, 'data': ufiles},
                 'V': {'lon': mesh_phy, 'lat': mesh_phy, 'data': vfiles}}
    #Variables
    variables = {'U': param['U_var'],
                 'V': param['V_var']}
    #Dimensions: C-grids need f nodes
    dimensions = {'U': {'lon': param['lon_phy'], 'lat': param['lat_phy'], 'time': param['time_var_phy']},
                  'V': {'lon': param['lon_phy'], 'lat': param['lat_phy'], 'time': param['time_var_phy']}} 
    if key_alltracers:
        tfiles = forcing_list(param['T_dir'], param['T_suffix'], param['ystart'], param['ndays_simu'], t_init)
        ffiles = forcing_list(param['food_dir'], param['food_suffix'], param['ystart'], param['ndays_simu'], t_init)
        #
        filenames['T'] = {'lon': mesh_phy, 'lat': mesh_phy, 'data': tfiles}
        filenames['NPP'] = {'lon': mesh_food, 'lat': mesh_food, 'data': ffiles}
        #
        variables['T'] = param['T_var']
        variables['NPP'] = param['food_var']
        #
        dimensions['T'] = {'lon': param['lon_phy'], 'lat': param['lat_phy'], 'time': param['time_var_phy']}
        dimensions['NPP'] = {'lon': param['lon_food'], 'lat': param['lat_food'], 'time': param['time_var_food']}
 
    #Fieldset
    if param['time_periodic']:
        fieldset = FieldSet.from_netcdf(filenames, variables, dimensions,\
                                    time_periodic=timedelta(days=param['time_periodic']),\
                                    field_chunksize=(400, 400)) 
    else:
        fieldset = FieldSet.from_netcdf(filenames, variables, dimensions, field_chunksize=(400, 400)) 
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


def forcing_list(f_dir, f_suffix, ystart, ndays_simu, t_init, print_date = False):
    """
    Return a list with data needed for simulation.
    """
    date_start = date(ystart, 1, 1) + timedelta(days=1+int(np.min(t_init)))
    date_end = date_start + timedelta(days=ndays_simu+1)
    #
    list_years = np.arange(ystart, date_end.year + 1)
    files = []
    for yr in list_years:
        files += sorted(glob(f_dir + '/*' + str(yr) + '*' + f_suffix))
    #remove useless files
    del(files[:int(np.min(t_init))+1])
    del(files[ndays_simu+2:])    
    #
    if files == []:
        print('   Years not in file names, loading the whole dataset.')
        files = sorted(glob(f_dir + '/*' + f_suffix))
        del(files[:int(np.min(t_init))+1])
    #
    if print_date:
        print('   Starting day: ', date_start)
        print('   Last day:     ', date_end)
        print('   First U file used: \n  ', files[0])
        print('   Last U file used: \n  ', files[-1])
        print('\n')
    return files


def create_particleset(fieldset, pclass, lon, lat, t_init):
    print('\n')
    print('****************************************************')
    print("Building ParticleSet...")
    print('****************************************************')
    #
    t0 = time.time()
    #
    pset = ParticleSet(fieldset, pclass=pclass, lon=lon, lat=lat, time = (t_init - int(np.min(t_init))) * 86400)
    #Time
    tt=time.time()-t0
    print('\n')
    print(' => ParticleSet created in: '+ str(timedelta(seconds=int(tt))))
    print('\n')
    return pset


def initialization(pset, param):
    """
    Links constant parameters to particle in order to use them within kernels.
    """
    for p in pset:  
        p.vscale = param['vscale']
        p.P0 = param['P0']
        p.grad_dx = param['grad_dx']
        p.alpha = param['alpha']
        if param['mode'] == 'active':
            p.mode = 1
        elif param['mode'] == 'passive':
            p.mode = 0
        p.tstep = param['tstep']

    