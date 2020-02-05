#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 31 15:56:33 2020

@author: pgiffard
"""
from glob import glob
from parcels import FieldSet

def build_fieldset(param):
    """
    Build fieldset from data files, dimensions and variables.
    """
    print('****************************************************')
    print("Checking data...")
    print('****************************************************')
    #
    if param['key_alltracers']:
        #Mesh masks
        mesh_phy = param['mesh_phy']
        mesh_food = param['mesh_food']
        if mesh_food == 'mesh_phy':
            mesh_food = mesh_phy
            param['lon_food'] = param['lon_phy']
            param['lat_food'] = param['lat_phy']
        #Forcings
        ufiles = sorted(glob(param['U_files']))
        vfiles = sorted(glob(param['V_files']))
        tfiles = sorted(glob(param['T_files']))
        ffiles = sorted(glob(param['food_files']))
        #Filenames
        filenames = {'U': {'lon': mesh_phy, 'lat': mesh_phy, 'data': ufiles},
                     'V': {'lon': mesh_phy, 'lat': mesh_phy, 'data': vfiles},
                     'T': {'lon': mesh_phy, 'lat': mesh_phy, 'data': tfiles},
                     'NPP': {'lon': mesh_food, 'lat': mesh_food, 'data': ffiles}}
        #Variables
        variables = {'U': param['U_var'],
                     'V': param['V_var'],
                     'T': param['T_var'],
                     'NPP': param['food_var']}
        #Dimensions: need f nodes
        dimensions = {'U': {'lon': param['lon_phy'], 'lat': param['lat_phy'], 'time': param['time_var_phy']},
                      'V': {'lon': param['lon_phy'], 'lat': param['lat_phy'], 'time': param['time_var_phy']},
                      'T': {'lon': param['lon_phy'], 'lat': param['lat_phy'], 'time': param['time_var_phy']},
                      'NPP': {'lon': param['lon_food'], 'lat': param['lat_food'], 'time': param['time_var_food']}} 
    
    else:
        #Mesh masks
        mesh_phy = param['mesh_phy']
        #Forcings
        ufiles = sorted(glob(param['U_files']))
        vfiles = sorted(glob(param['V_files']))
        #Filenames
        filenames = {'U': {'lon': mesh_phy, 'lat': mesh_phy, 'data': ufiles},
                     'V': {'lon': mesh_phy, 'lat': mesh_phy, 'data': vfiles}}
        #Variables
        variables = {'U': param['U_var'],
                     'V': param['V_var']}
        #Dimensions: need f nodes
        dimensions = {'U': {'lon': param['lon_phy'], 'lat': param['lat_phy'], 'time': param['time_var_phy']},
                      'V': {'lon': param['lon_phy'], 'lat': param['lat_phy'], 'time': param['time_var_phy']}} 
    #Fieldset
    print('****************************************************')
    print("Building FieldSet and ParticleSet...")
    print('****************************************************')
    fieldset = FieldSet.from_netcdf(filenames, variables, dimensions)  #time_periodic=delta(days=365), allow_time_extrapolation=True
    
    
    
    return fieldset



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

    