#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
In/out library.
Functions to read namelist, initial positions
"""

import numpy as np
import csv
import sys
from glob import glob
from datetime import datetime, timedelta
import netCDF4
import xarray as xr
import pandas as pd

def read_namelist(filename, display=True):
    """
    Read the namelist file and return a dictionnary containing the name of the items (str), and its corresponding value (str)
    read_namelist(filename) -> items
    - filename : str, name of the namelist file
    - items : dictionnary
    """

    items = {'init_file':'',
             'nturtles':'',
             'tstep':'',
             'ndays_simu':'',
             'time_periodic':'auto',
             't_output':'',
             'adv_scheme':'RK4',
             'species':'',
             'mode':'active',
             'alpha':'',
             'P0':'',
             'grad_dx':'',
             'periodicBC':'True',
             'key_alltracers':'True',
             'halo':'False',
             'U_dir':'',
             'U_suffix':'',
             'V_dir':'',
             'V_suffix':'',
             'mesh_phy':'',
             'lon_phy':'',
             'lat_phy':'',
             'time_var_phy':'',
             'U_var':'',
             'V_var':'',
             'mesh_food':'',
             'lon_food':'',
             'lat_food':'',
             'time_var_food':'',
             'T_dir':'',
             'T_suffix':'',
             'food_dir':'',
             'food_suffix':'',
             'T_var':'',
             'food_var':'',
             'ystart':'',
             'cold_death':'False',
             'cold_resistance':'30',
             'growth':'VGBF',
             'grid_phy':'',
             'lon_T':'',
             'lat_T':'',
             'SCL0':'',
             'tactic_factor':'1',
             'frenzy':'False',
             'wave_swim':'False',
             'wave_dir':'',
             'wave_suffix':'',
             'Ust_var':'',
             'Vst_var':'',
             'time_extrapolation':'False',
             'dt_swim':''
             }

    namelist = open(filename,'r')
    text = csv.reader(namelist,delimiter='=')
    if display:
        print("    ************")
        print("    * NAMELIST *")
        print("    ************")
    for line in text:
        if display:
            print(str(line).replace('[','').replace(']','').replace('\'','').replace(',',':'))
        for key in items.keys():
            if line[0].__contains__(key):
                items[key]=line[1].replace('\'','').replace(' ','').replace(',','')
    if display:           
        print("    ****************")
        print("    * END NAMELIST *")
        print("    ****************")
    
    """
    All items are read as strings, they must be converted to correct type
    """
    #Convert integers
    for key in ['nturtles','ndays_simu','t_output','tstep','ystart','cold_resistance']:
        try:
            items[key] = int(items[key])
        except ValueError:
            sys.exit("ERROR : %s must be integer" %(key))
   
    #Convert booleans
    for key in ['periodicBC','key_alltracers', 'cold_death','halo','frenzy','wave_swim','time_extrapolation']:
        if items[key] == '':
            print("\n WARNING: %s not found, set to False \n"%key)
        try:
            items[key] = items[key].__contains__('T')
        except ValueError:
            sys.exit("ERROR : %s must be boolean" %(key))
    
    # Species name to lower caracters
    items['species'] = items['species'].lower()

    #Check mode
    items['mode'] = items['mode'].lower()
    if items['mode'] not in ['passive','active']:
        sys.exit("ERROR : mode must be 'passive' or 'active'")
        
    if items['growth'] not in ['VGBF','Gompertz']:
        sys.exit("ERROR : mode must be 'Gompertz' or 'VGBF'")
        
    if items['adv_scheme'] not in ['RK4','Euler']:
        sys.exit("ERROR : mode must be 'RK4' or 'Euler'")
    
    #SCL
    try:
        items['SCL0'] = float(items['SCL0'])
    except ValueError:
        sys.exit("ERROR : %s must be float" %('SCL0'))
   
    
    #Active items
    if items['mode']=='active':
        for key in ['alpha','P0','grad_dx','tactic_factor','dt_swim']:        
            try:
                items[key] = float(items[key])
            except ValueError:
                sys.exit("ERROR : %s must be float" %(key))
                            
    
    #Time periodic
    if items['time_periodic'] == 'False':
        items['time_periodic'] = False
    elif items['time_periodic'] != 'auto':
        try:
            items['time_periodic'] = int(items['time_periodic'])
        except ValueError:
            sys.exit("ERROR: time_periodic must be integer or set to False")
    

    return items



def check_param(param,output_file):
    """
    Checks if all needed arguments are present.
    """
    param_check = {'init_file', 'nturtles', 'ndays_simu', 't_output', 'mode', 'key_alltracers',\
                   'periodicBC', 'tstep', 'adv_scheme','halo',\
                   'U_dir', 'V_dir', 'mesh_phy', 'lon_phy', 'lat_phy', 'time_var_phy', 'U_var','V_var',
                   'U_suffix', 'V_suffix','ystart', 'grid_phy'}
    
    if param['key_alltracers'] == True:
        param_check = set(list({'T_dir','food_dir', 'T_var', 'food_var','T_suffix', 'food_suffix',\
                                'mesh_food', 'lon_food', 'lat_food', 'time_var_food'}) + list(param_check))
    if param['mode'] == 'active':
        param_check = set(list({'species','P0', 'alpha', 'grad_dx', 'SCL0','dt_swim'}) + list(param_check))

    if param['grid_phy'] == 'C':
        param_check = set(list({'lat_T', 'lon_T'}) + list(param_check))
                   
    for key in param.keys():
        value = param[key]
        if key in param_check and value == '':
            raise ValueError('Please give a value to %s' %key)
    
    if param['mode'] == 'active' and param['key_alltracers'] == False:
        raise ValueError('In active mode key_alltracers has to be True')
    
    if param['cold_death'] and param['key_alltracers'] == False:
        raise ValueError('To compute cold_death key_alltracers has to be True')
    
    if param['grid_phy'] != 'A' and param['grid_phy'] != 'C':
        raise ValueError("Set grid_phy to A or to C")
    
    if param['frenzy'] and param['wave_swim']:
        raise ValueError('Choose between swimming frenzy and swimming against waves')
        
    if type(param['time_periodic']) == int and param['time_extrapolation']:
        raise ValueError('Choose between time_extrapolation and time_periodic')
    

def read_positions(param):
    """
    Function to read initial positions in file.
    - param : dictionnary containing all items of the namelist, output of read_namelist
    Returns: lon_init, lat_init, t_init initial positions and initial times relative to first data date.
    """
    print('****************************************************')
    print("Read initial positions in file", param['init_file'])
    print('****************************************************')
    nturtles = param['nturtles']

    #Count number of lines and check if there are enough positions 
    try:
        init = open(param['init_file'],'r')
    except IOError:
        sys.exit("Initial positions file does not exist")
    lines = 0
    for line in init:
        lines += 1
    init.close()

    if nturtles > lines:
        print("   ERROR - There are not enough initial positions for intended simulation. Add more lines to initial positions file or choose a lower nturtles.")
        sys.exit(1)
    elif nturtles < lines:
        print('   WARNING - There are more initial positions than turtles. \n')

    x_init = np.zeros(nturtles,dtype='float32')
    y_init = np.zeros(nturtles,dtype='float32')
    t_init = np.zeros(nturtles,dtype='float32')


    init = open(param['init_file'],'r')
    x_init, y_init, t_init = np.loadtxt(init,usecols=(0,1,3),unpack=True)
    x_init, y_init, t_init = x_init[:nturtles], y_init[:nturtles], t_init[:nturtles] 
    init.close()
        

    return x_init, y_init, t_init



def find_last_date(param):
    """
    Returns the last date it is possible to use for computation (datetime).
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
    t_value = int(file_U.variables[time_var_phy][-1].data)
    time_U = netCDF4.num2date(t_value, t_unit)
    file_U.close()
    print(last_U)
    #
    last_V = sorted(glob(V_dir + '/*' + V_suffix))[-1]
    file_V = netCDF4.Dataset(last_V)
    t_unit = file_V.variables[time_var_phy].units
    t_value = int(file_V.variables[time_var_phy][-1].data)
    time_V = netCDF4.num2date(t_value, t_unit)
    file_V.close()
    #
    if param['key_alltracers']:
        last_T = sorted(glob(T_dir + '/*' + T_suffix))[-1]
        file_T = netCDF4.Dataset(last_T)
        t_unit = file_T.variables[time_var_phy].units
        t_value = int(file_T.variables[time_var_phy][-1].data)
        time_T = netCDF4.num2date(t_value, t_unit)
        file_T.close()
        #
        last_food = sorted(glob(food_dir + '/*' + food_suffix))[-1]
        file_food = netCDF4.Dataset(last_food)
        t_unit = file_food.variables[time_var_food].units
        t_value = int(file_food.variables[time_var_food][-1].data)
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
    time_extra = param['time_extrapolation']
    #
    date_start = datetime(ystart, 1, 1) + timedelta(days=np.min(t_init))
    
    #
    if type(time_periodic) ==  int and time_periodic > ndays_simu:
        time_periodic = False
    #
    if time_periodic == False:
        date_end = date_start + timedelta(days=ndays_simu)
    #
    elif time_periodic == 'auto':
        if time_extra:
            date_end = date_start + timedelta(days=ndays_simu)
            time_periodic = False
        elif date_start + timedelta(days=ndays_simu) > last_date:
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
        time_periodic += 1
    
    if date_end > last_date and not time_extra:
        raise ValueError("Simulation ends after the date of last data file available. Please check parameter time_periodic or set it to auto. \n \
                          last_date: ", last_date, "date_end: ", date_end)
    
    print('   Date of first file: ', date_start)
    print('   Date of last file:  ', date_end)    
    print('\n')
    return date_start, date_end, time_periodic



def forcing_list(f_dir, f_suffix, date_start, date_end):
    """
    Return a list with data needed for simulation.
    It is important that the first file is the first day of release.
    Need file to have a time variable.
    """
    
    files = sorted(glob(f_dir + '/*' + f_suffix))
    t0 = pd.to_datetime(xr.open_dataset(files[0]).time.data[0])
    t1 = pd.to_datetime(xr.open_dataset(files[1]).time.data[0])
    dt = (t1 - t0).days
    
    i0 = (date_start - t0).days // dt
    i1 = (date_end - t0).days // dt + 1
    
    if dt == 8: #consider it is vgpm data, last file beeing 361 day of the year = pb each year
        nyear0 = date_start.year - t0.year
        nyear1 = date_end.year - t0.year
        i0 += nyear0 - (5 * nyear0) // dt # 5 days from 361 to 001
        i1 += nyear1 - (5 * nyear1) // dt-3 # 5 days from 361 to 001

        
    
    
    files = files[i0:i1+1]
    print(i0, i1, files[0], files[-1])
    
    return files

