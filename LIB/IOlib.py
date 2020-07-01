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
             'vgpm':'False',
             'growth':'VGBF',
             'grid_phy':'',
             'lon_T':'',
             'lat_T':'',
             'SCL0':'',
             'tactic_factor':'1',
             'frenzy':'False'
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
    for key in ['periodicBC','key_alltracers', 'cold_death','vgpm','halo','frenzy']:
        if items[key] == '':
            print("\n WARNING: %s not found, set to False \n"%key)
        try:
            items[key] = items[key].__contains__('T')
        except ValueError:
            sys.exit("ERROR : %s must be boolean" %(key))

    #Species name to lower caracters
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
        for key in ['alpha','P0','grad_dx','tactic_factor']:        
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
        param_check = set(list({'species','P0', 'alpha', 'grad_dx', 'SCL0'}) + list(param_check))

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
    
    if param['mode'] != 'active' and param['frenzy'] == True:
        raise ValueError('Frenzy swimming is available only in active mode')

    

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
        time_periodic += 1
    #
    if date_end > last_date:
        raise ValueError("Simulation ends after the date of last data file available. Please check parameter time_periodic or set it to auto")
    #
    print('   Date of first file: ', date_start)
    print('   Date of last file:  ', date_end)    
    print('\n')
    return date_start, date_end, time_periodic



def forcing_list(f_dir, f_suffix, date_start, date_end, vgpm = False):
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
    if tmin.days < len(files): #pas propre...
        del(files[:tmin.days])
    del(files[(date_end-date_start).days+2:])
    #
    if vgpm:
        files = []
        for yr in list_years:
            files += sorted(glob(f_dir + '/*' + '.' + str(yr) + '*' + f_suffix))
        #remove useless files
        del(files[:tmin.days//8])
        del(files[((date_end-date_start).days)//8+4:])
        if (date_end - datetime(date_end.year, 1, 1)).days >= 360: #add first file of following year
            files.append(glob(f_dir + '/*' + '.' + str(date_end.year + 1) + '001' + f_suffix)[0])
    #
    if files == []:
        print('   Years do not appear in file names of '+f_dir+', considering first file is 01/01/%d. \n'%date_start.year)
        files = sorted(glob(f_dir + '/*' + f_suffix))
        del(files[:tmin.days])
        del(files[(date_end-date_start).days+2:])
    return files
