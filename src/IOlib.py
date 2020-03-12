#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Functions to read namelist and initial positions
"""

import numpy as np
import csv
import sys


def read_namelist(filename):
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
             'time_periodic':'False',
             't_output':'',
             'adv_scheme':'RK4',
             'species':'',
             'mode':'active',
             'alpha':'',
             'vscale':'',
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
             'SCL0':'',
             'K0':''
             }

    namelist = open(filename,'r')
    text = csv.reader(namelist,delimiter='=')
    
    print("    ************")
    print("    * NAMELIST *")
    print("    ************")
    for line in text:
        print(str(line).replace('[','').replace(']','').replace('\'','').replace(',',':'))
        for key in items.keys():
            if line[0].__contains__(key):
                items[key]=line[1].replace('\'','').replace(' ','').replace(',','')
                
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
    for key in ['periodicBC','key_alltracers', 'cold_death','vgpm','halo']:
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
    
   
    
    #Active items
    if items['mode']=='active':
        for key in ['alpha','vscale','P0','grad_dx']:        
            try:
                items[key] = float(items[key])
            except ValueError:
                sys.exit("ERROR : %s must be float" %(key))
                            
    
    #Time periodic
    if items['time_periodic'] == 'False':
        items['time_periodic'] = False
    else:
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
                   'U_suffix', 'V_suffix','ystart'}
    
    if param['key_alltracers'] == True:
        param_check = set(list({'T_dir','food_dir', 'T_var', 'food_var','T_suffix', 'food_suffix',\
                                'mesh_food', 'lon_food', 'lat_food', 'time_var_food'}) + list(param_check))
    if param['mode'] == 'active':
        param_check = set(list({'species','P0', 'alpha', 'vscale', 'grad_dx'}) + list(param_check))
          
    if param['growth'] == 'Gompertz':
        param_check = set(list({'SCL0','K0'}) + list(param_check))
                   
    for key in param.keys():
        value = param[key]
        if key in param_check and value == '':
            raise ValueError('Please give a value to %s' %key)
    
    if param['mode'] == 'active' and param['key_alltracers'] == False:
        raise ValueError('In active mode key_alltracers has to be True')
    
    if param['cold_death'] and param['key_alltracers'] == False:
        raise ValueError('To compute cold_death key_alltracers has to be True')    
    

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