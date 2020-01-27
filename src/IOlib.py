#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Module description : All input/output functions to read :
- input parameters in namelist file
- initial positions
- external variables files (U,V,T...)
- grid parameters
and write output file.


"""

import numpy as np
import netCDF4 as nc
import csv
import sys
import time
import datetime


def read_namelist(filename):
    """
    Read the namelist file and return a dictionnary containing the name of the items (str), and its corresponding value (str)
    read_namelist(filename) -> items
    - filename : str, name of the namelist file
    - items : dictionnary

    For more information on the items, read 'IBM-tutorial' file 
    """

    items = {'init_file':'',
             'type':'x/y',
             'nturtles':'',
             'tstep':'86400',
             'nsteps_simu':'',
             'time_periodic':'False',
             'species':'',
             'mode':'active',
             'alpha':'',
             'vscale':'',
             'P0':'',
             'nsteps_max':'',
             'key_periodic':'',
             'overlap':'0',
             'key_jfold':'',
             'pivot':'T',
             'c_dir_zo':'',
             'c_prefix_zo':'',
             'ind0_zo':'',
             'indn_zo':'',
             'maxsize_zo':'',
             'c_suffix_zo':'.nc',
             'nc_var_zo':'',
             'nc_lon_zo':'',
             'nc_lat_zo':'',
             'c_dir_me':'',
             'c_prefix_me':'',
             'ind0_me':'',
             'indn_me':'',
             'maxsize_me':'',
             'c_suffix_me':'.nc',
             'nc_var_me':'',
             'c_dir_te':'',
             'c_prefix_te':'',
             'ind0_te':'',
             'indn_te':'',
             'maxsize_te':'',
             'c_suffix_te':'.nc',
             'nc_var_te':'',
             'c_dir_pp':'',
             'c_prefix_pp':'',
             'ind0_pp':'',
             'indn_pp':'',
             'maxsize_pp':'',
             'c_suffix_pp':'.nc',
             'nc_var_pp':'',
             'dir_mesh':'',
             'fn_mesh':'',
             'nc_var_xx_tt':'',
             'nc_var_xx_uu':'',
             'nc_var_yy_tt':'',
             'nc_var_yy_vv':'',
             'nc_var_e2u':'',
             'nc_var_e1v':'',
             'nc_var_e1t':'',
             'nc_var_e2t':'',
             'nc_var_tmask':'',
             'key_CenteredGrid':'',
             'key_alltracers':'True',
             'time_origin':'',
             'key_bounce':'False',
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
    for key in ['nturtles','nsteps_simu','overlap','ind0_zo','indn_zo','maxsize_zo','ind0_me','indn_me','maxsize_me']:
        try:
            items[key] = int(items[key])
        except ValueError:
            sys.exit("ERROR : %s must be integer" %(key))
 
        
    #convert floats
    for key in ['tstep']:
        try:
            items[key] = float(items[key])
        except ValueError:
            sys.exit("ERROR : %s must be float" %(key))

    #convert booleans
    for key in ['key_periodic','key_jfold','key_CenteredGrid','key_alltracers','key_bounce']:
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
    if items['mode'] not in ['passive','active','diffusion']:
        sys.exit("ERROR : mode must be 'passive', 'active' or 'diffusion'")
        
    
    #Check pivot
    if items['pivot']!='T' and items['pivot']!='F':
        sys.exit("ERROR : pivot must be 'T' or 'F'")

    #optional items
    if items['key_alltracers'] == True:
        for key in ['ind0_te','indn_te','maxsize_te','ind0_pp','indn_pp','maxsize_pp']:
            try:
                items[key] = int(items[key])
            except ValueError:
                sys.exit("ERROR : %s must be integer" %(key))

    for key in ['alpha','vscale','P0']:
        if items['mode']=='active':
            try:
                items[key] = float(items[key])
            except ValueError:
                sys.exit("ERROR : %s must be float" %(key))
        else:
            items[key] = 0.
    
     #time_periodic
    if items['time_periodic'] == 'False':
        items['time_periodic'] = False
    else:
        try:
            items['time_periodic'] = int(items['time_periodic'])
        except ValueError:
            sys.exit("ERROR: time_periodic must be integer or set to False")
            
   
    #time_origin
    try:
        items['time_origin'] = datetime.datetime.strptime(items['time_origin'], '%Y/%m/%d')
    except:
        print("\n WARNING: time_origin has to be in format YYYY/MM/DD, time won't be printed \n")
        items['time_origin'] = False
    
    #nsteps_max
    if items['key_alltracers'] == True:
        items['nsteps_max'] = min(items['indn_zo'], items['indn_me'], items['indn_te'], items['indn_pp'])
    else:
        items['nsteps_max'] = min(items['indn_zo'], items['indn_me'])
    
    
    return items



def check_param(param,output_file):
    """
    Checks if all needed arguments are present.
    """
    param_check = {'init_file', 'nturtles', 'nsteps_simu', 'species', 'mode', 'key_alltracers', 'overlap',\
                   'pivot', 'key_CenteredGrid', 'key_periodic', 'key_jfold',\
                   'c_dir_zo', 'c_prefix_zo', 'ind0_zo', 'indn_zo', 'maxsize_zo', 'c_suffix_zo', 'nc_var_zo',\
                   'c_dir_me', 'c_prefix_me', 'ind0_me', 'indn_me', 'maxsize_me', 'c_suffix_me', 'nc_var_me'}
    
    if param['key_alltracers'] == True:
        param_check = set(list({'c_dir_te', 'c_prefix_te', 'ind0_te', 'indn_te', 'maxsize_te', 'c_suffix_te', 'nc_var_te',\
                  'c_dir_pp', 'c_prefix_pp', 'ind0_pp', 'indn_pp', 'maxsize_pp', 'c_suffix_pp', 'nc_var_pp'}) + list(param_check))
    if param['mode'] == 'active':
        param_check = set(list({'P0', 'alpha', 'vscale'}) + list(param_check))
        
        
    for key in param.keys():
        value = param[key]
        if key in param_check and value == '':
            raise ValueError('Please give a value to %s' %key)
    
    if param['mode'] == 'active' and param['key_alltracers'] == False:
        raise ValueError('In active mode key_alltracers has to be True')
    


def fx_inv(lon,lat,lon_mat,lat_mat):
    """
    This function gets the x,y position on grid corresponding to longitude and latitude
        lon : f, longitude
        lat : f, latitude
        lon_mat : matrix containing the value of the longitude at each grid point
        lat_mat : matrix containing the value of the latitude at each grid point
    """

    lon_min = np.min(lon_mat)
    lon_max = np.max(lon_mat)
    lat_min = np.min(lat_mat)
    lat_max = np.max(lat_mat)

    # Compute position of nearest point on grid
    if ((lon> lon_max) or (lon< lon_min) or (lat>lat_max) or (lat<lat_min)):
        print("Longitude or latitude out of range, must be between ", np.min(lon_mat), " & ", np.max(lon_mat), " for longitude and ", np.min(lat_mat), " & ", np.max(lat_mat), " for latitude")
        return 
    else:
        distance = (lon-lon_mat)**2+(lat-lat_mat)**2
        inds = np.argmin(distance)     
        #argmin compute index of the flattened array (array.flat()), unravel_index gives the corresponding indexes of the 2D array
        i2,i1 = np.unravel_index(inds,distance.shape) 

    # Compute position of 4 neighbours on grid
        #   nw ---- ne  i2+1
        #   |        |
        #   |        |
        #   |        |
        #   sw ---- se  i2
        #   i1    i1+1

    if (lon>=lon_mat[i2,i1]):
        if (lat>=lat_mat[i2,i1]):
            sw = i2,i1
            nw = i2+1,i1
            se = i2,i1+1
            ne = i2+1,i1+1
        else:
            nw = i2,i1
            sw = i2-1,i1
            se = i2-1,i1+1
            ne = i2,i1+1
    else:
        if (lat>=lat_mat[i2,i1]):
            sw = i2,i1-1
            nw = i2+1,i1-1
            se = i2,i1
            ne = i2+1,i1
        else:
            nw = i2,i1-1
            sw = i2-1,i1-1
            se = i2-1,i1
            ne = i2,i1

    a = lon - lon_mat[sw]
    b = lat - lat_mat[sw]
    xsize = lon_mat[se] - lon_mat[sw]
    ysize = lat_mat[nw] - lat_mat[sw]
        #Weighted barycentre of the 4 neighbour points (cell can be trapezoidal in Orca)
    i = (xsize-a)/xsize*(ysize-b)/ysize*sw[1] + (xsize-a)/xsize*b/ysize*nw[1] + a/xsize*(ysize-b)/ysize*se[1] + a/xsize*b/ysize*ne[1]+1
    j = (xsize-a)/xsize*(ysize-b)/ysize*sw[0] + (xsize-a)/xsize*b/ysize*nw[0] + a/xsize*(ysize-b)/ysize*se[0] + a/xsize*b/ysize*ne[0]+1
    return i,j


def read_positions(param):
    """
    Function to read initial positions in file.
    - param : dictionnary containing all items of the namelist, output of read_namelist
    Returns: lon_init, lat_init, t_init initial positions and initial times relative to first data date.
    """

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
        print("ERROR - There are not enough initial positions for intended simulation. Add more lines to initial positions file or choose a lower nturtles.")
        sys.exit(1)

    x_init = np.zeros(nturtles,dtype='float32')
    y_init = np.zeros(nturtles,dtype='float32')
    t_init = np.zeros(nturtles,dtype='float32')

    
    print('****************************************************')
    print("Read initial positions in file: \n ", param['init_file'])
    print('****************************************************')
    init = open(param['init_file'],'r')
    x_init, y_init, t_init = np.loadtxt(init,usecols=(0,1,3),unpack=True)
    x_init, y_init, t_init = x_init[:nturtles], y_init[:nturtles], t_init[:nturtles] 
    init.close()
        


    #Check that number of days to be simulated does not exceed max number of input files
    ind_max = param['nsteps_simu']+np.max(t_init)
    if ind_max > param['nsteps_max']:
        if param['time_periodic'] == False:
            print('WARNING - There are not enough input files: Loop over time when last one is reached')
            print('WARNING - max(t_init) = %d ' % np.max(t_init))
            print('WARNING - nsteps_max  = %d ' % (param['nsteps_max']))
            print('WARNING - nsteps_simu = %d, should be less than %d ' % (param['nsteps_simu'],param['nsteps_max']-np.max(t_init)))
            param['time_periodic'] = param['nsteps_max']
        
        elif param['time_periodic'] > ind_max:
            print('WARNING - time_periodic is greater than the number of files: Loop over time when last one is reached')
            param['time_periodic'] = param['nsteps_max']
            
    time.sleep(3) #pause so that user can read warnings    


    # #Initial cell cannot be on land
    # i0 = np.int32(x_init) + 1
    # j0 = np.int32(y_init) + 1

    # position_on_mask = mask[j0,i0]
    
    # land = np.where(position_on_mask==0)[0]
    
    # if len(land) != 0:
    #     print("ERROR : found positions on land:")
    #     print(str(land))
    #     quit()

    return x_init, y_init, t_init