#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 18 10:44:31 2019

@author: tcandela
"""
# =============================================================================
# IMPORTS
# =============================================================================
import netCDF4 as nc
import numpy as np
import os
import sys
import IOlib as IO
from pathlib import Path
import datetime as dt
import glob

# =============================================================================
# FUNCTIONS
# =============================================================================
def read_nc(infile, dict_keys):
    """ Read data from a NetCDF file to python dictionnary format. """

    nc_file=nc.Dataset(infile, 'r')
    data = {}.fromkeys(dict_keys)
    for key in data.keys():
        data[key] = (np.array(nc_file.variables[key][:]))
    nc_file.close()
    return data



def get_name(infile):
    """
    Get infile name without path and .nc (return str)
    """
    name = Path(infile).stem.replace('.nc','')
    return name


def get_directory(infile):
    """
    Get directory of infile (return str).
    """
    directory = str(Path(infile).parents[0]) + '/'
    return directory
    

def age_to_date(traj_time, t_init, lat, lon) :
    """
    Order lat and lon in time instead of age. They can be any other array.
    Returns lon and lat, arrays of size ndays_simu + (t_init[-1] - t_init[0]).
    Returns date_mat, an array giving date since 01/01/ystart at each tstep (same for all turtles)
    nan values when no value
    """
    # Load parameters and initialize variables.

    nturtle = len(t_init)
    duration = np.max(traj_time)
    n_steps = np.shape(traj_time)[0]
    t_step = traj_time[1,0]-traj_time[0,0]
    start = t_init.min()
    end = t_init.max() + duration
    final_duration = end - start
    final_n_steps = int(final_duration / t_step)+2
    # Build universal date vector.
    # (calendar zero 01/01/ystart)
    date = traj_time[:,0] + t_init.min()
    after = np.arange(date[-1] + t_step, end, t_step)
    date = np.concatenate((date, after))
    date_mat = np.zeros((len(date),nturtle))
    new_lon = np.zeros((final_n_steps, nturtle))
    new_lat = new_lon.copy()
    # Add initial and final positions to fit with universal calendar.
    for turtle in range(nturtle) :
        date_mat[:,turtle] = date
        steps_before = int((t_init[turtle]-start)/t_step)
        steps_after = final_n_steps - (n_steps + steps_before)
        lat_before = np.zeros(steps_before) + float('nan') 
        lat_after = np.zeros(steps_after) + float('nan')
        lon_before = np.zeros(steps_before) + float('nan')
        lon_after = np.zeros(steps_after) + float('nan')
        new_lon[:,turtle] = np.concatenate((lon_before, lon[:,turtle], lon_after))
        new_lat[:,turtle] = np.concatenate((lat_before, lat[:,turtle], lat_after))
    return new_lat, new_lon, date_mat


def age_to_date2(duration, init_t,traj,var_age):
    var_date = np.ones([duration+int(np.max(init_t)-1),np.shape(traj)[0]])*np.nan
    k = 0
    for i in traj:
        var_date[int(np.ceil(init_t[k])):int(duration+np.ceil(init_t[k])-1),k] = var_age[1:,k]
        k = k+1
    return var_date


def time_interpolation(array0, array1, date0, date, date1):
    """
    Time interpolation at date 'date' of two arrays with dates 
    'date0' (before') and 'date1' (after).
    Returns the interpolated array.
    """
    w = (date-date0)/(date1-date0) #Weights
    array = (1 - w )  * array0 + w * array1
    return array

def interpolate_vgpm(current_date, param):
    """
    Returns NPP at date 'current_date' based on files in food_dir.
    Work for 8 days NPP.
    #
    current_date: datetime object
    param: namelist parameters
    """
    year0 = current_date.year
    year1 = year0
    days = (current_date - dt.datetime(year0,1,1)).days + 1
    dt0 = 8
    #
    list_days = np.arange(1,365,dt0)
    if days in list_days:
        name = str(year0) + str(("%03d") %days) + param['food_suffix']
        file = glob.glob(param['food_dir'] + '/*' + name)[0]
        ncfile = nc.Dataset(file)
        array = np.squeeze(ncfile.variables[param['food_var']])
        return array
    else:
        days0 = days//dt0 * dt0 + 1
        days1 = days0 + dt0
        if days1 > 365:
            days1 = 1
            year1 += 1
        name0 = str(year0) + str(("%03d") %days0) + param['food_suffix']
        name1 = str(year1) + str(("%03d") %days1) + param['food_suffix']
        file0 = glob.glob(param['food_dir'] + '/*' + name0)[0]
        file1 = glob.glob(param['food_dir'] + '/*' + name1)[0]
        ncfile0 = nc.Dataset(file0)
        ncfile1 = nc.Dataset(file1)
        array0 = np.squeeze(ncfile0.variables[param['food_var']])
        array1 = np.squeeze(ncfile1.variables[param['food_var']])
        date0 = dt.datetime(year0,1,1) + dt.timedelta(days=days0-1)
        date1 = dt.datetime(year0,1,1) + dt.timedelta(days=days1-1)
        array = time_interpolation(array0, array1, date0, current_date, date1)
        return array

def extract_x_y(ncfile,duration, turtles, alive):


    print('#### Debut extract xy')


    x  = np.array(ncfile.variables['traj_lon'])
    print('   => traj_lon loaded')
    y = np.array(ncfile.variables['traj_lat'])
    print('   => traj_lat loaded')
    x[np.where(x<200)]=x[np.where(x<200)]+360
    
    days = np.arange(duration)
    traj = np.arange(turtles)
    
    if alive == True:
        temp = np.squeeze(ncfile.variables['traj_temp'])
    else:
        temp = None
    
    t_init = np.array(ncfile.variables['init_t'])

    return x,y,temp,t_init

    print('#### Fin extract xy')
          
def find_traj_area(data,duration, turtles, correct_lon,area):


    days = np.arange(duration)
    traj = np.arange(turtles)

    x = data.lon[:,traj][days,:]
    y = data.lat[:,traj][days,:]

    if correct_lon == True :
        x = np.array([l+(((1-np.sign(l))/2)*360) for l in x])
        x[np.where(x<200)]=x[np.where(x<200)]+360

    ind_area=[]

    #Select tracks drifted into defined areas
    if area is 'all':
        area = np.arange(turtles)

    elif area =='NATL':
        for j in range (0,turtles):
            if np.shape(np.where((y[:,j] >20) ))[1]>0 : 
                ind_area.append(j)

        traj = np.array(ind_area)
    
    return traj
          
def classement(vect_class,val,wrap=False):
    nbclass = len(vect_class)-1
    X = np.concatenate((vect_class,np.array([val])))
    test = np.where(vect_class!=np.sort(X)[:-1])[0]
    if len(test)>0:
        ind = test[0]-1
    else :
        if wrap:
            ind = len(vect_class)
            ind %= (nbclass+1)
        else :
            ind = len(vect_class)-1
    return ind          

def find_nearest(navlat, navlon, LATid, LONid):
    distance = (navlat - LATid)**2 +(navlon-LONid)**2
    ind = np.where((distance == np.min(distance)))
    indlat = ind[0][0]
    indlon = ind[1][0]
    return indlat, indlon

class data:
    """
    class data contains all results from simulation
    contains also all routines for data analysis
    """

    def __init__(self,filename=str()):
        """
        called with the path to the netcdf output file as an argument
        load all data

        filename : string, path to netcdf file
        """

        # file opening
        if os.path.exists(filename):
            infile = nc.Dataset(filename)
        else :
            sys.exit(str(filename)+"\n\n\n file does not exist\n\n")

        # loading data
        self.lon=infile.variables['traj_lon'][:,:]
        self.lon=np.float32(self.lon)
        print('   => longitude loaded',self.lon.shape)
        #
        self.lat=infile.variables['traj_lat'][:,:]
        self.lat=np.float32(self.lat)
        print('   => latitude loaded',self.lat.shape)
        #

        self.traj_time = infile.variables['traj_time'][:,:]
        self.traj_time = np.int32(self.traj_time)
        print('   => traj_time loaded')
        #
        self.init_t=np.float32(infile.variables['init_t'][:])
        #Shape of data
        self.nb_output,self.nturtles=self.lon.shape
        infile.close()

        # assigning data
        self.filename=filename
        self.nb_output,self.nturtles=self.lat.shape
        self.nb_year=(self.nb_output-1)/365. # for daily outputs
        self.lat0=np.mean(self.lat[0,:])
        self.lon0=np.mean(self.lon[0,:])
        print('Zone de depart : ' + str(self.lat0) + ', ' + str(self.lon0))

        #Center longitudes
        self.lon = np.where(self.lon<=self.lon0-180,self.lon+360,self.lon)
        self.lon = np.where(self.lon>self.lon0+180,self.lon-360,self.lon)

        #Ariane days
        days = np.zeros(self.lon.shape,dtype='int32')
        for k in range(self.nturtles):
            days[:,k] = np.int32((self.traj_time[:,k]-self.traj_time[0,k])+self.init_t[k])
        self.days=days
        
        
def data_lists(param, end_day, t_init):
    """
    Return a list of list of data names with the first name being the first simulation day 
    and the last name being the last file needed for simulation.
    Passive mode: returns U_list and V_list
    Active mode: also returns T_list and Food_list
    Update time_periodic.
    """
    mode = param['mode']
    ndays_simu = param['ndays_simu']
    last_date_simu = IO.find_last_date(param)       
    date_start, date_end, param['time_periodic'] = IO.define_start_end(ndays_simu, param, t_init, last_date_simu)
    #
    data_list = []
    U_list = IO.forcing_list(param['U_dir'], param['U_suffix'], date_start, date_end)
    V_list = IO.forcing_list(param['V_dir'], param['V_suffix'], date_start, date_end)
    data_list.append(U_list)
    data_list.append(V_list)
    #
    if mode == 'passive':    
        return data_list
    elif mode == 'active':
        T_list = IO.forcing_list(param['T_dir'], param['T_suffix'], date_start, date_end)
        Food_list = IO.forcing_list(param['food_dir'], param['food_suffix'], date_start, date_end)
        data_list.append(T_list)
        data_list.append(Food_list)
    #
    return data_list
    
def classify_lon_init(dico, nb_cat):
    """
    Give a group number to each turtle depending on their release longitude (used for jeanette project)
    """
    init_lon = dico['traj_lon'][0, :]
    
    nb_turtles = len(init_lon)

    init_lon_min = np.min(init_lon)  
    init_lon_max = np.max(init_lon)

    sort_lon = np.linspace(init_lon_min, init_lon_max + 0.001, nb_cat + 1)

    group = np.zeros(nb_turtles)
    for turtle in np.arange(nb_turtles):
        sort = False
        for cat in np.arange(len(sort_lon) - 1):
            if (init_lon[turtle] >= sort_lon[cat]) and (init_lon[turtle] < sort_lon[cat + 1]) and (sort == False):
                group[turtle] = cat
                sort = True
    return group
 
    
    
    
    
    
