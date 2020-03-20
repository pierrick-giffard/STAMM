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




def age_to_date(traj_time, t_init, lat, lon,start_age) :
    """ """
    # Load parameters and initialize variables.
    nturtle = len(t_init)
    traj_time[0,:] = start_age
    traj_time = traj_time - start_age
    duration = traj_time.flatten().max()
    n_steps = np.shape(traj_time)[0]
    t_step = traj_time[1,0]-traj_time[0,0]
    start = t_init.min()
    end = t_init.max() + duration
    final_duration = end - start
    final_n_steps = int(final_duration / t_step)+2

    # Build universal date vector.
    # (calendar zero 01/01/2002)
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