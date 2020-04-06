#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Change time value in netcdf and create it if it doesn't exist.
File name = data_dir + prefix + date + suffix
"""


import datetime
import netCDF4

# =============================================================================
# INPUTS
# =============================================================================
#Files names
data_dir = '/data/rd_exchange2/tcandela/STAMM/ressources/VGPM_seawifs/'
prefix = 'npp.'
suffix = '.nc'
time_var = 'time'
#
t0 = datetime.datetime(1997, 9, 14) #date of first file. VGPM: datetime.datetime(ystart,1,1)+datetime.timedelta(nbdays-1)
nb_files = 500
dt = 8 #dt in days between 2 files
#
digit = False #True for Ariane type files
maxsize = 4 #For Ariane type files, number of digits
#
vgpm = True #if instead of YYYY/MM/DD the day number appear, like in VGPM files.

# =============================================================================
# LOOP
# =============================================================================
time_origin = datetime.datetime(1950, 1, 1)
time = (t0 - time_origin).days * 24
if digit:
    for t in range(1,nb_files+1):
        ind = str(("%0" + str(maxsize)+ "d") %t)
        file = data_dir + prefix + ind + suffix
        print(file)
        nc = netCDF4.Dataset(file,'r+')
        time_var = nc.variables[time_var]
        time_var.setncattr("units","hours since 1950-01-01")
        time_var[:] = time
        time += dt*24
        nc.close()

else:
    current_date = t0
    date_end = current_date + datetime.timedelta(days=nb_files*dt)
    while (current_date <= date_end):
        
        if vgpm:
            day = (current_date - datetime.datetime(current_date.year,1,1)).days + 1
            day_str = str(("%03d") %day)
            date = str(current_date.year) + day_str
        else:
            date = current_date.strftime("y%Ym%md%d")
            try:
                file = data_dir + prefix + date + suffix
                nc = netCDF4.Dataset(file,'r')
                nc.close()
            except:
                date = current_date.strftime("%Y%m%d")   
        file = data_dir + prefix + date + suffix
        print(file)
        nc = netCDF4.Dataset(file,'r+')
        try:
            nc.createDimension(time_var, None)
            Time = nc.createVariable(time_var, 'i4', time_var)
            Time.units = "hours since 1950-01-01"
        except:
            Time = nc.variables[time_var]
        Time[:] = time
        
        if (current_date + datetime.timedelta(days=dt)).year != current_date.year:
            current_date = datetime.datetime(current_date.year + 1, 1, 1)
            
        else:
            current_date += datetime.timedelta(days=dt)
        time = (current_date - time_origin).days * 24
        nc.close()


