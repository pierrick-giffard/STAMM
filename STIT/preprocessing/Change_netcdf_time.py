#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Change time value in netcdf and create it if it doesn't exist.
File name = data_dir + prefix + date + suffix
For VGPM, date is the middle of the period of 8 days.

USE: python $stit/preprocessing/Change_netcdf_time.py path_data prefix year0 year1 start end
"""


from datetime import datetime, timedelta
import netCDF4
import sys

# =============================================================================
# INPUTS
# =============================================================================
#Files names
data_dir = sys.argv[1]
prefix = sys.argv[2]
year0 = int(sys.argv[3])
year1 = int(sys.argv[4])
start = int(sys.argv[5])
end = int(sys.argv[6])
suffix = '.nc'
time_var = 'time'

# Dates
t0 = datetime(year0,1,1) + timedelta(start-1) # date of first file
t1 = datetime(year1,1,1) + timedelta(end-1) # date of last file
nb_files = (t0 - t1).days
dt = 8 #dt in days between 2 files
#
digit = False #True for Ariane type files
maxsize = 4 #For Ariane type files, number of digits
#
vgpm = True #if instead of YYYY/MM/DD the day number appear, like in VGPM files.

# =============================================================================
# LOOP
# =============================================================================
time_origin = datetime(1950, 1, 1)
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
    date_end = t1
    while (current_date <= date_end):
        
        if vgpm:
            day = (current_date - datetime(current_date.year,1,1)).days + 1
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
        Time[:] = time + dt / 2 * 24 #middle of the period
        
        if (current_date + timedelta(days=dt)).year != current_date.year:
            current_date = datetime(current_date.year + 1, 1, 1)
            
        else:
            current_date += timedelta(days=dt)
        time = (current_date - time_origin).days * 24
        nc.close()


