#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Rename FillValue of a variable.
Here used for VGPM data.
"""


import datetime
import netCDF4

# =============================================================================
# INPUTS
# =============================================================================
#Files names
data_dir = '/data/rd_exchange2/tcandela/STAMM/ressources/VGPM/data/'
prefix = 'npp.'
suffix = '.nc'
#
variable = 'npp'
current_FillValue = 'Hole_Value'
rename_FillValue = '_FillValue'
t0 = datetime.datetime(1997, 9, 14) #date of first file  # datetime.datetime(2002,1,1)+datetime.timedelta(185-1)
nb_files = 1000 #better too high
dt = 8 #dt between 2 files


# =============================================================================
# LOOP
# =============================================================================
time_origin = datetime.datetime(1950, 1, 1)
time = (t0 - time_origin).days * 24
current_date = t0
date_end = current_date + datetime.timedelta(days=(nb_files-1)*dt)
while (current_date <= date_end):
    day = (current_date - datetime.datetime(current_date.year,1,1)).days + 1
    day_str = str(("%03d") %day)
    file = data_dir + prefix + str(current_date.year) + day_str + suffix
    print(file)
    nc = netCDF4.Dataset(file,'r+')
    npp = nc.variables[variable]
    npp.renameAttribute(current_FillValue,rename_FillValue)
    
    if (current_date + datetime.timedelta(days=dt)).year != current_date.year:
        current_date = datetime.datetime(current_date.year + 1, 1, 1)
        
    else:
        current_date += datetime.timedelta(days=dt)
    time = (current_date - time_origin).days * 24
    nc.close()


