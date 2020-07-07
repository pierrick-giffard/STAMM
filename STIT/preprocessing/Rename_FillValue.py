#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Rename FillValue of a variable.
Here used for VGPM data.

USE: python $stit/preprocessing/Rename_FillValue.py path_data prefix year0 year1 start end
"""


from  datetime import datetime, timedelta
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

# var names
variable = 'npp'
current_FillValue = 'Hole_Value'
rename_FillValue = '_FillValue'


# Dates
t0 = datetime(year0,1,1) + timedelta(start-1) # date of first file
t1 = datetime(year1,1,1) + timedelta(end-1) # date of last file
nb_files = (t0 - t1).days
dt = 8 #dt in days between 2 files


# =============================================================================
# LOOP
# =============================================================================
time_origin = datetime(1950, 1, 1)
time = (t0 - time_origin).days * 24
current_date = t0
date_end = t1
while (current_date <= date_end):
    day = (current_date - datetime(current_date.year,1,1)).days + 1
    day_str = str(("%03d") %day)
    file = data_dir + prefix + str(current_date.year) + day_str + suffix
    print(file)
    nc = netCDF4.Dataset(file,'r+')
    npp = nc.variables[variable]
    npp.renameAttribute(current_FillValue,rename_FillValue)
    
    if (current_date + timedelta(days=dt)).year != current_date.year:
        current_date = datetime(current_date.year + 1, 1, 1)
        
    else:
        current_date += timedelta(days=dt)
    time = (current_date - time_origin).days * 24
    nc.close()


