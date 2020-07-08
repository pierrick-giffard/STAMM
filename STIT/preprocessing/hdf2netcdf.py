#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Convert hdf to netcdf. Use egutknecht converter.

USE: python $stit/preprocessing/hdf2netcdf.py path_data prefix year0 year1 start end
"""


from datetime import datetime, timedelta
import sys
import subprocess

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


# Dates
t0 = datetime(year0,1,1) + timedelta(start-1) # date of first file
t1 = datetime(year1,1,1) + timedelta(end-1) # date of last file
nb_files = (t0 - t1).days
dt = 8 #dt in days between 2 files

# Converter
tool='/home/egutknecht/bin/h4tonccf_nc4'



# =============================================================================
# LOOP
# =============================================================================

current_date = t0
date_end = t1
while (current_date <= date_end):
    day = (current_date - datetime(current_date.year,1,1)).days + 1
    day_str = str(("%03d") %day)
    file = data_dir + prefix + str(current_date.year) + day_str + '.hdf'

    # Convert
    cmd = [tool, file]
    subprocess.run(cmd)
    
    # Increment time
    if (current_date + timedelta(days=dt)).year != current_date.year:
        current_date = datetime(current_date.year + 1, 1, 1)
        
    else:
        current_date += timedelta(days=dt)




