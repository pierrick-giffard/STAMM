# -*- coding: utf-8 -*-
"""
USE: Quick_stats.py   Outputfile.nc

Give some quick statistics about STAMM Outputs: mean, min, max and std of all variables.
"""


# =============================================================================
# IMPORTS
# =============================================================================
import numpy as np
import sys, os
import pandas as pd

#Personal librairies
sys.path.insert(1, os.path.join(sys.path[0], '../../LIB'))
#sys.path.append('C:\\Users\\pgiffard\\Desktop\\CODES\\stamm\LIB') #tmp
import turtle_lib as tul
import netCDF_lib as ncl



# =============================================================================
# DATA
# =============================================================================
infile = sys.argv[1]
mode = 'active'
key_alltracers = True
#infile = 'C:/Users/pgiffard/Desktop/test_ref.nc'#sys.argv[1]#  #tmp
OutFile = infile.replace(".nc","_Quick_stats.csv")
variables = ['traj_lat',
             'traj_lon',
             'u_current',
             'v_current']

if key_alltracers:
    variables.append('traj_temp')
    variables.append('traj_pp')
    
if mode == 'active':
    variables.append('u_swim')
    variables.append('v_swim')
    variables.append('PPmax')
    variables.append('habT')
    variables.append('habPP')
    variables.append('hab')
    variables.append('xgrad')
    variables.append('ygrad')
    variables.append('SCL')
dic = ncl.read_nc(infile, variables)


# =============================================================================
# STATISTICS
# =============================================================================

disabled_turtles, disabled_time = tul.find_disabled(infile)
print('Inactive turtles are not taken into account for diagnostics.')

stats = {}.fromkeys(variables)
stats['disabled_turtles'] = len(disabled_turtles)
for v in variables:
    print(v)
    array = tul.insert_nan(dic[v], disabled_turtles, disabled_time)[1:,:] #insert nan and remove first value
    Mean = np.nanmean(array)
    Min = np.nanmin(array)
    Max = np.nanmax(array)
    std = np.nanstd(array)
    stats[v] = [Mean, Min, Max, std]


# =============================================================================
# TABLE
# =============================================================================
columns=['Mean', 'Min', 'Max', 'std']
df = pd.DataFrame(stats, index=columns)
df = df.transpose()
df.to_csv(OutFile)
print(df)
print('\n')
print('*******************************************************************')
print("Wrote", OutFile)
print('*******************************************************************')    
       