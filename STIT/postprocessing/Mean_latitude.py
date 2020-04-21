# -*- coding: utf-8 -*-
"""
Plot mean latitude. 
It is possible to plot only a time subset (see parameters).
Select a latitude strip (latmin/latmax)
Only active turtles in the strip are considered.

USE:    python Mean_latitude.py output1.nc output2.nc output3.nc
        work with one or several output files
"""


# =============================================================================
# IMPORTS
# =============================================================================
import numpy as np
import netCDF4
import matplotlib.pyplot as plt
import sys, os
from pathlib import Path

#Personal librairies
sys.path.insert(1, os.path.join(sys.path[0], '../../LIB'))
#sys.path.append('C:\\Users\\pgiffard\\Desktop\\CODES\\stamm\LIB') #tmp
import netCDF_lib as ncl
import turtle_lib as tul

# =============================================================================
# PARAMETERS
# =============================================================================
#Define latitude min and max to select turtles only in this strip
latmin = 30
latmax = 50
#
outfile = 'auto' #set to 'auto' or choose name
#Time subset
t0 = 0 #first time to plot (in days)
tmax = 6500 #last time to plot (in days); -1 for last time
#
zone = 'auto' #'auto' or 'manual'. If manual, enter latmin_plot and latmax_plot

 
# =============================================================================
# INPUTS
# =============================================================================
infiles = sys.argv[1:]
#infile = 'C:/Users/pgiffard/Desktop/test_ref.nc'  #tmp
 



# =============================================================================
# ZONE
# =============================================================================
if zone == 'manual':
    latmin_plot = 30
    latmax_plot = 44



# =============================================================================
# FIGURE
# =============================================================================
plt.rc('ytick',labelsize=16);plt.rc('xtick',labelsize=16)
#
fig=plt.figure(facecolor='w')
fig.set_size_inches(10,8)
plt.title('Mean latitude',Fontsize=18,Fontweight='bold')
plt.xlabel('Years',fontsize=16)
plt.ylabel('Degrees North',fontsize=16)
plt.grid(linestyle='--')
if zone == 'manual':
    plt.ylim(latmin_plot, latmax_plot)
     
        



# =============================================================================
# LOOP OVER FILES
# =============================================================================
for file in infiles:
    namef = ncl.get_name(file)
    print(file)
    #data
    nc = netCDF4.Dataset(file)
    lat = np.squeeze(nc.variables['traj_lat'])
    lon = np.squeeze(nc.variables['traj_lon'])
    t_init = np.squeeze(nc.variables['init_t'])
    traj_time = np.squeeze(nc.variables['traj_time'])
    active = np.squeeze(nc.variables['active'])
    
    #1- select latitude strip and active turtles
    lat = tul.latitude_strip(lat, latmin, latmax)
    disabled_turtles, disabled_time = tul.find_disabled(file)
    lat = tul.insert_nan(lat, disabled_turtles, disabled_time)

    #2- reorder lat
    lat, lon, date_mat = ncl.age_to_date(traj_time, t_init, lat, lon)
    
    #3- compute mean latitude
    lat_m = np.nanmean(lat,axis=1)

# =============================================================================
# PLOT
# =============================================================================
    dt = max(t_init) - min(t_init)
    x = (min(t_init) + np.arange(len(lat_m)))/365
    xmax = int(x[t0:tmax][-1])
    #
    if namef == 'REF':
        plt.plot(x[t0:tmax], lat_m[t0:tmax], 'k--', linewidth = '2', label=namef)
    else:
        plt.plot(x[t0:tmax], lat_m[t0:tmax], linewidth = '1', label=namef)
    #
ax=plt.gca()
ax.xaxis.set_ticks(range(0,xmax+1,2))
ax.xaxis.set_ticklabels(range(0,xmax+1,2))
plt.legend(fontsize=14)



# =============================================================================
# OUTPUT
# =============================================================================
if outfile == 'auto':
    if len(infiles) == 1:
        outfile = file.replace(".nc","_Mean-Latitude.png")
    else:
        directory = str(Path(file).parents[0])
        outfile = directory + '/Mean-Latitude.png'
plt.savefig(outfile, format='png',bbox_inches='tight',dpi=300)
print('\n')
print('*******************************************************************')
print("Wrote", outfile)
print('*******************************************************************')

#plt.show() #tmp



