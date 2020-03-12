#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 22 14:49:18 2019

@author: pgiffard

Calculates P0 in the ocean Atlantic, Pacific or Indian.
There are two options to select the area of interest:
    1-Select a box
    2-Use a mask with delimited basins. In this case, choose a basin and a latitude min/max.

Calculation method:
NPP is discretized in nb_bin ranges between 0 and npp_max. The number of grid points with npp in each range is calculated over all files,
and then the range for which the number of grid points is equal to 90% of the total number of grid points gives P0.
"""


########################## LIBRARIES ##########################################
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import datetime
import operator

########################## PARAMETERS #########################################
#Files names
DATADIR = '/data/rd_exchange2/tcandela/STAMM/ressources/VGPM/data/'
prefix = 'vgpm.'
suffix = '.nc'
#
t0 = datetime.datetime(2003, 1, 1) #date of first file. VGPM: datetime.datetime(ystart,1,1)+datetime.timedelta(nbdays-1)
dt = 8 #dt in days between 2 files
#
digit = False #True for Ariane type files
maxsize = 4 #For Ariane type files, number of digits
#
vgpm = True
#
nb_files = 42
var_name = 'npp'
unit = 'mgC/m2/day'

##Select the area of interest: 2 options
option = 1
#
#Option 1: rectangular area
if option == 1:
    lat_min = 0 #°N
    lat_max = 45  #°N
    lon_min = -180  #°E
    lon_max = -80   #°E
    #
    fillvalue = '_FillValue'
    #Mesh for lon/lat
    mesh = '/data/rd_exchange2/tcandela/STAMM/ressources/VGPM/VGPM_083_mesh.nc'
    lat_name = "latitude"
    lon_name = "longitude"

    
    


#Option 2: Mask basins
if option == 2:
    lat_min = -45 #°N
    lat_max = 45  #°N
    ocean = 'Pacific' #Atlantic, Pacific or Indian
    mask2 = '/homelocal/pgiffard/DATA/Mask_basins_regular_GLORYS025.nc'
    mask_var = 'mask'
    lat_name = "latitude"
    lon_name = "longitude"

#Choose percentile
percentile = 0.9 #90th percentile

#Calculation parameters
nb_bin = 100 #the higher the more accurate, calculation can be slow if too high
npp_max = 5000 #highest expected npp value




#################### Convert lon,lat to grid indices ##########################
if option == 1:
    nc = Dataset(mesh)
    lat = np.squeeze(nc.variables[lat_name])
    lon = np.squeeze(nc.variables[lon_name])
    lon = np.where(lon > 180, lon - 360, lon)
    #
    try:
        lat_min_i = np.abs(lat[:,0] - lat_min).argmin()
        lat_max_i = np.abs(lat[:,0] - lat_max).argmin()
        lon_min_i = np.abs(lon[0,:] - lon_min).argmin()
        lon_max_i = np.abs(lon[0,:] - lon_max).argmin()
        nblat = lat.shape[0]
        nblon = lat.shape[1]
    except:
        lat_min_i = np.abs(lat - lat_min).argmin()
        lat_max_i = np.abs(lat - lat_max).argmin()
        lon_min_i = np.abs(lon - lon_min).argmin()
        lon_max_i = np.abs(lon - lon_max).argmin()
        nblat = len(lat)
        nblon = len(lon)

        

elif option == 2:
    nc = Dataset(mask2)
    mask0 = np.squeeze(nc.variables[mask_var])
    lat = np.squeeze(nc.variables[lat_name])
    lon = np.squeeze(nc.variables[lon_name])
    lat_min_i = np.abs(lat[:,0] - lat_min).argmin()
    lat_max_i = np.abs(lat[:,0] - lat_max).argmin()

else:    
    raise ValueError("Option has to be set as 1 or 2")
    
if lat_min_i > lat_max_i:
    temp = lat_min_i
    lat_min_i = lat_max_i
    lat_max_i = temp
    



########################## MASK: area of interest  ############################
"""Mask has to be set to 1 in the area of interest and to NaN outside."""

if option == 1:  
    mask_rectangular = np.zeros((nblat, nblon))
    mask_rectangular[lat_min_i:lat_max_i+1, lon_min_i:lon_max_i+1] = 1
   

else:
    mask_lat = np.nan * np.ones((mask0.shape))
    mask_lat[lat_min_i:lat_max_i+1, :] = 1
    mask0 = mask0*mask_lat
    if ocean == 'Atlantic':
        mask = np.where((mask0 ==1 ) | (mask0 ==2) | (mask0 ==17) | (mask0 ==18) | (mask0 ==20),1,np.NaN)
    if ocean == 'Pacific':
        mask = np.where((mask0 == 4) | (mask0 == 5),1,np.NaN)
    if ocean == 'Indian':
        mask = np.where(mask0 ==3,1,np.NaN)
              

    


########################## NPP discretization #################################
discrete_npp = np.zeros(nb_bin)
current_date = t0
date_end = current_date + datetime.timedelta(days=nb_files*dt)
k=0
while (current_date <= date_end):
    k += 1
    #Filename
    if digit:
        ind = str(("%0" + str(maxsize)+ "d") %k)
        file = DATADIR + prefix + ind + suffix
    else:
        if vgpm:
            day = (current_date - datetime.datetime(current_date.year,1,1)).days +1
            day_str = str(("%03d") %day)
            date = str(current_date.year) + day_str
        else:
            date = current_date.strftime("y%Ym%md%d")
            try:
                file = DATADIR + prefix + date + suffix
                nc = Dataset(file,'r')
                nc.close()
            except:
                date = current_date.strftime("%Y%m%d")
        file = DATADIR + prefix + date + suffix
    ######### DATA ###########
    print(file)
    nc = Dataset(file)
    npp = nc.variables[var_name]
    f = operator.attrgetter(fillvalue)
    filval = f(npp)
    npp = np.squeeze(npp)
    mask = np.where(npp == filval, 0, 1)
    mask = mask * mask_rectangular
    mask = np.where(mask==0,np.NaN,1)
    npp = npp * mask  #npp is set to NaN outside the area of interest
    
    npp_area = npp[lat_min_i:lat_max_i+1, :]  
    flat_npp = npp_area.flatten()
    idx = np.where(np.isnan(flat_npp)==False)[0]
    flat_npp_area = flat_npp[idx]
    occurences, bins = np.histogram(flat_npp_area,bins=nb_bin, range=(0,npp_max))
    discrete_npp += occurences

    ######### Increment time ###########
    if (current_date + datetime.timedelta(days=dt)).year != current_date.year:
        current_date = datetime.datetime(current_date.year + 1, 1, 1)
        
    else:
        current_date += datetime.timedelta(days=dt)

discrete_npp_sum = np.cumsum(discrete_npp)       
nb_grid_points = discrete_npp_sum[-1] #nb_files * np.count_nonzero(mask[lat_min_i:lat_max_i+1, :] == 1)

########################## P0 CALCULATION #####################################
for i in range(nb_bin):
    if discrete_npp_sum[i] > percentile*nb_grid_points:
        P0 = (i+0.5) * npp_max/nb_bin       
        break
    if i == nb_bin-1:
        raise ValueError('Try with a higher value of npp_max')

########################## Plot area ##########################################
if len(lon.shape) == 1:
    lon,lat = np.meshgrid(lon,lat)    
fig=plt.figure(facecolor='w')
fig.set_size_inches(8,6)
plt.title('Area of interest',Fontsize=16,Fontweight='bold')
map = Basemap()
map.contourf(lon,lat,mask,colors='r')        
map.fillcontinents(color='grey')
map.drawcoastlines()
parallels = np.arange(-80,80,30.)
map.drawparallels(parallels,labels=[True,False,True,False],linewidth=0.5,fontsize=16,dashes=[1, 4])
meridians = np.arange(-180.,180,60.)
map.drawmeridians(meridians,labels=[True,False,False,True],linewidth=0.5,fontsize=16,dashes=[1, 4]) 
        
        
########################## Plot histogram #####################################
npp_plot = np.linspace(0.5*npp_max/nb_bin, npp_max-0.5*npp_max/nb_bin, nb_bin)
plt.rc('ytick',labelsize=13);plt.rc('xtick',labelsize=13)
fig=plt.figure()
fig.set_size_inches(8,4)
hist = plt.bar(npp_plot,discrete_npp/nb_grid_points,width=nb_bin/10,linewidth=4,color='g')
plt.xlabel('NPP (%s)'%unit, Fontsize=14, fontweight='bold')
plt.ylabel('%', Fontsize=14, fontweight='bold')
plt.show()
        
print('\n')
if option == 1:
    print('Option 1: rectangular area')
    print('Latitude min: %s°N' %lat_min)
    print('Latitude max: %s°N' %lat_max)
    print('Longitude min: %s°E' %lon_min)
    print('Longitude max: %s°E'%lon_max)

if option == 2:
    print('Option 2: ocean mask')    
    print('Ocean:',ocean)
    print('Latitude min: %s°N' %lat_min)
    print('Latitude max: %s°N' %lat_max)
    
print('Number of files:',nb_files)
print('Percentile:',percentile)
print('\n')
print('***************************')
print('==> P0 =', P0,unit)
print('***************************')
