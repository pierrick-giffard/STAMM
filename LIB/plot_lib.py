#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 17 14:32:28 2019

@author: tcandela
"""

"""
Librairies to read nesting beach features and to compute initial positions
"""
# =============================================================================
# IMPORTS
# =============================================================================
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid.inset_locator import inset_axes
from matplotlib import gridspec
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
import time
import os
from os.path import isfile, join
import cv2
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import cm
import random

#Personal librairies
import librumeau as brum
import netCDF_lib as ncl
import turtle_lib as tul
# =============================================================================
# FUNCTIONS
# =============================================================================
def complete_release_map(infile, path, gridfile, lonlat1D, nturtles, xmin, xmax, ymin, ymax, lat_space, lon_space):

    figname = 'fig/Release_Map.png'

    #Zone to plot in minimap
    xmin_sl = -180
    xmax_sl =  180
    ymin_sl = -80
    ymax_sl =  85

    #yticks parameters
    ymax_pos = 11.5
    dy = 0.1
    
    #To annotate country names and change map scale, go to plot code#
    
    #Loading initial positions
    initfile = open(infile,'r')
    x_init, y_init, t_init = np.loadtxt(initfile,usecols=(0,1,3),unpack=True)
    x_init, y_init, t_init = x_init[:nturtles], y_init[:nturtles], t_init[:nturtles]
    print('\nInitial positions loaded')
    
    #Loading grid file
    grid = nc.Dataset(gridfile)
    if lonlat1D == True:
        lon_mat = np.squeeze(grid['glamt'])
        lat_mat = np.squeeze(grid['gphit'])
        
    else:
        lon_mat = np.squeeze(grid['glamt'])[0,:]
        lat_mat = np.squeeze(grid['gphit'])[:,0]
    print('\nGrid file loaded')
    
    #Convert grid point into lon/lat
    lon, lat = brum.grid_to_geo(x_init, y_init, lon_mat, lat_mat)
    
    for i in np.arange(len(lon)):
        if lon[i] > 180:
            lon[i] = lon[i]-360
            
    # PLOT CODE
    
    #Figure parameters
    fig = plt.figure(figsize=(7/2.54,5/2.54))
    ax_pos = fig.add_subplot(111)
    ax_pos_small = inset_axes(ax_pos, width="50%", height="25%", loc='upper right')
    
    ax_pos.set_yticks(np.arange(ymin,ymax_pos,dy))
    ax_pos.spines['right'].set_linewidth(0.5)
    ax_pos.spines['left'].set_linewidth(0.5)
    ax_pos.spines['bottom'].set_linewidth(0.5)
    ax_pos.spines['top'].set_linewidth(0.5)
    
    #Plot map    
    m = Basemap(ax=ax_pos, projection='merc',lon_0=5,lat_0=90.0,llcrnrlon=xmin,urcrnrlon=xmax,llcrnrlat=ymin,urcrnrlat=ymax,resolution='h')
    m.fillcontinents(color='0.65',alpha=1, lake_color='w')
    m.drawcoastlines(color='0.3',linewidth=0.2)
    m.drawcountries(color='w',linewidth=0.1)
    m.drawparallels(np.arange(ymin,ymax_pos,lat_space), labels=[1,0,0,0], fontsize=5, linewidth=0.1)
    m.drawmeridians(np.arange(xmin,xmax,lon_space), labels=[0,0,0,1], fontsize=5, linewidth=0.1)
    
    #Plot annotations
    m.drawmapscale(-55.5, 4.25,-55.5,4.25, length=75,fontsize=4.5, barstyle='simple',labelstyle='simple')
    ax_pos.annotate(u'French Guiana', xy=(0.60, 0.25), fontsize=5,fontweight='normal',xycoords='axes fraction' ,horizontalalignment='right', verticalalignment='center',color='k')
    ax_pos.annotate(u'Suriname', xy=(0.25, 0.4), fontsize=5,fontweight='normal',xycoords='axes fraction' ,horizontalalignment='right', verticalalignment='center',color='k')
    
    #PLot turtles
    lon_m,lat_m = m(lon,lat)
    m.scatter(lon_m,lat_m,color='b', marker="o",edgecolors='None',s = 0.1)#,alpha=0.8,marker=".")
        
    #Plot minimap
    m_pt = Basemap(ax=ax_pos_small, projection='cyl',lon_0=-180,lat_0=90.0,llcrnrlon=xmin_sl,urcrnrlon=xmax_sl,llcrnrlat=ymin_sl,urcrnrlat=ymax_sl,resolution='l')
    m_pt.fillcontinents(color='0.65',alpha=1, lake_color='w')
    m_pt.scatter([-53.85],[6.1], color='b', marker="o", edgecolors='None', s=4,zorder=12)
    
    #Save figure
    plt.subplots_adjust(left=0.1, right=0.95, bottom=0.02, top=0.98)
    plt.savefig(path+figname,dpi = 300)
    
    print('\nMission accomplie ! (plot saved at ' + path + figname +')\n')
    



def plot_habitat(ax,hab_mode, gridfile, numday,latlim,lonlim, SCL, To, food_max,dmin,dmax,param,data_lists,current_date,log=False) :
    """ Plot habitat on map."""
    # Read Temp end Mnk data.
    lonmax = max(lonlim)
    lat = param['lat_phy']
    lon = param['lon_phy']
    temp = param['T_var']
    U_var = param['U_var']
    V_var = param['V_var']
    food_path = param['food_dir'] + '/'
   
    #Feeding habitat
    if hab_mode == 'food' or hab_mode == 'tot':
        PP = ncl.interpolate_vgpm(current_date, param)
        if hab_mode == 'tot':
            PP = PP[::-1,:] #reverse lat
            PP = PP[119:,:] #remove first 10 degrees !!!!!! il faut sélectionner les indices pour que les grilles correspondent mais il peut y avoir un décalage d'une demi maille. La résolution doit être la même
        Food_hab = tul.food_hab(PP,food_max)
        #
        if hab_mode == 'food':
            lat = param['lat_food']
            lon = param['lon_food']
            latlon = ncl.read_nc(gridfile,[lat,lon])
            latmat = np.asarray(latlon[lat])
            lonmat = np.asarray(latlon[lon])

    
    #Temperature habitat
    if hab_mode == 'temp' or hab_mode == 'tot':
        T_files = data_lists[2]
        current_T_file = T_files[numday]
        T_dict = ncl.read_nc(current_T_file, [lat,lon,temp])
        T = np.squeeze(T_dict[temp])
        T_hab = tul.t_hab(T,SCL,To,param['species'])
        latmat = np.asarray(T_dict[lat])
        lonmat = np.asarray(T_dict[lon])
    
        
    
    #Ocean currents
    if hab_mode == 'current':
        U_files = data_lists[0]
        current_U_file = U_files[numday]
        U_dict = ncl.read_nc(current_U_file, [lat, lon, U_var])
        #
        V_files = data_lists[1]
        current_V_file = V_files[numday]
        V_dict = ncl.read_nc(current_V_file, [V_var])   
        #
        latmat = np.asarray(U_dict[lat])
        lonmat = np.asarray(U_dict[lon])
        #
        U = np.squeeze(U_dict[U_var])
        V = np.squeeze(V_dict[V_var])
        norm = np.sqrt((U**2)+(V**2))
        

    
    if hab_mode == 'food':
        hab = Food_hab
        legend = u"Foraging Habitat suitability index"
        cmap = 'pink_r'
        levels = np.arange(0.,1.1,0.1)
        ticks = levels
    elif hab_mode == 'temp':
        hab = T_hab
        legend = u"Thermal Habitat suitability index"
        cmap = 'pink_r'
        levels = np.arange(0.,1.1,0.1)
        ticks = levels
    elif hab_mode == 'tot':
        hab = T_hab*Food_hab
        legend = u"Habitat suitability index"
        cmap = 'pink_r'
        levels = np.arange(0,1.1,0.1)
        ticks = levels
    elif hab_mode == 'current':
        hab = norm
        legend = u'Current velocity [m/s]'
        cmap = 'pink_r'
        levels = np.arange(0,2.1,0.1)
        ticks = np.arange(0,2.2,0.2)
        hab = np.where(hab>levels[-1], levels[-1], hab)
   
    
    hab[hab<0.001] = float('nan')
    hab2 = np.column_stack((hab[:,max(np.where(lonmat[lonmat<lonmax])[0]):-2],hab[:,:max(np.where(lonmat[lonmat<lonmax])[0])]))
    lonmat2 = np.hstack((lonmat[max(np.where(lonmat[lonmat<lonmax])[0]):-2]-360,lonmat[:max(np.where(lonmat[lonmat<lonmax])[0])]))
    #
    im = ax.contourf(lonmat2,latmat,hab2,levels,cmap=cmap, alpha = 0.9,zorder=0)
    cbar = plt.colorbar(im, orientation='horizontal',pad = 0.1, shrink=0.87, ticks = ticks)#, shrink=0.9)#, shrink=0.45, pad=0.03, fraction=0.25)
    cbar.ax.tick_params(labelsize=12)
    cbar.set_label(legend, labelpad=5, size=16)
    #cbar.outline.set_linewidth(0.5)
    #cbar.ax.xaxis.set_tick_params(width=0.5)
    






def display_fig(frame_title=''):
    c=1
    fig=plt.figure(num=1,figsize=(11*c,6.21*c), facecolor='w', edgecolor='w')
    # Display frame title
    ax=fig.add_subplot(111)
    #cax=plt.axes([0.85, 0.1, 0.075, 0.8])
    #cax = fig.add_axes([0.85, 0.09, 0.045, 0.8])

    ax.text(0.15, 1.06,frame_title, ha='center',va='center', transform=ax.transAxes,fontweight = 'bold', color='k',fontsize=16,)
    '''
    cb_ax = fig.add_axes([0.05, 0.09, 0.045, 0.8])
    cb = pl.colorbar(im,cax=cb_ax,orientation='vertical')
    cb.solids.set(alpha=1)
    cb.ax.tick_params(labelsize=12)
    cb.set_label(u"habitat suitability index", labelpad=3, size=12)
    cb.outline.set_linewidth(0.5)
    cb.ax.xaxis.set_tick_params(width=0.5)
    '''
    return ax#,cax
       
def display_colorbar(f,im, ax_cb, label):
    cb=f.colorbar(im,cax=ax_cb,orientation='horizontal')
    cb.solids.set(alpha=1)
    cb.ax.tick_params(labelsize=8)
    cb.set_label(label, labelpad=1, size=8)
    cb.outline.set_linewidth(0.5)
    cb.ax.xaxis.set_tick_params(width=0.5)

def display_tracks(ax, lat='NA',lon='NA',dates='NA',ms=0.00,col='b', marker= 'o',alpha=0.5) :
    """ """
    ax.scatter(lon, lat, marker=marker,s=ms, edgecolor='none',c=col, alpha=alpha)

def plot_map(ax, latmin, latmax, lonmin, lonmax,value=0.6,res=0.25,alpha=1, lon_space=20,lat_space=10) :
    """ Plot continents. """

    map=Basemap(ax=ax,llcrnrlon=lonmin,llcrnrlat=latmin,urcrnrlon=lonmax,urcrnrlat=latmax,projection='cyl',resolution='l')
    #Get meridians and parallels spaces
    
    def getTicks(lmin, lmax, step):
        lmaxabs = (int(max(abs(lmax), abs(lmin)))/step+1)*step
        return np.intersect1d(np.arange(lmin+1, lmax), np.arange(-lmaxabs, lmaxabs, step))
    #Draw parallels & meridians
    map.drawparallels(np.arange(latmin, latmax, lat_space),labels=[1,0,0,0], fontsize=8,zorder=0.2,linewidth=0.1)
    map.drawmeridians(np.arange(lonmin, lonmax, lon_space),labels=[0,0,0,1], fontsize=8,zorder=0.2,linewidth=0.1)
    map.drawcountries(color='k',linewidth=0.01,zorder=0.3)
    map.drawcoastlines(color='grey',linewidth=0.2,zorder=0.3)
    map.fillcontinents(color='0.35')

def getTicks(lmin, lmax, step):
    lmaxabs = (int(max(abs(lmax), abs(lmin)))/step+1)*step
    return np.intersect1d(np.arange(lmin+1, lmax), np.arange(-lmaxabs, lmaxabs, step))

def show_start_point(ax, lat,lon) :
   """ """
   ax.plot((np.mean(lon[0,:]),),(np.mean(lat[0,:]),),markerfacecolor='w',
            markeredgecolor='k',marker='o',ms=6,mew=0.3,zorder=999)   
   
def plot_animation_frames(gridfile, dico,hab_mode,To,lethargy,coef_SMR,start_day,end_day,h,latlim,lonlim, save_path, param, data_lists, last_turtle, mortality, group, nb_cat, colors, dpi=100):
    """ Plot animation frames with turtles positions and approximate habitat. """  
    species = param['species']
    nturtles = param['nturtles'] - 1 if last_turtle == -1 else last_turtle
    #
    latmin = min(latlim)
    latmax = max(latlim)
    lonmin = min(lonlim)
    lonmax = max(lonlim)

    dmin = 80.
    dmax = 200. 
    lat = dico['traj_lat'][:,:last_turtle]          
    lon = dico['traj_lon'][:,:last_turtle]
    init_t = dico['init_t'][:last_turtle]
    traj_time = dico['traj_time'][:,:last_turtle]
    group = group[:last_turtle]
    #
    #Correction de certaines longitudes
    #lon[np.where(lon<200)]=lon[np.where(lon<200)]+360
    #if np.mean(lon[0,:]) < 36 :
        #Si point de départ dans l'Atlantique, on corrige l'affichage des longitudes > 180°E(cad toutes puisque 180°E est dans le Pacifique et que l'on se trouve dans l'Atlantique)
        # Lorsque les particules dépassent greenwich on ajoute 360 (elles passent de 359 à 361 plutot que de 359 à 1, par exemple en mediterranée)
    lon[lon>=lonmax]-=360
    
    if hab_mode != 'void' and mortality:
        temp = dico['traj_temp'][start_day:end_day,:last_turtle]
        date_death = tul.find_date_death(nturtles,temp,To,coef_SMR,lethargy,init_t, end_day-start_day)
    
    date_start_physfile = dt.datetime(param['ystart'],1,1)
    date_start_physfile_entier= date_start_physfile.toordinal()
    
    if hab_mode != 'void' and mortality:
        date_death_entier = date_death + date_start_physfile_entier
        
    month_names = ['Jan.','Feb.','Mar.','Apr.','May','Jun.','Jul.','Aug.','Sep.','Oct.','Nov.','Dec.']
    #
    SCL = tul.compute_SCL_VGBF(param['SCL0'], species, start_day)
    for step in range(start_day,end_day,h):
        print('\n')
        print(step, 'of', end_day-h)
        days_since_ref = int(step+init_t.min())
        date_title = date_start_physfile + dt.timedelta(days_since_ref)
        if param['time_periodic']:
            days_since_ref = int(days_since_ref%(param['time_periodic'] + init_t.min()))
            print("verifier date")
        date = date_start_physfile + dt.timedelta(days_since_ref)
        # Frame title
        date_today_entier = date.toordinal()
        m = '00'
        month = month_names[date_title.month - 1]
        day = str(("%02d") %date_title.day)
        year = str(date_title.year)
        title ='| '+day+' '+month+' '+year+' |'
        print(title)
        #
        newlat,newlon,date_mat = ncl.age_to_date(traj_time,init_t,lat,lon)
        #
        ax = display_fig(frame_title=title)
        # Display habitat.
        if hab_mode != 'void':
            # Calcul des paramètre relatifs à la nage active et à l'habitat
            SCL = tul.compute_SCL_VGBF(SCL, species, h) #increment SCL of h days
            food_max = tul.compute_Fmax(step+start_day,species,SCL,param['P0'])
            numday = days_since_ref - int(init_t.min())
            plot_habitat(ax, hab_mode, gridfile, numday, [latmin, latmax], [lonmin,lonmax], SCL, To, food_max, dmin, dmax, param, data_lists,date)


        # Find alive and dead turtles
        # Blue dots : alive turtles
        # Black dots: dead turtles
        # Dead turtles are removed from the animation 90 days after they died
        if hab_mode != 'void' and mortality and len(group)==0:
            index_dead_at_date = np.where((date_death_entier<=date_today_entier)&(date_death_entier+90>date_today_entier)) #+90 > dead disappear after 90 days
            index_alive_at_date = np.where(date_death_entier>date_today_entier)
        if hab_mode == 'void' and mortality and len(group)==0:
            index_dead_at_date=[]
            index_alive_at_date=np.arange(lat.shape[1])
            
        # Display position (scatter)
        if mortality and len(group)==0:
            display_tracks(ax, lat=newlat[step,index_dead_at_date],lon=newlon[step,index_dead_at_date],ms=11,col='k', marker = 'o',alpha=0.6)
            display_tracks(ax, lat=newlat[step,index_alive_at_date],lon=newlon[step,index_alive_at_date],ms=11,col='#1f78b4', marker = 'o',alpha=0.6)        
        elif len(group)==0:
            display_tracks(ax, lat=newlat[step,:],lon=newlon[step,:],ms=11,col='#1f78b4',alpha=0.6)
        
        if hab_mode != 'void' and len(group)>0:
            for cat in np.arange(nb_cat):
                if mortality:
                    index_dead_at_date = np.where((date_death_entier <= date_today_entier) & (date_death_entier + 90 > date_today_entier) & (group == cat)) #+90 > dead disappear after 90 days
                    index_alive_at_date = np.where((date_death_entier > date_today_entier) & (group == cat))

                    display_tracks(ax, lat=newlat[step,index_dead_at_date], lon=newlon[step,index_dead_at_date], ms=5, col=colors[cat], marker = 'x', alpha=0.6)
                    display_tracks(ax, lat=newlat[step,index_alive_at_date], lon=newlon[step,index_alive_at_date], ms=5, col=colors[cat], marker = 'o', alpha=0.6)
                else:
                    idx = np.where(group == cat)
                    display_tracks(ax, lat=newlat[step,idx], lon=newlon[step,idx], ms=5, col=colors[cat], marker = 'o', alpha=0.6)   
        # Plot starting point
        #show_start_point(ax, lat,lon)

        lon_space = (lonmax - lonmin)/7
        lat_space = (latmax - latmin)/7
        # Display map.
        plot_map(ax, latmin, latmax, lonmin, lonmax, lon_space,lat_space)
        plt.xlim([lonmin,lonmax])
        plt.ylim([latmin,latmax])
        
        #save figure
        m = str(("%04d") %step)
        plt.savefig(save_path + 'frame_' + m + '.png', bbox_inches='tight', dpi=dpi)
        plt.close()
 
       
def plot_animation_frames_tuned(gridfile, dico,hab_mode,To,lethargy,coef_SMR,start_day,end_day,h,latlim,lonlim,save_path, param, data_lists, last_turtle, mortality = True, dpi=100):
    """ 
    Plot animation frames for 1 turtle with a dt < 24h (for example 24 dt / day)
    Also plot 4 points at a distance grad_dx to see where gradients are computed
    Might work with several turtles
    """  
    #Tuned parameters
    delta = 24 #24 positions for 1 data file (dt = 1h)
    grad = False #to plot points where gradient is computed
    deg = 111195 #1degree = 111,195 km approx
    grad_dx = param['grad_dx']
    #
    species = param['species']
    nturtles = param['nturtles'] - 1 if last_turtle == -1 else last_turtle
    #
    latmin = min(latlim)
    latmax = max(latlim)
    lonmin = min(lonlim)
    lonmax = max(lonlim)

    dmin = 80.
    dmax = 200. 
    lat = dico['traj_lat'][:,:last_turtle]          
    lon = dico['traj_lon'][:,:last_turtle]
    init_t = dico['init_t'][:last_turtle]
    traj_time = dico['traj_time'][:,:last_turtle]
    
    date_start_physfile = dt.datetime(param['ystart'],1,1) #à modifier éventuellement
        
    month_names = ['Jan.','Feb.','Mar.','Apr.','May','Jun.','Jul.','Aug.','Sep.','Oct.','Nov.','Dec.']
    #
    SCL = param['SCL0'] + tul.age_to_SCL(start_day,species) #not exact if (SCL0 is not hatchling SCL and start_day > 0)
    for step in range(0,end_day,h): #here not days but time_steps
        print('\n')
        print(step, 'of', end_day-h)
        days_since_ref = int(init_t.min()) + 1 + step//delta #increment days each delta time steps
        date_title = date_start_physfile + dt.timedelta(days_since_ref)
        date = date_start_physfile + dt.timedelta(days_since_ref)
        # Frame title
        m = '00'
        month = month_names[date_title.month-1]
        day = str(("%02d") %date_title.day)
        year = str(date_title.year)
        title ='| '+day+' '+month+' '+year+' |'
        print('  ',title)
        print('File date : ',date.strftime("%d-%m-%Y"))
        #
        newlat,newlon,date_mat = ncl.age_to_date(traj_time,init_t,lat,lon)
        #
        ax = display_fig(frame_title=title)
        # Display habitat.
        if hab_mode != 'void':
            # Calcul des paramètre relatifs à la nage active et à l'habitat
            SCL = tul.compute_SCL_VGBF(SCL, species, 1)
            food_max = tul.compute_Fmax(step+start_day,species,SCL,param['P0'])
            numday = days_since_ref - int(init_t.min())
            plot_habitat(ax, hab_mode, gridfile, numday, [latmin, latmax], [lonmin,lonmax], SCL, To, food_max, dmin, dmax, param, data_lists,date)
            print(numday)

        display_tracks(ax, lat=newlat[step,:],lon=newlon[step,:],ms=11,col='#1f78b4',alpha=0.6)
        #For gradients points
        if grad:
            dx_lon =  grad_dx / (deg  * np.cos(newlat[step,0] * np.pi / 180))
            dx_lat =  grad_dx / deg
            display_tracks(ax, lat=newlat[step,:]-dx_lat,lon=newlon[step,:],ms=11,col='k',alpha=0.6)
            display_tracks(ax, lat=newlat[step,:]+dx_lat,lon=newlon[step,:],ms=11,col='k',alpha=0.6)
            display_tracks(ax, lat=newlat[step,:],lon=newlon[step,:]-dx_lon,ms=11,col='k',alpha=0.6)
            display_tracks(ax, lat=newlat[step,:],lon=newlon[step,:]+dx_lon,ms=11,col='k',alpha=0.6)
        
        # Plot starting point
        show_start_point(ax, lat,lon)

        lon_space = (lonmax - lonmin)/7
        lat_space = (latmax - latmin)/7
        # Display map.
        plot_map(ax, latmin, latmax, lonmin, lonmax, lon_space,lat_space)
        plt.xlim([lonmin,lonmax])
        plt.ylim([latmin,latmax])
        
        #save figure
        m = str(("%04d") %step)
        plt.savefig(save_path + 'frame_' + m + '.png', bbox_inches='tight', dpi=dpi)
        plt.close()
        
def convert_frames_to_video(pathIn, pathOut, fps):
    print('\n')
    print('****************************************************')
    print("Converting frames to video...")
    print('****************************************************')
    print('\n')
    frame_array = []
    files = [f for f in os.listdir(pathIn) if (isfile(join(pathIn, f)) and os.path.splitext(join(pathIn, f))[-1] == '.png')]
    
    #for sorting the file names properly
    #only png files should be in the directory
    try:
        files.sort(key = lambda x: int(x[6:10])) #work if name = frame_****.png
    except:
        files = sorted(files)
 
    for i in range(len(files)):
        time.sleep(0.01)
        filename = pathIn + files[i]
        #reading each files
        img = cv2.imread(filename)
        height, width, layers = img.shape
        size = (width,height)
        print(filename)
        #inserting the frames into an image array
        frame_array.append(img)
 
    out = cv2.VideoWriter(pathOut,cv2.VideoWriter_fourcc('M','J','P','G'), fps, size)
 
    for i in range(len(frame_array)):
        # writing to a image array
        out.write(frame_array[i])
    out.release()
    
def compute_presence_map(dx,x,y,xmin,xmax,ymin,ymax,seuil_numpos=30000,lat_space=5,lon_space=20):

    x = np.array(x)  
    x[x>=xmax]-=360
    y = np.array(y)
    
    x_ravel = x.ravel()
    y_ravel = y.ravel()
    x_ravel[np.where(y_ravel==0)] = 0
    y_ravel[np.where(x_ravel==0)] = 0
    x_ravel_nonzero = np.delete(x_ravel,np.where(x_ravel==0))
    y_ravel_nonzero = np.delete(y_ravel,np.where(y_ravel==0))

    x_ravel_nonone = np.delete(x_ravel_nonzero,np.where(x_ravel_nonzero==1.0))
    y_ravel_nonone = np.delete(y_ravel_nonzero,np.where(y_ravel_nonzero==1.0))
    x_ravel_nonone= np.append(x_ravel_nonone,0)
    x_ravel_nonone = np.append(x_ravel_nonone,410)
    y_ravel_nonone = np.append(y_ravel_nonone,-90)
    y_ravel_nonone = np.append(y_ravel_nonone,90)
    
    x,y = x_ravel_nonone, y_ravel_nonone

    heatmap, xedges, yedges = np.histogram2d(x,y,bins=(np.arange(xmin,xmax+dx,dx),np.arange(ymin,ymax+dx,dx)))
    extent = [xedges[0],xedges[-1], yedges[0],yedges[-1]]

    heatmap[np.where(heatmap>seuil_numpos)]=seuil_numpos

    #CMAP IMOS
    pimos = {'red':   ((0.00,  1.0,   1.0),
#                       (0.03,  0.624, 0.624),
                       (0.08,  0.000, 0.000),
                       (0.25,  0.000, 0.000),
#                       (0.33,  0.000, 0.000),
                       (0.42,  0.000, 0.000),
#                      (0.50,  0.000, 0.000),
                       (0.58,  1.000, 1.000),
                       (0.67,  1.000, 1.000),
#                       (0.75,  1.000, 1.000),
                       (0.83,  1.000, 1.000),
#                       (0.92,  0.651, 0.651),
                       (1.00,  0.471, 0.471)),
             'green': ((0.00,  1.0,   1.0),
#                      (0.03,  0.624, 0.624),
                       (0.08,  0.500, 0.500),
                       (0.25,  0.000, 0.000),
#                       (0.33,  0.000, 0.000),
                       (0.42,  0.749, 0.749),
#                       (0.50,  0.498, 0.498),
                       (0.58,  1.000, 1.000),
                       (0.67,  0.498, 0.498),
#                       (0.75,  0.000, 0.000),
                       (0.83,  0.000, 0.000),
#                       (0.92,  0.325, 0.325),
                       (1.00,  0.000, 0.000)),
             'blue':  ((0.00,  1.0,   1.0),
#                      (0.03,  0.624, 0.624),
                       (0.08,  1.000, 1.000),
                       (0.25,  1.000, 1.000),
#                       (0.33,  0.498, 0.498),
                       (0.42,  0.000, 0.000),
#                       (0.50,  0.000, 0.000),
                       (0.58,  0.000, 0.000),
                       (0.67,  0.000, 0.000),
#                       (0.75,  0.749, 0.749),
                       (0.83,  0.000, 0.000),
#                       (0.92,  0.235, 0.235),
                       (1.00,  0.000, 0.000))}

    cmap_imos = LinearSegmentedColormap('cmap_imos', pimos)
    cmap_imos.set_bad(color='w', alpha=None)

    return heatmap,extent,cmap_imos

def plot_presence_map(heatmap,extent,cmap,xmin,ymin,xmax,ymax,lat_space,lon_space,imageName,mode,traj_lat, traj_lon,dx,seuil_numpos):
    c = 0.89
    f = plt.figure(figsize = (12*c/2.54,8*c/2.54))
    gs = gridspec.GridSpec(2,1,height_ratios=[11,1],left=0.08, right=0.98, bottom=0.08, top=0.95)
    ax = plt.subplot(gs[0])
    ax_cb = plt.subplot(gs[1])
    
    if mode == 'Density':
        im = ax.imshow(heatmap.T/1000, extent=extent, origin='lower',interpolation='None',cmap=cmap, vmax=seuil_numpos/1000)
    else:
        im = ax.imshow(heatmap.T, extent=extent, origin='lower',interpolation='None',cmap=cmap)
    
    map=Basemap(ax=ax,llcrnrlon=xmin,llcrnrlat=ymin,urcrnrlon=xmax,urcrnrlat=ymax,projection='cyl',resolution='l')
    map.drawparallels(getTicks(ymin, ymax, lat_space), labels=[1,0,0,0], fontsize=8, linewidth=0.2)
    map.drawmeridians(getTicks(xmin, xmax, lon_space), labels=[0,0,0,1], fontsize=8, linewidth=0.2)
    map.drawcoastlines(color='grey',linewidth=0.2)
    map.drawcountries(color='k',linewidth=0.01)
    map.fillcontinents(color='0.35')
    
    ax.spines['right'].set_linewidth(0.5)
    ax.spines['left'].set_linewidth(0.5)
    ax.spines['bottom'].set_linewidth(0.5)
    ax.spines['top'].set_linewidth(0.5)
    
    show_start_point(ax, traj_lat, traj_lon)
    cb=f.colorbar(im,cax=ax_cb,orientation='horizontal')
    cb.solids.set(alpha=1)
    cb.ax.tick_params(labelsize=8)
    
    if mode == 'Density':
        cb.set_label(u"Thousands of turtle days per %s°x %d° box"%(dx,dx), labelpad=1, size=8)
    else:
        cb.set_label(u"Number of cold-induced death events %s°x %d° box"%(dx,dx), labelpad=1, size=8)

    cb.outline.set_linewidth(0.5)
    cb.ax.xaxis.set_tick_params(width=0.5)
    
    plt.savefig(imageName,bbox_inches='tight', dpi = 800)  



def display_trajectories(dataFile,f,ax,xmin,lw=0.005, ms=0.04, col='b', alpha=0.5,show_start=True):

    """
    Highly parametrized display function, some settings however are still
    hard-coded (colorbar).
    
    var : variable to display with color axis, if None, all trajectories will be displayed using a uniform color
    lw: linewidth
    ms: markersize
    col: color
    alpha: transparency

    show_cb : set to True to plot colorbar
    show_start : set to True to plot average start point of plotted trajectories
    """

    nb_days,nb_traj=dataFile.lon.shape

    #########################
    # DISPLAYS TRAJECTORIES #
    #########################

    x = np.array(dataFile.lon)  
    x[x<=xmin]+=360
    y = np.array(dataFile.lat)

#    if np.mean(dataFile.lon[0,:]) < 36 :
#        #Si point de départ dans l'Atlantique, on corrige l'affichage des longitudes > 180°E(cad toutes puisque 180°E est dans le Pacifique et que l'on se trouve dans l'Atlantique)
#        x = np.array([l+(((1-np.sign(l))/2)*360) for l in x]) 
#        # Lorsque les particules dépassent greenwich on ajoute 360 (elles passent de 359 à 361 plutot que de 359 à 1, par exemple en mediterranée)
#        x[np.where(x<200)]=x[np.where(x<200)]+360
        

    colmin = None
    colmax = None
    
    cb_title=None
    colorpalette = cm.get_cmap('jet')

    time = np.array(dataFile.traj_time)
    if np.max(time)>3*365:
        col = (time/365.).flat
    else:
        col = time.flat
        
    xbis = x.flat
    ybis = y.flat
    
    #Sample points if too many to avoid long computations
    nsample = int(1e6)
    if len(xbis) > nsample:
        rand_ind = random.sample(range(len(xbis)),nsample)
        sort_ind = rand_ind.sort()
        x = xbis[rand_ind]
        y = ybis[rand_ind]
        col = col[rand_ind]
        p=ax.scatter(x, y,c=col, s=ms, edgecolor='none',cmap=colorpalette, vmin=colmin, vmax=colmax, alpha=alpha)

    else:
        p=ax.scatter(xbis, ybis, c=col, s=ms, edgecolor='none',cmap=colorpalette, vmin=colmin, vmax=colmax, alpha=alpha)


    #Shows starting point
    dataFile.lon = np.array([l+(((1-np.sign(l))/2)*360) for l in dataFile.lon])
    ax.plot((np.mean(dataFile.lon[0,:]),),(np.mean(dataFile.lat[0,:]),),markerfacecolor='w',markeredgecolor='k',marker='o',ms=6,mew=0.3,zorder=999)
    
    return p,time


def display_trajectories_particular(lon, lat, traj_time, f, ax, lw=0.005, ms=0.04, col='b', alpha=0.5, show_start=True):

    """
    Highly parametrized display function, some settings however are still
    hard-coded (colorbar).
    
    var : variable to display with color axis, if None, all trajectories will be displayed using a uniform color
    lw: linewidth
    ms: markersize
    col: color
    alpha: transparency

    show_cb : set to True to plot colorbar
    show_start : set to True to plot average start point of plotted trajectories
    """

    nb_days,nb_traj = lon.shape
#        #Si point de départ dans l'Atlantique, on corrige l'affichage des longitudes > 180°E(cad toutes puisque 180°E est dans le Pacifique et que l'on se trouve dans l'Atlantique)
#        x = np.array([l+(((1-np.sign(l))/2)*360) for l in x]) 
#        # Lorsque les particules dépassent greenwich on ajoute 360 (elles passent de 359 à 361 plutot que de 359 à 1, par exemple en mediterranée)
#        x[np.where(x<200)]=x[np.where(x<200)]+360
        

    colmin = None
    colmax = None
    
    cb_title=None
    colorpalette = cm.get_cmap('jet')

    time = np.array(traj_time)
    if np.max(time)>3*365:
        col = (time/365.).flat
    else:
        col = time.flat
        
    xbis = x.flat
    ybis = y.flat
    
    #Sample points if too many to avoid long computations
    nsample = int(1e6)
    if len(xbis) > nsample:
        rand_ind = random.sample(range(len(xbis)),nsample)
        sort_ind = rand_ind.sort()
        x = xbis[rand_ind]
        y = ybis[rand_ind]
        col = col[rand_ind]
        p=ax.scatter(x, y,c=col, s=ms, edgecolor='none',cmap=colorpalette, vmin=colmin, vmax=colmax, alpha=alpha)

    else:
        p=ax.scatter(xbis, ybis, c=col, s=ms, edgecolor='none',cmap=colorpalette, vmin=colmin, vmax=colmax, alpha=alpha)


    #Shows starting point
    lon = np.array([l+(((1-np.sign(l))/2)*360) for l in lon])
    ax.plot((np.mean(lon[0,:]),),(np.mean(lat[0,:]),),markerfacecolor='w',markeredgecolor='k',marker='o',ms=6,mew=0.3,zorder=999)
    
    return p,time


def swimming_rose(theta, speed,  bins_classes,  save, nb_directions=32, title=False, pmax='auto', mode='height'):
    """
    theta: velocity directions
    speed: velocity norm
    bins_classes: velocity categories. Ex: [0,0.05,0.1,0.2,0.3,0.5,2]
    save: False or figure name
    nb_directions: divide 360° in nb_directions parts
    title: figure title
    pmax: maximum percentage for yaxis
    mode: area or height
    """
    rose, directions = brum.compute_directions(speed, theta, nb_directions, bins_classes)
    
    #useful variables for plot
    delta = 2 * np.pi / nb_directions - 0.05 # Angle du segment de cercle
    width = np.ones(nb_directions) * 2 * np.pi /nb_directions - 0.05
    colors = ['blue', 'deepskyblue', 'greenyellow', 'yellow', 'orange', 'red', 'purple', 'k'] #max 8 classes
    labels = ['%.2f - %.2f m/s'%(bins_classes[k], bins_classes[k+1]) for k in range(len(bins_classes)-1)]
    labels[-1] = '> %.1f m/s'%bins_classes[-2]
    
    #decoration
    fontsize = 20
    plt.rc('ytick',labelsize=fontsize); plt.rc('xtick',labelsize=fontsize)
    plt.figure(figsize=(10,12)) 
    ax = plt.subplot(111, projection='polar')    
    
    #plot
    bottom=0
    nb_data = np.sum(~np.isnan(theta))
    for k in range(len(bins_classes) - 1):
        height = brum.height_trapezoid(rose[:, k], bottom, delta) if mode == 'area' else rose[:, k]
        height *= 100 / nb_data 
        ax.bar(directions, height, width=width, bottom=bottom, edgecolor='k', color = colors[k], label=labels[k])
        bottom += height
    
    #maximum % for yticks
    if pmax == 'auto':
        pmax = np.max([np.sum(rose[k, :]) for k in range(nb_directions)])
    
    if mode == 'height':
        #plot mean direction
        mean_dir = brum.mean_direction(theta)
        ax.plot([mean_dir, mean_dir],[0,pmax], '-.',color='k',linewidth =5)
    
    #legend, title and ticks
    plt.legend(fontsize=fontsize, loc=(1,0))
    if title:
        plt.title(title, fontsize=fontsize+4, fontweight='bold')
    ax.set_xticklabels(('E','NE','N','NW','W','SW','S','SE'),fontsize=fontsize)    
    yticks = np.arange(2, pmax + 2, 2, dtype='int')
    ax.set_yticks(yticks)
    ax.set_yticklabels([str(i) + '%' for i in yticks], fontsize=fontsize, fontweight = 'bold')

    if save:
        save_fig(save)


def save_fig(fname, dpi=200):
    plt.savefig(fname, bbox_inches='tight', dpi=dpi)
    print ('Saved ', fname)
    plt.show()
    plt.close()
