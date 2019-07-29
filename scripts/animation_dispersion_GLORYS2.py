#!/usr/bin/env python
#-*- coding:utf-8 -*-
import pylab as pl
import numpy as np
import time
import sys
import os
import datetime as dt
import matplotlib.cm as cm
from mpl_toolkits.basemap import Basemap, shiftgrid
import netCDF4 as nc
import matplotlib.pyplot as plt
#personnal libraries


# ||===========================================================================
# ||                                                                          || 
# || module description: plot frames with turtles positions over habitat map  ||
# ||                     from a netCDF IBM2D simulation file                  ||
# ||                                                                          ||
# ||                                                                          ||
# ============================================================================||

########################
# Fonction utilitaires
########################

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

def find_date_death(turtles,temp,To,coef_SMR,lethargy,init_t,days=6572.):

    print('Calcul des dates de mort')
    #Age et Tmin fonction de lage
    age_year = np.arange(days)/365. 

    Tmin = To - coef_SMR*0.21*np.sqrt(0.000214*(1.43*(1-np.exp(-0.226*(age_year+0.17)))*100)**2.86)
    
    #On vire la première donnée de temp
    temp[0,:] = temp[1,:]


    date_death = []
    a = 0
    d = 0
    #On boucle sur les individus passés dans la zone donnée (Nord v. Sud)
    for i in range(turtles):
        temp[:,i]=temp[:,i]-Tmin
        temp_serie = temp

        c=1
        #When T>Tmin, tdiff>0 ==> NaN for this day
        #When T<Tmin, tdiff<0 ==> Number of days in row where T<Tmin
        for j in range(0,np.shape(temp)[0]):
            if temp[j,i]>0:
                temp_serie[j,i] = np.nan
                c=1
            else : 
                #print("Tdiff",k,i,c)
                temp_serie[j,i] = c
                c=c+1

        #Search the first day when T<Tmin (after a defined period of time, e.g 1 or 10 days)
        if np.shape(np.where(temp_serie[:,i]==lethargy))[1]>0:
            #print('Day death', np.min(np.where(temp_serie[:,i]==lethargy)))
            date_death.append(np.min(np.where(temp_serie[:,i]==lethargy)[0][:])+init_t[i])
            d=d+1
        #Ou s'il est resté dans des eaux "chaudes"
        else :
            date_death.append(1e34)
            a =a+1
  
    print('     alive', a)
    print('     dead',d)
    print('End fin dead  turtles')

    date_death = np.array(date_death)
    return date_death






##################
#Biologie
##################
def age_to_SCL(age,species) :

    """ Compute SCL for a given age using a model proposed in F. Perham et
    al., Age and growth of Loggerhead sea turtle of coastal Georgia, 1997"""
    
    if species == 'loggerhead':
        #######Loggerhead
        ## VB Parham & Zug 1997
        A = 1.088
        B = 0.9649
        k = 0.0739
        SCL = A*(1 - B*np.exp(-k*age/365.))
    
    elif species == 'leatherback':
        #######Leatherback
        SCL = 1.43*(1-np.exp(-0.226*(age/365.+0.17)))
    
    return SCL

def compute_M(species, SCL):
    """
    Compute all turtles mass(kg)
    """
    if species == 'leatherback':
        """
        Ref : Jones, T.T., Hastings, M.D., Bostrom, B.L., Pauly, D., Jones, D.R., 2011. Growth of captive leatherback turtles, Dermochelys coriacea, with inferences on growth in the wild: Implications for population decline and recovery. Journal of Experimental Marine Biology and Ecology 399, 84–92.
        """
        M = 0.000214*(SCL*100)**2.86
        return M

    elif species == 'loggerhead':
        A= 0.00034
        b = 2.89
        M = A*(SCL*100)**b
        return M

def food_hab(F,F_max) :
    """ Compute food habitat"""
    Food_hab = F/F_max
    Food_hab[Food_hab>1] = 1
    return Food_hab

def compute_Fmax(age,tracer,species,SCL,Fa):
    """
    Compute food threshold
    Ref : Jones, T.T., Bostrom, B.L.,  Hastings, M.D., Van Houtan K.S., Pauly, D., 2012. Resource requirements of the Pacific Leatherback Turtle Population
    """

    if tracer == "mnk" :
        Fmax= (SCL/0.325)

    elif tracer == "PP":
        if species == 'loggerhead' : 
            #########Loggerhead - 75
            Fmax =0.195*(((1-np.exp(-0.0981*(age/365.+0.36)))**(1.5))*(np.exp(-0.0981*(age/365.+0.36))))/(1-(1-np.exp(-0.0981*(age/365.+0.36)))**(0.195))
            Fmax = Fmax*float(Fa)
            print(Fmax)

        elif species == 'leatherback':
            #########Leatherback
            Fmax = 312*2.86*0.299*(((1-np.exp(-0.299*(age/365.+0.17)))**(2.86-1))*(np.exp(-0.299*(age/365.+0.17))))/(1-(1-np.exp(-0.299*(age/365.+0.17)))**(2.86*0.0328))
            Fmax = Fmax/(2835.24/float(Fa))
    return Fmax



def t_hab(T,SCL,To,species) :
    """ Compute temperature habitat """
    T_hab = np.ones([np.shape(T)[0],np.shape(T)[1]])

    if species == 'leatherback' :
        Mass = 0.000214*(SCL*100)**2.86
        #########################
        ## GAUSSIENNE ASYMETRIQUE
        #########################
        
        Topt = 24 - 0.21*np.sqrt(Mass)
        Tmin = 24 - 1.05*np.sqrt(Mass)
        sigma = (Topt-Tmin)/2
        inf = np.where(T <= Topt)
        sup = np.where(T>Topt)
        T_hab[inf]=np.exp(-2*((T[inf]-Topt)/(Topt-Tmin))**2)
        T_hab[sup]=1.0
        
    elif species == 'loggerhead':

        T1 = 14.
        T2 = 18.
        T3 = 29.
        T4 = 32.

        mid = np.where((T>=T2)&(T<=T3))
        sup = np.where(T > T3)
        inf = np.where(T < T2)
        T_hab[inf]=np.exp(-2*((T[inf]-T2)/(T2-T1))**2)
        T_hab[mid]=1.0
        T_hab[sup]=np.exp(-2*((T[sup]-T3)/(T4-T3))**2)


    return T_hab

##################
#Visualisation
##################
#def plot_habitat(ax,cax,GLORYS_path,numday,latlim,lonlim, SCL, To, food_max,dmin,dmax,tracer,species,log=False) :
def plot_habitat(ax,GLORYS_path,numday,latlim,lonlim, SCL, To, food_max,dmin,dmax,tracer,species,log=False) :
    """ Plot habitat on map."""
    # Read Temp end Mnk data.
    lonmin = min(lonlim)
    lonmax = max(lonlim)
    latmin = min(latlim)
    latmax = max(latlim)
    numday = '0000'+str(numday)
    numday = numday[-4:]
    #mask_path = GLORYS_path+'/mesh_grid/mesh_reg.nc'
    mask_path = GLORYS_path+'/sad/glorys2v4/mesh_reg.glorys2v4.nc'

    #Feeding habitat
    if tracer == "mnk":
        food_path = '/data/MEMMS/IBM/GLORYS/mnk/mnk_'+str(numday)+'.nc'
        food_dict = read_nc(food_path,['latitude','longitude','mnk'])
        mnk = np.asarray(food_dict['mnk'])[0,:,:]
        Food_hab = food_hab(mnk,food_max)

    elif tracer == "PP":
        #food_path = GLORYS_path+'/forcings/PP/PP_'+str(numday)+'.nc'
        food_path = GLORYS_path+'/forcings/dddd/bgch/npp/cmems_R2018_npp_'+str(numday)+'.nc'
        #food_dict = read_nc(food_path,['latitude','longitude','pp'])
        food_dict = read_nc(food_path,['latitude','longitude','npp_vgpm'])
        PP = np.asarray(food_dict['npp_vgpm'])[0,:,:]
        Food_hab = food_hab(PP,food_max)

    #T_path = GLORYS_path+'/forcings/reg_sosie/T/GLORYS_'+str(numday)+'_gridT.nc'
    T_path = GLORYS_path+'/forcings/dddd/phys/mercatorglorys2v4_T_'+str(numday)+'.nc'
    mask_dict = read_nc(mask_path,['mask','x','y'])
    #T_dict = read_nc(T_path,['votemper']) 
    T_dict = read_nc(T_path,['temperature']) 
    mask = np.asarray(mask_dict['mask'])[:,:]
    mask[mask==1][:] = float("nan")
    x_mask = np.asarray(mask_dict['x'])
    y_mask = np.asarray(mask_dict['y'])
    latmat = np.asarray(food_dict['latitude'])
    lonmat = np.asarray(food_dict['longitude'])

    #Flip de la matrice de latitude
    #latmat = np.flipud(latmat)

    # Compute temperature and food habitats.
    #T = np.asarray(T_dict['votemper'])[0,:,:]
    T = np.asarray(T_dict['temperature'])[0,:,:]
    T_hab = t_hab(T,SCL,To,species)
    
    #Habitat total
    hab = T_hab*Food_hab

    hab[hab<0.001] = float('nan')
    hab2 = np.column_stack((hab[:,max(np.where(lonmat[lonmat<lonmax])[0]):-2],hab[:,:max(np.where(lonmat[lonmat<lonmax])[0])]))
    lonmat2 = np.hstack((lonmat[max(np.where(lonmat[lonmat<lonmax])[0]):-2]-360,lonmat[:max(np.where(lonmat[lonmat<lonmax])[0])]))
    #levels = np.arange(0,1.1,0.1)
    levels = np.arange(0,1.1,0.1)
    im = ax.contourf(lonmat2,latmat,hab2,levels,cmap='pink_r',alpha = 0.9,zorder=0)
    cbar = plt.colorbar(im, orientation='horizontal',pad = 0.1, shrink=0.87)#, shrink=0.9)#, shrink=0.45, pad=0.03, fraction=0.25)
    cbar.ax.tick_params(labelsize=12)
    cbar.set_label(u"Habitat suitability index", labelpad=5, size=16)
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
    


    

def display_tracks(ax, lat='NA',lon='NA',dates='NA',ms=0.00,col='b',alpha=0.5) :
    """ """
    p=ax.scatter(lon, lat, marker='o',s=ms, edgecolor='none',c=col, alpha=alpha)

    

def plot_map(ax, latmin, latmax, lonmin, lonmax, lon_space=20,lat_space=10) :
    """ Plot continents. """

    map=Basemap(ax=ax,llcrnrlon=lonmin,llcrnrlat=latmin,urcrnrlon=lonmax,
                urcrnrlat=latmax,projection='cyl',resolution='l')
    #Get meridians and parallels spaces
    def getTicks(lmin, lmax, step):
        lmaxabs = (int(max(abs(lmax), abs(lmin)))/step+1)*step
        return np.intersect1d(np.arange(lmin+1, lmax), np.arange(-lmaxabs, lmaxabs, step))
    #Draw parallels & meridians
    map.drawparallels(getTicks(latmin, latmax, lat_space),
                     labels=[1,0,0,0], fontsize=14,zorder=1)
    map.drawmeridians(getTicks(lonmin, lonmax, lon_space),
                        labels=[0,0,0,1], fontsize=14,zorder=1)
    map.drawcoastlines(linewidth=0.2)
    map.fillcontinents()

def show_start_point(ax, lat,lon) :
   """ """
   ax.plot((np.mean(lon[0,:]),),(np.mean(lat[0,:]),),markerfacecolor='w',
            markeredgecolor='k',marker='o',ms=12,mew=1.,zorder=999)




###############
# Production des images pour l'animation
###############
    
def plot_animation_frames(path,dico,To,lethargy,coef_SMR,Fa,start_day,end_day,turtles,h,latlim,lonlim,lat_space,lon_space,place,tracer,mode,species,start_age) :
    """ Plot animation frames with turtles positions and approximate habitat. """  
    latmin = min(latlim)
    latmax = max(latlim)
    lonmin = min(lonlim)
    lonmax = max(lonlim)

    dmin = 80.
    dmax = 200. 
    lat = dico['traj_lat'][0:-1,:]          
    lon = dico['traj_lon'][0:-1,:]

    #Correction de certaines longitudes
    lon[np.where(lon<200)]=lon[np.where(lon<200)]+360
    #if np.mean(lon[0,:]) < 36 :
        #Si point de départ dans l'Atlantique, on corrige l'affichage des longitudes > 180°E(cad toutes puisque 180°E est dans le Pacifique et que l'on se trouve dans l'Atlantique)
        # Lorsque les particules dépassent greenwich on ajoute 360 (elles passent de 359 à 361 plutot que de 359 à 1, par exemple en mediterranée)
        #lon[np.where(lon<200)]=lon[np.where(lon<200)]+360

    temp = dico['traj_temp'][0:-1]
    init_t = dico['init_t']
    traj_time = dico['traj_time'][0:-1,:]
    date_death=find_date_death(turtles,temp,To,coef_SMR,lethargy,init_t, 6572.)#, end_day-start_day)
    
    #date_start_physfile = dt.datetime(1998,5,17,)
    date_start_physfile = dt.datetime(2002,1,1,)
    date_start_physfile_entier= date_start_physfile.toordinal()
    date_death_entier = date_death+date_start_physfile_entier
    
    month_names = ['Jan.','Feb.','Mar.','Apr.','May','Jun.','Jul.','Aug.','Sep.','Oct.','Nov.','Dec.']
    for step in range(start_day,end_day,h) :      
        days_since_ref1 = int(step+init_t.min()) # for title
        days_since_ref = days_since_ref1%2557 + 1 # looping for GLORYS
        t1=time.clock()
        # Frame title.
        date = date_start_physfile + dt.timedelta(days_since_ref1)
        date_today_entier = date.toordinal()
        m = '00'
        #month = m + str(date.month)
        #month = month[-2:]
        print(date.month-1)
        month = month_names[date.month-1]
        print(month)
        day = m + str(date.day)

        day = day[-2:]
        year = str(date.year-2001)
        #title ='| '+month+' / '+day+'  year '+year+' |'
        title ='| '+month+' '+day+', year '+year+' |'
        print(title)
        newlat,newlon,date_mat = age_to_date(traj_time,init_t,lat,lon,start_age)

        #Compute num. of physical files
        #GLORYS_path = '/Users/baptistemourrain/Desktop/ModeleStamm/glorys1'       
        GLORYS_path = '/data/FSHML/Tortues/LMTL-WP4/mercatorglorys2v4_1998_2015'       
        n_day =  '0000'+str(step+1)
        n_day = n_day[-4:]
        ax = display_fig(frame_title=title)
        # Display habitat.
        if mode == 'active':
            # Calcul des paramètre relatifs à la nage active et à l'habitat
            SCL = age_to_SCL(step+start_age,species)
            M   = compute_M(species,SCL)
            food_max = compute_Fmax(step+start_age,tracer,species,SCL,Fa)
            #ax,cax = display_fig(frame_title=title)
            #plot_habitat(ax,cax,GLORYS_path, days_since_ref,[latmin,latmax],[lonmin,lonmax],SCL,To,food_max,dmin,dmax,tracer,species)
            plot_habitat(ax,GLORYS_path, days_since_ref,[latmin,latmax],[lonmin,lonmax],SCL,To,food_max,dmin,dmax,tracer,species)

        
        # Find alive and dead turtles
        # Blue dots : alive turtles
        # Black dots: dead turtles
        # Dead turtles are removed from the animation 90 days after they died
        index_dead_at_date = np.where((date_death_entier<=date_today_entier)&(date_death_entier+90>date_today_entier))
        index_alive_at_date = np.where(date_death_entier>date_today_entier)
        
        #print(np.shape(index_dead_at_date), np.shape(index_alive_at_date))
        
        # Display position (scatter)
        display_tracks(ax, lat=newlat[step,index_dead_at_date],lon=newlon[step,index_dead_at_date],ms=11,col='k',alpha=0.6)
        display_tracks(ax , lat=newlat[step,index_alive_at_date],lon=newlon[step,index_alive_at_date],ms=11,col='#1f78b4',alpha=0.6)

        # Plot starting point
        show_start_point(ax, lat,lon)


        # Display map.
        plot_map(ax, latmin, latmax, lonmin, lonmax, lon_space,lat_space)
        plt.xlim([lonmin,lonmax])
        plt.ylim([latmin,latmax])

        # Saving frame (if non existent, the specified directory is created)
        #if os.path.isdir(path+place) == False : 
            #os.mkdir(path+place)
        m = '0000'+str(step)
        m = m[-4:]
        
        #plt.savefig(path+place+'/'+place+'_'+m+'.png',
                   #bbox_inches='tight', dpi=85)
        plt.savefig('/data/FSHML/Tortues/LMTL-WP4/mercatorglorys2v4_1998_2015/run/actif_vgpm_cor_param_2002_2008_GLO2_F0_alpha_adjusted/animation/actif_vgpm_cor_param_2002_2008_GLO2_F0_alpha_adjusted_'+m+'.png',
        bbox_inches='tight', dpi=85)
        plt.close()


def main() :

    """ Main function where parameters for the plotting and files paths 
    can be entered """
    #lecture des arguments

    path='/'
    file_path=path+sys.argv[1]+'.nc'
    file_name=sys.argv[1] 
    
    # Lecture du fichier d'entrée
    nc_dico=read_nc(file_path,['traj_pp'])
    nsteps,turtles=np.shape(nc_dico['traj_pp'])

    #Age au départ
    start_age = 0
    #Date de départ
    start_day = 0
    end_day = 6570
    #Température optimale pour le calcul d'habitat
    To = 24.
    # Temps de léthargie (nombre de jours maximum dans Tw<Tmin)
    lethargy=10.
    # Variables pour le calcul de l'habitat thermique et alimentaire
    coef_SMR=5.
    Fa = 55.

    place = file_name #'Yalimapo_summer_2002_5000indiv_18y_PPmax80_alpha2.86e6_vscale1.2_speedrecord_yesrebond_init_pos_dbleG_10d_85dpi'
    tracer = "PP"
    mode = 'active'
    species = 'leatherback'


    ##Zone de plot
    #Pacifique

    #Atlantique
    lonmin = 255.
    lonmax = 396.
    latmin = -7
    latmax = 62
    # Espacement des méridiens et des parallèles sur les figures
    lon_space = 30
    lat_space = 20
    # Espacement des jours pris en compte pour l'animation 
    h = int(10)
    # Read nc file
    dico = read_nc(file_path, ['traj_lat','traj_lon','init_t', 'traj_time','traj_temp'])
    plot_animation_frames(path,dico,To,lethargy,coef_SMR,Fa,start_day,end_day,turtles,h,[latmin,latmax],[lonmin,lonmax],lat_space,lon_space,place,tracer,mode,species,start_age)  

    
if __name__  ==  '__main__':
    main()
    print("\nTerminé.\n")

