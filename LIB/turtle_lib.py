#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 18 10:46:43 2019

@author: tcandela
"""
# =============================================================================
# IMPORTS
# =============================================================================
import numpy as np
import netCDF_lib as ncl
# =============================================================================
# FUNCTIONS
# =============================================================================
def find_date_death(turtles,temp,To,coef_SMR,lethargy,init_t,days=6570.):
    """
    Returns 1d-array with 1 date per turtle (number of days since 01/01/ystart)
    1e34 if turtle is still alive at the end of the simulation.
    """
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

    date_death = np.array(date_death)
    return date_death


def find_disabled(infile):
    """
    Find disabled turtles.
    Returns 2 arrays:
    - disabled_turtles: 1D-array (size nb_disabled) containing the index of disabled turtles at the end of simulation.
    - disabled_time: 1D-array (size nb_disabled) containing for each turtle the time index where it was disabled.
    """
    print('Calculating disabled turtles...')
    dic = ncl.read_nc(infile, ['active'])
    active = dic['active']
    #
    disabled_turtles = np.where(active[-1] == 0)[0]
    nb_disabled = len(disabled_turtles)
    print('       ', nb_disabled, 'turtles disabled of', active.shape[1])
    #
    disabled_time = np.zeros(nb_disabled, dtype='int32')
    for k in range(nb_disabled):
        t = disabled_turtles[k]
        tmp = list(active[:,t])
        disabled_time[k] = tmp.index(0)
    #
    return disabled_turtles, disabled_time



def latitude_strip(lat, latmin, latmax):
    """ 
    Set lat to nan where turtles are not between latmin and latmax.
    """
    lat = np.where((lat > latmax) | (lat < latmin), np.nan, lat)
    return lat


def insert_nan(array, turtles_idx, time_idx):
    """
    Insert nan to turtles 'turtles_idx' from 'time_idx' to the end.
    """
    for k in range(len(turtles_idx)):
        t = turtles_idx[k]
        time = time_idx[k]
        array[time:,t] = np.nan
    return array





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

def remove_dead_turtles(x,y,temp,t_init,To,lethargy,coef_SMR, days=6571):
    
    #Age et Tmin fonction de lage
    age_year = np.arange(days-1)/365.
    Tmin = To - coef_SMR*0.21*np.sqrt(0.000214*(1.43*(1-np.exp(-0.226*(age_year+0.17)))*100)**2.86)

    #On vire la première donnée de temp
    temp[0,:] = temp[1,:]
    dayDeath = []    
    status=[]
    k=0
    a=0
    d=0
    
    #l'animal considéré a t-il expérimenté des températures dangereuses ?
    for i in range(np.shape(temp)[1]):
        temp[:days-1,i]=temp[:days-1,i]-Tmin
        temp_serie = temp
        c=1
        #When T>Tmin, tdiff>0 ==> NaN for this day
        #When T<Tmin, tdiff<0 ==> Number of days in row where T<Tmin
        for j in range(0,days):
            if temp[j,i]>0:
                temp_serie[j,i] = np.nan
                c=1
            else : 
                temp_serie[j,i] = c
                c=c+1

        #Search the first day when T<Tmin (after a defined period of time, e.g 1 or 10 days)
        if np.shape(np.where(temp_serie[:,i]==lethargy))[1]>0:
            status.append(0)
            status_temp = 0
            dayDeath = (np.min(np.where(temp_serie[:,i]==lethargy)))
            x[dayDeath:,i] = 0
            y[dayDeath:,i] = 0
            d=d+1
            
        #Ou s'il est resté dans des eaux "chaudes"
        else :
            status.append(1)
            status_temp=1
            
            a=a+1

        k=k+1

    return x, y

def identify_dead_indiv(coef_SMR,temp_area_date,init_t,lethargy=10,days=6570):

    age_year = np.arange(days-1)/365.

    Tmin = 24 - coef_SMR*0.21*np.sqrt(0.000214*(1.43*(1-np.exp(-0.226*(age_year+0.17)))*100)**2.86)

    #print 'SHAPE Tmin',np.shape(Tmin)
    temp_area_date[0,:] = temp_area_date[1,:]

    index_dead = []
    index_alive = []
    date_death = []
    a = 0
    d = 0
    for i in np.arange(np.shape(temp_area_date[0,:])[0]):

        c=1

        temp_diff = np.ones([np.shape(temp_area_date[:,0])[0]])*np.nan
        temp_diff[int(np.ceil(init_t[i])):int(days+np.ceil(init_t[i])-1)]=temp_area_date[int(np.ceil(init_t[i])):int(days+np.ceil(init_t[i])-1),i]-Tmin
        temp_serie = temp_diff
        #When T>Tmin, tdiff>0 ==> NaN for this day
        #When T<Tmin, tdiff<0 ==> Number of days in row where T<Tmin
        for j in np.arange(np.shape(temp_diff)[0]):
            if temp_diff[j]>0 or np.isnan(temp_diff[j])==True:
                temp_serie[j] = np.nan
                c=1
            else :
                #print "Tdiff",k,i,c
                temp_serie[j] = c
                c=c+1

        if np.shape(np.where(temp_serie[:]==lethargy))[1]>0:
            index_dead.append(i)
            #print np.min(np.where(temp_serie[:]==lethargy))
            date_death.append(np.min(np.where(temp_serie[:]==lethargy)))
            d = d+1
        #Ou s'il est resté dans des eaux "chaudes"
        else :
            index_alive.append(i)
            date_death.append(1e34)
            a = a+1

    return index_dead, index_alive, date_death

def num_dead(duration,tseuil):
    sum_dead = 0
    death = []

    for t in np.arange(duration-1) :

        ind_dead = np.where(tseuil[:]==int(t))
        sum_dead = sum_dead + np.shape(ind_dead)[1]
        death.append(sum_dead)
    
    return np.array(death)