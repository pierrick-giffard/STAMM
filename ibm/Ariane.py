#!/bin/env python
# -*- coding: utf-8 -*-


"""@package docstring
A useful tool to get, print(and manipulate Ariane output.)
Used for leatherback turtle drift simulation during summer 2009 (A. Reveiller),
green turtle during summer 2010 (G. Jacob) and
leatherback and green turtles Feb-June 2011 (C. Meetoo)

Note: Most of these functions depend on data input and output frequencies, here
we assume daily input and ouputs. If you change those values (through Ariane 
namelist settings), some functions might give wrong results you modified
them accordingly.
"""

from mpl_toolkits.basemap import Basemap
import numpy as np
from pylab import *
from math import *
import matplotlib.pyplot as plt
from matplotlib import mpl
import Scientific.IO.NetCDF as nc
import matplotlib.cm as cm
import os, sys
import matplotlib.colors
from MA import *
import matplotlib.patches as patches
import matplotlib.path as path


__author__ = "A. Reveillere, G. Jacob, C. Meetoo"
__version__ = "2011.06"


# Color vector used for multi-color plots
color=['b','y','r','g', 'c', 'k']
colormany= ['#FFCC00', #light orange
            '#993300', #brown
            '#CC0066', #dark pink
            '#FF0000', #red
            '#006600', #green
            '#0000FF', #blue
            '#33FF33', #green fluo
            '#3300CC', #violet
            '#FF00FF', #fuchsia
            '#663300', #dark brown
            '#FF9933', #orange
            '#FF99CC', #light pink
            '#000066', #dark blue
            '#660066', #magenta
            '#FFFF99', #yellow
            '#00CCFF', #turquoise
            '#003300', #dark green
            '#FF6600', #reddish orange
            '#999900' ] #kaki



class data:
    """
    class data contains all results from Ariane simulation
    contains also all routines for data analysis
    """

    def __init__(self,filename=str()):
        """
        called with the path to the netcdf Ariane output file as an argument
        load all data

        filename : string, path to netcdf file
        """

        # file opening
        try :
            infile = nc.NetCDFFile(filename,'r')
        except :
            print(str(filename)+"\n\n\n file does not exist\n\n")
            exit()

        # loading data
        lon=infile.variables['traj_lon'].getValue()
        print('   => longitude loaded',lon)
        lat=infile.variables['traj_lat'].getValue()
        print('   => latitude loaded',lat)
        
        try : # key_alltracers = .TRUE. in namelist
            temp=infile.variables['traj_temp'].getValue()
            print('   => temperature loaded')
            geotracer=infile.variables['traj_dens'].getValue()
            print('   => geotracer loaded')
            pp=infile.variables['traj_salt'].getValue()
            print('   => primary production loaded')
        except :
            temp=np.zeros(lat.shape)
            geotracer=np.zeros(lat.shape)
            pp=np.zeros(lat.shape)
        init_t=infile.variables['init_t'].getValue()
        init_x=infile.variables['init_x'].getValue()
        init_y=infile.variables['init_y'].getValue()
        #final_age=infile.variables['final_age'].getValue()

        # closing file
        infile.close()

        # assigning data (variables from class data)
        self.nb_output,self.ntraj=lat.shape
        self.nb_year=(self.nb_output-1)/365. # for daily outputs
        self.nb_season=int(init_t[-1]-0.5)/365+1 ## 1j (daily inputs)

        self.mask=(np.array(geotracer)>1E10)
        self.init_t=np.array(init_t,dtype=np.float32)
        self.init_x=np.array(init_x,dtype=np.float32)
        self.init_y=np.array(init_y,dtype=np.float32)
        #self.final_age=np.array(final_age,dtype=np.float32) #### CM ###

        self.lat=var(lat,self.mask,np.float32,'latitude',-90,90)
        self.lat0=np.mean(lat[0,:])
        self.lon=var(lon,self.mask,np.float32,'longitude',-180,180)
        self.lon0=np.mean(lon[0,:])
        self.temp=var(temp,self.mask,np.float32,'temperature',0,35)
        self.geotracer=var(geotracer,self.mask,np.int,'geotracer',1,np.max(geotracer))
        self.pp=var(pp,self.mask|(np.array(pp)<0),np.float32,'primary production',0,80)
        self.filename=filename

        print('Zone de depart : ' + str(self.lat0) + ', ' + str(self.lon0))


    def mask_zone(self,zone):
        """
        Returns a mask for all data
        mask=False only where geotracer==zone

        zone: integer between 1 and 4
        """
        mask = (self.geotracer.table!=zone) | self.mask
        return mask


    def list_traj_in_zone_at(self,zone,year=1.,traj='all'):
        """
        Returns a list of Ariane trajectories that are in the specified zone at age=year

        zone: integer between 1 and 4
        year: float, age of the particles for the test, default is 1
        """
        if traj=='all':
            traj=range(self.ntraj)


        ind=int(year*(self.nb_output-1)/self.nb_year)
        l=np.ma.where(np.array(self.geotracer.table[ind,traj]+0.5,dtype=int)==zone)
        return np.array(l[0])


    def list_traj_in_box_at(self,lat1,lat2,lon1,lon2,year=1.):
        """
        Returns the list of indexes for Ariane trajectories that are in a zone at age=year

        lat1, lat2: floats between -90 and 90
        lon1, lon2: floats between -180 and 180
        year : float, age of the particles for the test, default is 1
        """
        ind=int(year*(self.nb_output-1)/self.nb_year)

        traj=np.where((self.lat.table[ind,:] > min(lat1, lat2))&
                      (self.lat.table[ind,:] < max(lat1, lat2))&
                      (self.lon.table[ind,:] > min(lon1, lon2))&
                      (self.lon.table[ind,:] < max(lon1, lon2)))
        return np.array(traj[0])


    def list_traj_in_box(self,lat1,lat2,lon1,lon2,inverse=False):
        """
        Returns the list of indexes for Ariane trajectories that have been through a zone
        and the mean age of particle positions in this zone

        lat1, lat2: floats between -90 and 90
        lon1, lon2: floats between -180 and 180
        """

        aux=np.where((np.array(self.lat.table) > min(lat1, lat2))&
                     (np.array(self.lat.table) < max(lat1, lat2))&
                     (np.array(self.lon.table) > min(lon1, lon2))&
                     (np.array(self.lon.table) < max(lon1, lon2)))
        #print('aux', aux)
        age_mean=np.mean(aux[0])
        traj=np.unique(aux[1])
        #print('traj', traj)

        if inverse==False:
            return np.array(traj) #,age_mean peut etre renvoye

        if inverse==True:
            traj = list(traj)
            newtraj = []
            for i in range(self.ntraj):
                if traj.count(i)==0:
                    newtraj.append(i)
            return np.array(newtraj)


    def list_traj_through_zone(self,zone=3,traj='all',tstay=1,time='all',last=False):
        """
        Returns the list of indexes for trajectories that have been through a zone

        tstay: number of days the particle must have been in zone (not consecutively)
        time: number of first days of simulation to be considered (ex: 61 = 2 first months), if 'all' entire simulation
        last: if True, checks if particle in zone at time
        """

        print('  Zone ',zone)
        if traj=='all':
            traj=range(self.ntraj)
        
        # converting traj into int
        traj = [int(i) for i in traj]

        if last==True:
            table = np.array(self.geotracer.table[time,:])
            ind = np.where(table==zone)
            return ind[0]

        # getting zones table (geotracer)
        if time == 'all':
            table = np.array(self.geotracer.table[:,:])
        else:
            table = np.array(self.geotracer.table[:time,:])

        trajex = []
     
        """
        old version ?
        for ind in range(len(table[:,0])) : # first 2 months
            l=np.where(table[ind,:]==zone)
            l=l[0]
           
            if (len(l)>0): # particle has ever been in zone
                for i in range(len(l)):
                    trajex.append(l[i])
        trajex = np.array(trajex)
        trajex = np.unique(trajex)
        """

        # balayage des particules pour voir s'ils ont passe un temps min tstay dans une zone
        for part in range(len(traj)) :
            l=np.where(table[:,traj[part]]==zone)
            l=l[0]
            
           
            if (len(l)>=tstay): # particle has ever been in zone tstay days
                trajex.append(traj[part])

        trajex = np.array(trajex)
        trajex = np.unique(trajex)

        return trajex


    def list_traj_at_month(self, month):
        """
        Returns the list of indexes corresponding to particles launched during 
        a given month, no matter the season.
        
        month: integer, 1 <= month <= 12
        """

        ind=[]
        for s in range(self.nb_season):
            ind += np.where((self.init_t>=(365*s+30*(month-1)))&(self.init_t<(365*s+30*month)))[0].tolist()

        return ind


    def traj_closest_to_point(self,lat,lon):
        """
        Returns the index for the closest Ariane trajectory to a given point
        calculations for distance are to be improved (Dlat**2+Dlon**2)

        lat: float between -90 and 90
        lon: float between -180 and 180
        """
        mat_distance=(self.lat.table-lat)**2+(self.lon.table-lon)**2
        min_distance=np.amin(mat_distance,axis=0)
        ind_traj=np.argmin(min_distance)
        print('la trajectoire d\'indice '+str(ind_traj)+' est la plus proche des coordonnees : \nlat : '+str(lat)+'\nlon : '+str(lon))
        return ind_traj


    def traj_closest_to_traj(self,lat,lon):
        """
        Returns the index for the closest Ariane trajectory to a given trajectory
        calculations for distance are to be improved (Dlat**2+Dlon**2)

        lat: 1D vector of floats between -90 and 90
        lon: 1D vector of floats between -180 and 180
        lat and lon must have the same shape
        """
        if not(lat.shape.__len__()==1 and lon.shape.__len__()==1):
            print("trajectory is not 1D")
            exit()
        N=lat.shape[0]
        distance=np.zeros((N,self.ntraj),dtype=np.float32)
        p=-1
        for i in range(N):
            distance[i,:]=np.amin((self.lat.table[:,:]-lat[i])**2+(self.lon.table[:,:]-lon[i])**2,axis=0)
            if i*100/N!=p : print(str(i*100/N)+' %')
            p=i*100/N
        distance_cum=np.sum(distance,axis=0)
        ind_traj=np.argmin(distance_cum)
        print("la trajectoire la plus proche est la trajectoire numero "+str(ind_traj))
        return ind_traj


    def table_zones_season6(self,year=1):
        """
        6 different zones
        Calculates and displays the number of particle in each geographic zone for different nesting seasons
        return the results as a 2D np array

        year: float, age of the particles for the test, default is 1
        """
        table=np.zeros((7,self.nb_season+1),dtype=np.int)
        if year>self.nb_year or year <0 :
            print("year must be >0 and <= number of years simulated")
            year=self.nb_year
        # counters
        for zone in range(6):
            mask=self.mask_zone(zone+1)
            for season in range(self.nb_season):
                jmin=(season)*self.ntraj/self.nb_season
                jmax=(season+1)*self.ntraj/self.nb_season ## POURQUOI PAS PLUTOT self.nb_output
                table[zone,season]=(1-mask[year*365,jmin:jmax]).sum(0)##POURQUOI PAS SUM(1) pour compter les particules qui sont ds la zone.
        #sums
        for zone in range(6):
            table[zone,self.nb_season]=np.sum(table[zone,0:self.nb_season],0)##0 ==> axis where the sum is taken
        for season in range(self.nb_season):
            table[6,season]=np.sum(table[0:6,season],0)
        table[6,self.nb_season]=np.sum(table[0:6,self.nb_season],0)
        self.hist_zones(table)
        return table
        

    def table_zones_season4(self,year=1):
        """
        4 different zones
        Calculates and displays the number of particle in each geographic zone for different nesting seasons
        return the results as a 2D np array

        year: float, age of the particles for the test, default is 1
        """
        table=np.zeros((5,self.nb_season+1),dtype=np.int)
        if year>self.nb_year or year <0 :
            print("year must be >0 and <= number of years simulated")
            year=self.nb_year
        # counters
        for zone in range(4):
            mask=self.mask_zone(zone+1)
            for season in range(self.nb_season):
                jmin=(season)*self.ntraj/self.nb_season
                jmax=(season+1)*self.ntraj/self.nb_season
                table[zone,season]=(1-mask[year*365,jmin:jmax]).sum(0)
        #sums
        for zone in range(4):
            table[zone,self.nb_season]=np.sum(table[zone,0:self.nb_season],0)
        for season in range(self.nb_season):
            table[4,season]=np.sum(table[0:4,season],0)
        table[4,self.nb_season]=np.sum(table[0:4,self.nb_season],0)
        self.hist_zones(table)
        return table


    def table_zones_age(self,season=1):
        """
        6 different zones
        Calculate and displays the number of particle in each geographic zone for different ages
        return the results as a 2D np array

        season : integer between 1 and the number of seasons, default is 1
        """
        table=np.zeros((7,self.nb_year+1),dtype=np.int)
        if season>self.nb_season or season <0 :
            print("season must be >0 and <= number of seasons simulated")
            season=1
        jmin=(season-1)*self.ntraj/self.nb_season
        jmax=(season)*self.ntraj/self.nb_season

        # counters
        for zone in range(6):
            mask=self.mask_zone(zone+1)
            for age in range(self.nb_year):
                table[zone,age]=(1-mask[(age+1)*365,jmin:jmax]).sum(0)

        #sums
        for zone in range(6):
            table[zone,self.nb_year]=np.sum(table[zone,0:self.nb_year],0)
        for age in range(self.nb_year):
            table[6,age]=np.sum(table[0:6,age],0)
        table[6,self.nb_year]=np.sum(table[0:6,self.nb_year],0)
        print(table)
        self.hist_zones(table)
        return table


    def table_initial_zone_time_departure(self, traj=None, mintime=0, plot=False):
        """
        Calculates for each trajectory the time needed to leave it's spawning
        zone (given by geotracer[0,:]).

        mintime : minimum number of consecutive days that a particle must
                  spend outside it's spawning zone to be considered 'gone'
        """
        if traj is None:
            traj = range(self.ntraj)
        ntraj = np.size(traj)

        tracer = np.array(self.geotracer.table[:, traj])
        ind = np.zeros(ntraj)

        for tnum in range(ntraj):
            ind_zone = np.where(tracer[:,tnum] == tracer[0,tnum])[0] # indices des positions dans la zones
            hors_zone = np.where(ind_zone[1:] - ind_zone[:-1] - 1 > mintime)[0] # ecart entre ces indices = 1+temps passe hors de la zone
            if hors_zone.shape[0]==0: # si la trajectoire ne sort pas, retourne le temps max
                ind[tnum] = ind_zone[-1]+1
            else:
                ind[tnum] = np.min(ind_zone[hors_zone])+1

        if plot==True:
            fig = plt.figure()
            n, bins, patches = plt.hist(ind, bins=range(0, tracer.shape[0]+30, 30), facecolor='green', alpha=0.75)
            plt.xlabel('Temps de residence (jours)')
            plt.ylabel('Nombre de trajectoires (total = '+str(ntraj)+')')
            plt.grid(True)
            plt.xticks(bins[range(0, len(bins), 2)])
            return fig, ind, n

        return ind, n


    def print_zones_age(self,season=1):
        """
        4 different zones
        Displays the proportion of particle in each zone against time (age of particles in year)

        season: integer between 1 and the number of seasons, default is 1
        """
        if season>self.nb_season or season <=0 :
            print("season must be >0 and <= number of seasons simulated")
            season=1
        jmin=(season-1)*self.ntraj/self.nb_season
        jmax=(season)*self.ntraj/self.nb_season

        plt.figure(figsize=(8, 6), dpi=100, facecolor='w', edgecolor='w')

        vect_zone=np.zeros(self.nb_output,dtype=float)
        vect_zone_old=np.zeros(self.nb_output,dtype=float)
        tt=np.arange(self.nb_output,dtype=float)

        for zone in range(4):
            vect_zone_old[:]=vect_zone[:]
            mask=self.mask_zone(zone+1)
            vect_zone+=np.array((1-mask[:,jmin:jmax]).sum(1),dtype=float)/(jmax-jmin)
            plt.fill_between(tt,vect_zone_old,vect_zone,color=color[zone],lw=0)

        plt.ylim(0,1)
        plt.yticks(0.2*np.arange(6),('0 %','20 %','40 %','60 %','80 %','100 %'))
        plt.xlim(0,self.nb_output)
        plt.xticks(365*np.arange(self.nb_year),np.arange(self.nb_year))
        plt.grid(lw=2,ls='--')
        plt.show()

    def cold_selection(self,traj=None, trajmin = None, trajmax=None, age=None, agemin=None, agemax=None, cold=0, Tmin_adult=15.0,Tmin_hatch = 24.0, tolerance=False, traj_mortes=False):
        """
        CM
        Returns traj vector according to cold param

        traj: trajectories to be displayed ; [trajmin trajmax]
        age: restriction on the age of the particles ; [agemin agemax]

        cold: 0 if no cold_death, 1 if cold_death non evol, 2 if cold_death evol
        agemax: in years, age wanted
        Tmin_adult: if cold=1 Tmin ; if cold=2 Tmin_adult
        tolerance: if cold, number of days tolerated
        """

        if cold==0:
            print('No cold death')
            if traj is None:
                if trajmin is None:
                    trajmin = 0
                if trajmax is None:
                    trajmax = self.ntraj
                traj = np.arange(trajmin, trajmax) # vecteur [0, 1 ... , 2999]
            if age is None:
                if agemin is None:
                    agemin = 0
                if agemax is None:
                    agemax = self.nb_output
                age = np.arange(agemin, agemax)##agemin=0 && agemax=nb de jours de la simulation.
            cold_death = [] ## cold_death est définie comme une liste vide.

        elif cold==1:
            print('Cold death, fixed temperature =', Tmin_adult)
            print('   * Import du vecteur cold_death')
            if traj is None:
                traj = np.arange(0,self.ntraj)
                cold_death = self.noevol_temperature(plot=False, Tmin=Tmin_adult,tolerance=tolerance)
            else:
                cold_death = self.noevol_temperature(plot=False, Tmin=Tmin_adult, traj=traj, tolerance=tolerance) # recup vecteur contenu jour de mort
            if agemax is None:
                agemax = 1
            ind = np.where(cold_death >= agemax*365-1) # recup indices des part encore vivantes a agemax
            ind = ind[0]
            trajex = np.array(traj)
            traj = trajex[ind]
            age = np.arange(0,agemax*365)
            print('cold_death', cold_death)

        elif cold==2:
            print('Cold death, evolving temperature, Tmin_adult =', Tmin_adult ## Temperature minimum supportée par les adultes.)
            print('   * Import du vecteur cold_death')
            if traj is None:
                traj = np.arange(self.ntraj)## Indices des particules
                cold_death = self.evol_temperature(plot=False, Tmin_adult=Tmin_adult,Tmin_hatch=Tmin_hatch, tolerance=tolerance)
            else:
                cold_death = self.evol_temperature(plot=False, Tmin_adult=Tmin_adult, traj=traj, Tmin_hatch=Tmin_hatch, tolerance=tolerance) # recup vecteur contenu jour de mort
            if  traj_mortes==False:
                ind = np.where(cold_death > agemax*365-1) # recup indices des part encore vivantes a agemax
                ind = ind[0]
                traj = traj[ind]
                age = np.arange(0,agemax*365)
            print('cold_death', cold_death)

        return traj, age, cold_death


    def hist_pp(self, traj=None, age=None, agemin=None, agemax=None, cold=0, Tmin_adult=15.0,tolerance=False,precision=20,titre="", histo=False, month='all',traj_mortes=False):
        """
        Traitement des pp : histogram, integrales pp 

        traj: trajectories to be displayed ; [trajmin trajmax]
        age: restriction on the age of the particles ; [agemin agemax]

        cold: 0 if no cold_death, 1 if cold_death non evol, 2 if cold_death evol
        agemax: in years, age wanted
        Tmin_adult: if cold=1 Tmin ; if cold=2 Tmin_adult
        tolerance: if cold, number of days tolerated

        histo: if True, draw histogram
        precision: precision of histogram (number of subdivisions)
        titre: title of histogram

        month: in [1:12], month concerned, no histogram plotted

        calls self.evol_temperature, self.noevol_temperature, self.traj_month
        """

        # Selection des trajectoires
        traj, age, cold_death = self.cold_selection(traj=traj,age=age,agemin=agemin,agemax=agemax,cold=cold,Tmin_adult=Tmin_adult,Tmin_hatch=Tmin_hatch, tolerance=tolerance,traj_mortes=traj_mortes)

        print('traj', traj # survivantes, ok )

        nb_traj = np.size(traj)
        nb_output = np.size(age)

        # Calculate sum of pp for each particule p
        print('COMPT PP')

        compt_pp = np.zeros(nb_traj)
        sum_pp = 0                
        for p in range(nb_traj):
            for t in range(nb_output):
                compt_pp[p] = compt_pp[p] + self.pp.table[t,traj[p]]
            sum_pp = sum_pp + compt_pp[p]

        #compt_pp_norm = compt_pp/max(compt_pp)
        compt_pp_norm = compt_pp/sum_pp

        # If selection selon month
        if month != 'all':
            compt_ppmonth = []
            trajmonth = self.traj_month(month)
            traj = list(traj)
            sum_ppmonth = 0 

            for part in range(len(trajmonth)):
                if trajmonth[part] in traj:
                    ind = traj.index(trajmonth[part])
                    compt_ppmonth.append(compt_pp_norm[ind])
                    sum_ppmonth = sum_ppmonth + compt_pp_norm[ind]

            print('compt_pp_normalized month', len(compt_ppmonth)#, compt_ppmonth)
            print('sum_pp_month (ratio)', sum_ppmonth                        )
            return sum_ppmonth

        # If plot histogram
        if histo==True: 

            pp_pas = (max(compt_pp)-min(compt_pp))/precision # pas de l'histogramme

            # axis
            axis = np.zeros(precision+1)
            axis[0] = int(min(compt_pp))
            for i in range(0,precision-1):
                born_sup = min(compt_pp) + (i+1)*pp_pas
                axis[i+1] = int(born_sup)
            axis[-1] = int(max(compt_pp))+1
            print('axis', axis)

            # table
            table = np.zeros(precision)
            for i in range(0,precision-1):
                born_inf = min(compt_pp) + i*pp_pas
                born_sup = min(compt_pp) + (i+1)*pp_pas
                ind_sup = np.where(np.array(compt_pp) < born_sup+0.1)
                compt2 = compt_pp[ind_sup[0]]
                ind_inf = np.where(np.array(compt2) >= born_inf)       
                table[i] = len(ind_inf[0])
            
            table = [int(table[i]) for i in range(len(table))]
            table = np.array(table)
            print('table histo', table)

            # histogram
            n, bins, patches = plt.hist(compt_pp, precision, facecolor='#FFCC00')
            plt.xticks(axis, rotation=45)
            plt.xlim(min(compt_pp)-1000,max(compt_pp+1000))
            plt.ylabel('Number of particles')
            plt.xlabel('PP')
            plt.title(titre)
            plt.show()


    def hist_temp(self, traj=None, age=None, agemin=None, agemax=None, cold=0, Tmin_adult=15.0, tolerance=False, histo=False, precision=20, titre="", month='all', traj_mortes=False):
        """
        Traitement des temperatures de tortues survivantes : histogram, integrales temp 

        traj: trajectories to be displayed ; [trajmin trajmax]
        age: restriction on the age of the particles ; [agemin agemax]

        cold: 0 if no cold_death, 1 if cold_death non evol, 2 if cold_death evol
        agemax: in years, age wanted
        Tmin_adult: if cold=1 Tmin ; if cold=2 Tmin_adult
        tolerance: if cold, number of days tolerated

        histo: if True, draw histogram
        precision: precision of histogram (number of subdivisions)
        titre: title of histogram

        month: in [1:12], month concerned

        calls self.evol_temperature, self.noevol_temperature, self.traj_month
        """
        
        
        # Selection des trajectoires
        traj, age, cold_death = self.cold_selection(traj=traj,age=age,agemin=agemin,agemax=agemax,cold=cold,Tmin_adult=Tmin_adult,Tmin_hatch=Tmin_hatch, tolerance=tolerance, traj_mortes=traj_mortes)

        print('traj', traj # survivantes, ok )

        nb_traj = np.size(traj)
        nb_output = np.size(age)

        # Calculate sum of temp for each particule p
        print('COMPT TEMP')
        compt_temp = np.zeros(nb_traj)
        sum_temp = 0                
        for p in range(nb_traj):
            #print(' p ', p)
            for t in range(nb_output):
                #print('    t ', t)
                compt_temp[p] = compt_temp[p] + self.temp.table[t,traj[p]]
            sum_temp = sum_temp + compt_temp[p]

        print('compt_temp', len(compt_temp)#, compt_temp)
        #compt_temp_norm = compt_temp/max(compt_temp)
        compt_temp_norm = compt_temp/365
        print('compt_temp_normalized', len(compt_temp_norm), compt_temp_norm)

        if month != 'all':

            print('Month ', month)

            compt_tempmonth = []
            trajmonth = self.traj_month(month)
            print('trajmonth', len(trajmonth))#, trajmonth)
            traj = list(traj)
            sum_tempmonth = 0 

            for part in range(len(trajmonth)):
                if trajmonth[part] in traj:
                    ind = traj.index(trajmonth[part])
                    compt_tempmonth.atempend(compt_temp_norm[ind])
                    sum_tempmonth = sum_tempmonth + compt_temp_norm[ind]

            print('compt_temp_normalized month', compt_tempmonth)
            print('sum_temp_month (ratio)', sum_tempmonth         )
                       
            return sum_tempmonth

        if histo==True: 

            temp_pas = (max(compt_temp_norm)-min(compt_temp_norm))/precision # pas de l'histogramme

            # axis
            axis = np.zeros(precision+1)
            axis[0] = int(min(compt_temp_norm))
            for i in range(0,precision-1):
                born_sup = min(compt_temp_norm) + (i+1)*temp_pas
                axis[i+1] = int(born_sup)
            axis[-1] = int(max(compt_temp_norm))+1
            print('axis', axis)

            # table
            table = np.zeros(precision)
            for i in range(0,precision-1):
                born_inf = min(compt_temp_norm) + i*temp_pas
                born_sup = min(compt_temp_norm) + (i+1)*temp_pas
                ind_sup = np.where(np.array(compt_temp_norm) < born_sup+0.1)
                compt2 = compt_temp_norm[ind_sup[0]]
                ind_inf = np.where(np.array(compt2) >= born_inf)       
                table[i] = len(ind_inf[0])
            
            table = [int(table[i]) for i in range(len(table))]
            table = np.array(table)
            print('table histo', table)

            # histogram
            n, bins, patches = plt.hist(compt_temp_norm, precision, facecolor='#FFCC00')
            plt.xticks(axis, rotation=45)
            plt.xlim(min(compt_temp_norm)-0.5,max(compt_temp_norm)+0.5)
            plt.ylabel('Number of particles')
            plt.xlabel('Temp')
            plt.title(titre)
            plt.show()

    def plot_temp_age(self, traj=None, titre='Profil de temperature des tortues', cut=False, Tmin_adult = 15.0, Tmin_hatch=24.0, showcurve=False, tolerance=1, mean=False):
        """
        Plots temperature profiles for some trajectories
        Ex: approaching california, not entering the Tomini Bay, below 30degS or above 30degN

        cut:  - if True, in the beginning all trajectories, and as time goes by,
                trajectories which have a too low temperature (below tolerated = f(Tmin_adult)) are removed...
                In the end: only survivors left
                returns number of survivors
              - if False, all trajectories are plotted

        tolerance: number of days of tolerance
        showcurve : if True, shows death curve at Tmin_adult               
        
        mean = - if False, temperature profiles as lines, fouillis
               - if True, calculates mean profile and displays colored +2STD and -2STD
        """
        
        if traj==None:
            traj = np.arange(0, self.ntraj) # vecteur [0, 1 ... , 4999]
        #print(traj)

        #Temperature loading
        t = self.temp.table[:,traj]
        #print(t)
         
        #Variables : converting age into SCL
        Linf = 1.43
        L0 = 0.053 
        k = 0.226
        age0 = -0.17

        #Temperature en fonction de SCL - necessary variables
        a = (Tmin_hatch - Tmin_adult)/(L0-Linf)
        b = Tmin_hatch - a*L0

        #Local variables
        temp = 0.0*np.zeros(t.shape[0]) # critical temperatures depending on age

        #Critical temperatures for different ages
        #temp[10] => age = 9 days
        print('  - Calcul des temperatures tolerees')
        ages = []

        for time in range(len(temp)): #nb de lignes => nb de jours
            age = float(time) / 365.0
            ages.append(age)
            SCL = Linf*(1-exp(-k*(age-age0)))
            temp[time] = a*SCL + b

        if cut==False :
            print('  - Plot, cut=False')
            if mean==False:
                print('   - Mean = False => courbes temperature')
                for n in range(len(traj)): 
                    ### Courbes des temperatures en fonction de l'age 
                    print(n, t[-1,n]      )
                    p = plt.plot(ages, t[:,n], '-', lw=0.5)
                    
                if (showcurve==True):
                    print('   - Plot curve = True')
                    print('?    ', len(ages))
                    print('?    ', len(temp))
                    p4 = plt.plot(ages, temp, '--',lw=8, c='k')
                    #legend((p,p4),('Temperature for each traj', 'Death curve'), loc='best') 
            else:
                print('   - Mean = True => mean temperature ; - 2STD ; +2STD')
                p3stdtemp = np.zeros(len(ages))
                m3stdtemp = np.zeros(len(ages))
                meantemp = np.zeros(len(ages))
                for ti in range(len(ages)):
                    meantemp[ti]=np.mean(t[ti,:])
                    p3stdtemp[ti]=meantemp[ti]+2*np.std(t[ti,:])
                    m3stdtemp[ti]=meantemp[ti]-2*np.std(t[ti,:])
                p1 = plt.plot(ages,meantemp,'-', lw =4, c='k')
                p2 = plt.plot(ages,p3stdtemp,':', lw =1, c='k')
                p3 = plt.plot(ages,m3stdtemp,':', lw =1, c='k')
                plt.fill_between(ages,m3stdtemp,p3stdtemp,color='#CCCCCC',alpha=0.6)
                if (showcurve==True):
                    print('   - Plot curve = True')
                    print('?    ', len(ages))
                    print('?    ', len(temp))
                    p4 = plt.plot(ages, temp, '--',lw=4, c='r')
                #legend((p1,p2,p3,p4), ('Mean', 'Mean + 2 Std', 'Mean - 2 Std', 'Death curve')) 
                #legend((p1,p4), ('Mean SST', 'Minimum tolerated temperature')) 

            xlabel('Age (years)')
            ylabel('Temperature ( C)')
            xlim(0,6)
            ylim(0,40)
            title(titre)
            grid(True)

        else:
            print('  - Plot, cut=True')
            dead =[]
            for n in range(len(traj)):   
                ### Courbes des temperatures en fonction de l'age
                #print('Particle ', n)
                tempn = t[:,n]

                ind = np.where(np.array(tempn) < np.array(temp))

                ind = ind[0]

                if (len(ind) >= tolerance):
                    for i in range(len(ind)-tolerance):
                        first = len(tempn)
                        if (ind[i+tolerance-1]==ind[i]+tolerance-1):
                            dead.append(n)
                            first = ind[i] # premier indice ou la temperature est trop faible
                            break
                # a partir de cet indice la temperature n'est plus plottee
                else:
                    first = len(tempn)

                #print('ind',len(ind[0]), ind)
                #print('temp',len(t[ind,n]))
                #print('tempn', len(temp[ind])               )

                tempn = tempn[0:first]
                ages_first = ages[0:first]

                plt.plot(ages_first, tempn, '-', lw=0.5)
            
            if (showcurve==True):
                print('   - Plot curve = True')
                print('?    ', len(ages))
                print('?    ', len(temp))
                plt.plot(ages, temp, '--',lw=8, c='k')
                
            xlabel('Age (years)')
            ylabel('Temperature ( C)')
            title(titre)
            grid(True)
            print(len(dead), ' dead particles')
            #print(dead)
            left = len(traj) - len(dead)
            return left

    def plot_lat_time(self,traj='all',titre='Bla'):
        age=range(self.nb_output)
        plot(age, self.lat.table[age, :][:, traj])
        xlabel('Age (day)')
        ylabel('Latitude (deg)')
        title(titre)
        #show()

    def plot_temp_time(self,traj='all',titre='Bla'):
        age=range(self.nb_output)
        plot(age, self.temp.table[age, :][:, traj])
        xlabel('Age (day)')
        ylabel('Temperature (deg)')
        title(titre)
        #show()
            

    def traj_month(self, month):
        """
        Returns vector containing particle indices, given a month in 1:12
        """

        alltraj=np.arange(0,self.ntraj)

        if month==1: #jan
            month=np.array([])
            for y in range(0,6):
                month = np.hstack((month,alltraj[0+3650*y:309+3650*y+1]))
        elif month==2: #fev
            month=np.array([])
            for y in range(0,6):
                month = np.hstack((month,alltraj[310+3650*y:589+3650*y+1]))
        elif month==3: #mar
            month=np.array([])
            for y in range(0,6):
                month = np.hstack((month,alltraj[590+3650*y:899+3650*y+1]))
        elif month==4: #apr
            month=np.array([])
            for y in range(0,6):
                month = np.hstack((month,alltraj[900+3650*y:1199+3650*y+1]))
        elif month==5: #may
            month=np.array([])
            for y in range(0,6):
                month = np.hstack((month,alltraj[1200+3650*y:1509+3650*y+1]))
        elif month==6: #jun
            month=np.array([])
            for y in range(0,6):
                month = np.hstack((month,alltraj[1510+3650*y:1809+3650*y+1]))
        elif month==7: #jul
            month=np.array([])
            for y in range(0,6):
                month = np.hstack((month,alltraj[1810+3650*y:2119+3650*y+1]))
        elif month==8: #aug
            month=np.array([])
            for y in range(0,6):
                month = np.hstack((month,alltraj[2120+3650*y:2429+3650*y+1]))
        elif month==9: #sep
            month=np.array([])
            for y in range(0,6):
                month = np.hstack((month,alltraj[2430+3650*y:2729+3650*y+1]))
        elif month==10: #oct
            month=np.array([])
            for y in range(0,6):
                month = np.hstack((month,alltraj[2730+3650*y:3039+3650*y+1]))
        elif month==11: #nov
            month=np.array([])
            for y in range(0,6):
                month = np.hstack((month,alltraj[3040+3650*y:3339+3650*y+1]))
        elif month==12: #dec
            month=np.array([])
            for y in range(0,6):
                month = np.hstack((month,alltraj[3340+3650*y:3649+3650*y+1]))

        # converting into int
        month = [int(i) for i in month]

        return month


    def hist_zones(self,table):
        """
        Display routine, creates histogram from table_zone result

        called by: table_zones_season, table_zones_age
        should not be called directly
        """
        plt.figure(figsize=(6, 6), dpi=100, facecolor='w', edgecolor='w')

        plt.grid()
        width=0.15

        ind=np.arange(table.shape[1]-1)

        
        p1=plt.bar(ind+0*width,table[0,:-1],width,color=color[0])
        p2=plt.bar(ind+1*width,table[1,:-1],width,color=color[1])
        p3=plt.bar(ind+2*width,table[2,:-1],width,color=color[2])
        p4=plt.bar(ind+3*width,table[3,:-1],width,color=color[3])

        self.autolabel(p1)
        self.autolabel(p2)
        self.autolabel(p3)
        self.autolabel(p4)
        plt.xticks(ind+2*width,[])
        plt.yticks(plt.yticks()[0],[])
        plt.ylim(0,plt.ylim()[1]*1.5)
        plt.xlim(0,ind[ind.shape[0]-1]+4*width)

        plt.legend((p1[0],p2[0],p3[0],p4[0]),('(IC) Indonesia and China','(NP) North Pacific','(SP) South Pacific','(IO) Indian Ocean'),loc='best')
        plt.show()
        
    def hist_zones_nourrissageIO(self,title='Blabla',percent=False,loca='upper center',month='all', time=61, hist=True,zone=0):
        """
        Display histogram of particles that have been in different food zones (used for the Indian Ocean ONLY)

        title: histogram title
        percent: if False, y-axis in number of Turtles. Else, number of total turtles, gives percentage
        loca: localization of legend on histogram
        month: if given number in [1:12], specific month and all years
        time: number of days considered (ex: time=61 => first 2 months of simulation)
        
        Food zones (18):
        11 = Mauritius
        16 = Reunion
        14 = Rodrigues
        25 = Tromelin
        7 = Seychelles
        4 = Aldabra
        2 = Mayotte
        8 = Comoros
        24 = Glorieuses
        10 = Europa
        15 = Juan de Nova
        23 = Madagascar East
        27 = Madagascar SW
        18 = Madagascar NW
        1 = South Africa
        29 = Mozambique
        22 = Tanzania/Kenya
        26 = Somalia

        calls self.traj_month, self.list_traj_through_zone

        """
        if month!='all':
            month = self.traj_month(month)

            print('No of particles', len(month))

        # testing if particles in month went through zone
        if zone == 0:
            traj11 = self.list_traj_through_zone(zone=11,traj=month,time=time)
            traj16 = self.list_traj_through_zone(zone=16,traj=month,time=time)
            traj14 = self.list_traj_through_zone(zone=14,traj=month,time=time)
            traj25 = self.list_traj_through_zone(zone=25,traj=month,time=time)
            traj07 = self.list_traj_through_zone(zone=07,traj=month,time=time)
            traj04 = self.list_traj_through_zone(zone=04,traj=month,time=time)
            traj02 = self.list_traj_through_zone(zone=02,traj=month,time=time)
            traj08 = self.list_traj_through_zone(zone=8,traj=month,time=time)
            traj24 = self.list_traj_through_zone(zone=24,traj=month,time=time)
            traj10 = self.list_traj_through_zone(zone=10,traj=month,time=time)
            traj15 = self.list_traj_through_zone(zone=15,traj=month,time=time)
            traj23 = self.list_traj_through_zone(zone=23,traj=month,time=time)
            traj27 = self.list_traj_through_zone(zone=27,traj=month,time=time)
            traj18 = self.list_traj_through_zone(zone=18,traj=month,time=time)
            traj01 = self.list_traj_through_zone(zone=01,traj=month,time=time)
            traj29 = self.list_traj_through_zone(zone=29,traj=month,time=time)
            traj22 = self.list_traj_through_zone(zone=22,traj=month,time=time)
            traj26 = self.list_traj_through_zone(zone=26,traj=month,time=time)
            """
            print('traj11',traj11)
            print('traj25',traj25)
            print('traj07',traj07)
            print('traj04',traj04)
            print('traj10',traj10)
            print('traj15',traj15)
            print('traj23',traj23)
            print('traj27',traj27)
            print('traj18',traj18)
            print('traj01',traj01)
            print('traj29',traj29)
            print('traj22',traj22)
            print('traj26',traj26    )
            """

            # Particles that have never been in any zones
            # alltrajzones = particle indexes that have been in zones
            # nullzone = particle indexes that have never been in any zones
            alltrajzones=np.hstack((traj11,traj16))
            alltrajzones=np.hstack((alltrajzones,traj14))
            alltrajzones=np.hstack((alltrajzones,traj25))
            alltrajzones=np.hstack((alltrajzones,traj07))
            alltrajzones=np.hstack((alltrajzones,traj04))
            alltrajzones=np.hstack((alltrajzones,traj02))
            alltrajzones=np.hstack((alltrajzones,traj08))
            alltrajzones=np.hstack((alltrajzones,traj24))
            alltrajzones=np.hstack((alltrajzones,traj10))
            alltrajzones=np.hstack((alltrajzones,traj15))
            alltrajzones=np.hstack((alltrajzones,traj23))
            alltrajzones=np.hstack((alltrajzones,traj27))
            alltrajzones=np.hstack((alltrajzones,traj18))
            alltrajzones=np.hstack((alltrajzones,traj01))
            alltrajzones=np.hstack((alltrajzones,traj29))
            alltrajzones=np.hstack((alltrajzones,traj22))
            alltrajzones=np.hstack((alltrajzones,traj26))

            print('alltrajzones ok')
            alltrajzones=list(np.unique(alltrajzones))

            if month=='all':
                nullzone=list(arange(self.ntraj))
            else:
                nullzone=list(month)

            
            for i in range(len(alltrajzones)):
                nullzone.remove(alltrajzones[i])
            
            nullzone = np.array(nullzone)
            print('nullzone', len(nullzone))
            print('alltrajzones', len(alltrajzones))

            if hist==False:
                if percent==False:
                    print('Mau', len(traj11))
                    print('Reu', len(traj16))
                    print('Rod', len(traj14))
                    print('Tro', len(traj25))
                    print('Sey', len(traj07))
                    print('Ald', len(traj04))
                    print('May', len(traj02))
                    print('Com', len(traj08))
                    print('Glo', len(traj10))
                    print('Eur', len(traj15))
                    print('Jdn', len(traj23))
                    print('ME', len(traj27))
                    print('MSW', len(traj18))
                    print('MNW', len(traj01))
                    print('SA', len(traj29))
                    print('Moz', len(traj22))
                    print('Ken', len(traj26))
                    print('Som', len(traj14))
                    print('Null', len(nullzone))
                else:
                    print('Mau', 100*len(traj11)/percent)
                    print('Reu', 100*len(traj16)/percent)
                    print('Rod', 100*len(traj14)/percent)
                    print('Tro', 100*len(traj25)/percent)
                    print('Sey', 100*len(traj07)/percent)
                    print('Ald', 100*len(traj04)/percent)
                    print('May', 100*len(traj02)/percent)
                    print('Com', 100*len(traj08)/percent)
                    print('Glo', 100*len(traj10)/percent)
                    print('Eur', 100*len(traj15)/percent)
                    print('Jdn', 100*len(traj23)/percent)
                    print('ME', 100*len(traj27)/percent)
                    print('MSW', 100*len(traj18)/percent)
                    print('MNW', 100*len(traj01)/percent)
                    print('SA', 100*len(traj29)/percent)
                    print('Moz', 100*len(traj22)/percent)
                    print('Ken', 100*len(traj26)/percent)
                    print('Som', 100*len(traj14)/percent)
                    print('Null', 100*len(nullzone)/percent)

            else:
                # STARTING PLOT #
                # opening histogram
                plt.figure(figsize=(12, 8), dpi=100, facecolor='w', edgecolor='w')
                plt.grid()
                plt.title(title)
                plt.xlabel('Zones')
                plt.xticks(np.arange(1.005,1.185,0.01), ('Mau','Reu','Rod','Tro','Sey','Ald','May','Com','Glo','Eur','Jdn','ME','MSW','MNW','SA','Moz','Ken','Som','Null') )
                plt.xlim(0.995,1.195)

                width=0.01
                ind=1

                if (percent==False):
                    print('percent=False')
                    p1=plt.bar(ind+0*width,len(traj11),width,color=colormany[0])
                    p15=plt.bar(ind+1*width,len(traj16),width,color=colormany[14])
                    p16=plt.bar(ind+2*width,len(traj14),width,color=colormany[15])
                    p2=plt.bar(ind+3*width,len(traj25),width,color=colormany[1])
                    p3=plt.bar(ind+4*width,len(traj07),width,color=colormany[2])
                    p4=plt.bar(ind+5*width,len(traj04),width,color=colormany[3])
                    p17=plt.bar(ind+6*width,len(traj02),width,color=colormany[16])
                    p18=plt.bar(ind+7*width,len(traj08),width,color=colormany[17])
                    p19=plt.bar(ind+8*width,len(traj24),width,color=colormany[18])
                    p5=plt.bar(ind+9*width,len(traj10),width,color=colormany[4])
                    p6=plt.bar(ind+10*width,len(traj15),width,color=colormany[5])
                    p7=plt.bar(ind+11*width,len(traj23),width,color=colormany[6])
                    p8=plt.bar(ind+12*width,len(traj27),width,color=colormany[7])
                    p9=plt.bar(ind+13*width,len(traj18),width,color=colormany[8])
                    p10=plt.bar(ind+14*width,len(traj01),width,color=colormany[9])
                    p11=plt.bar(ind+15*width,len(traj29),width,color=colormany[10])
                    p12=plt.bar(ind+16*width,len(traj22),width,color=colormany[11])
                    p13=plt.bar(ind+17*width,len(traj26),width,color=colormany[12])
                    p14=plt.bar(ind+18*width,len(nullzone),width,color=colormany[13])

                    plt.ylabel('Number of turtles')
                    if month=='all':
                        plt.ylim(0,10200)

                else:
                    print('percent=True')
                    p1=plt.bar(ind+0*width,100*len(traj11)/percent,width,color=colormany[0])
                    p15=plt.bar(ind+1*width,100*len(traj16)/percent,width,color=colormany[14])
                    p16=plt.bar(ind+2*width,100*len(traj14)/percent,width,color=colormany[15])
                    p2=plt.bar(ind+3*width,100*len(traj25)/percent,width,color=colormany[1])
                    p3=plt.bar(ind+4*width,100*len(traj07)/percent,width,color=colormany[2])
                    p4=plt.bar(ind+5*width,100*len(traj04)/percent,width,color=colormany[3])
                    p17=plt.bar(ind+6*width,100*len(traj02)/percent,width,color=colormany[16])
                    p18=plt.bar(ind+7*width,100*len(traj08)/percent,width,color=colormany[17])
                    p19=plt.bar(ind+8*width,100*len(traj24)/percent,width,color=colormany[18])
                    p5=plt.bar(ind+9*width,100*len(traj10)/percent,width,color=colormany[4])
                    p6=plt.bar(ind+10*width,100*len(traj15)/percent,width,color=colormany[5])
                    p7=plt.bar(ind+11*width,100*len(traj23)/percent,width,color=colormany[6])
                    p8=plt.bar(ind+12*width,100*len(traj27)/percent,width,color=colormany[7])
                    p9=plt.bar(ind+13*width,100*len(traj18)/percent,width,color=colormany[8])
                    p10=plt.bar(ind+14*width,100*len(traj01)/percent,width,color=colormany[9])
                    p11=plt.bar(ind+15*width,100*len(traj29)/percent,width,color=colormany[10])
                    p12=plt.bar(ind+16*width,100*len(traj22)/percent,width,color=colormany[11])
                    p13=plt.bar(ind+17*width,100*len(traj26)/percent,width,color=colormany[12])
                    p14=plt.bar(ind+18*width,100*len(nullzone)/percent,width,color=colormany[13])

                    plt.ylabel('% of turtles over '+str(int(percent)))
                    plt.ylim(0,110)

                self.autolabel(p1)
                self.autolabel(p2)
                self.autolabel(p3)
                self.autolabel(p4)
                self.autolabel(p5)
                self.autolabel(p6)
                self.autolabel(p7)
                self.autolabel(p8)
                self.autolabel(p9)
                self.autolabel(p10)
                self.autolabel(p11)
                self.autolabel(p12)
                self.autolabel(p13)
                self.autolabel(p14)
                self.autolabel(p15)
                self.autolabel(p16)
                self.autolabel(p17)
                self.autolabel(p18)
                self.autolabel(p19)          

                # legend only if all particles
                if month=='all':
                    plt.legend((p1[0],p15[0],p16[0],p2[0],p3[0],p4[0],p17[0],p18[0],p19[0],p5[0],p6[0],p7[0],p8[0],p9[0],p10[0],p11[0],p12[0],p13[0],p14[0]),('Mauritius', 'Reunion', 'Rodrigues', 'Tromelin', 'Seychelles', 'Aldabra', 'Mayotte', 'Comoros', 'Glorieuses', 'Europa', 'Juan de nova', 'Madagascar E', 'Madagascar SW', 'Madagascar NW', 'South Africa', 'Mozambique', 'Tanzania/Kenya', 'Somalia','Null'),loc=loca)

        else:
            trajzone = self.list_traj_through_zone(zone=zone,traj=month,time=time)
            
            # Compteur jours
            compt_days = np.zeros(len(trajzone))
            print('len(trajzone)', len(trajzone))
            print('trajzone', trajzone)
            table = np.array(self.geotracer.table)
            print('table.shape', table.shape)
            #set_printoptions(threshold=nan) # non truncated array
            print('table', table)

            for t in trajzone:
                print(t)
                tab = np.array(table[:][t])
                print(tab)
                ind = np.where(tab==zone)
                print(ind)
                compt_days[t] = len(ind[0]) 
                

            # axis
            precision = 10
            pas = (max(compt_days) - min(compt_days))/precision
            axis = np.zeros(precision+1)
            axis[0] = int(min(compt_days))
            for i in range(0,precision-1):
                born_sup = min(compt_days) + (i+1)*pas
                axis[i+1] = int(born_sup)
            axis[-1] = int(max(compt_days))+1
            print('axis', axis)

            # table
            table = np.zeros(precision)
            for i in range(0,precision-1):
                born_inf = min(compt_days) + i*pas
                born_sup = min(compt_days) + (i+1)*pas
                ind_sup = np.where(np.array(compt_days) < born_sup+0.1)
                compt2 = compt_days[ind_sup[0]]
                ind_inf = np.where(np.array(compt2) >= born_inf)       
                table[i] = len(ind_inf[0])
            
            table = [int(table[i]) for i in range(len(table))]
            table = np.array(table)
            print('table histo', table)

            # histogram
            n, bins, patches = plt.hist(compt_days, precision, facecolor='#FFCC00')
            plt.xticks(axis, rotation=45)
            plt.xlim(min(compt_days)-20,max(compt_days)+20)
            plt.ylabel('Number of particles')
            plt.xlabel('Number of days')
            plt.title(titre)
            plt.show()

    def hist_zones_oceanIO(self,title='Blabla',percent=False,loca='best',month='all', time=61, hist=True):
        """
        Display histogram of particles that have been in different oceanic zones (used for the Indian Ocean ONLY)

        title: histogram title
        percent: if False, y-axis in number of Turtles. Else, number of total turtles, gives percentage
        loca: localization of legend on histogram
        month: if given number in [1:12], specific month and all years
        time: number of days considered (ex: time=61 => first 2 months of simulation)
        
        Oceanic zones (8):
        19 = Others (Altlantic, Pacific)
        24 = Indian Ocean East
        30 = Indian Ocean South
        17 = Indian Ocean North
        28 = Mozambique South
        3 = Mozambique North
        21 = Mascareignes
        13 = Seychelles

        calls self.traj_month, self.list_traj_through_zone

        """
        if month!='all':
            month = self.traj_month(month)

            print('No of particles', len(month))

        # testing if particles in month went through zone
        traj19 = self.list_traj_through_zone(zone=19,traj=month,time=time)
        traj24 = self.list_traj_through_zone(zone=24,traj=month,time=time)
        traj30 = self.list_traj_through_zone(zone=30,traj=month,time=time)
        traj17 = self.list_traj_through_zone(zone=17,traj=month,time=time)
        traj28 = self.list_traj_through_zone(zone=28,traj=month,time=time)
        traj03 = self.list_traj_through_zone(zone=03,traj=month,time=time)
        traj21 = self.list_traj_through_zone(zone=21,traj=month,time=time)
        traj13 = self.list_traj_through_zone(zone=13,traj=month,time=time)

        if hist==False:
            if percent==False:
                print('Oth', len(traj19))
                print('IOE', len(traj24))
                print('IOS', len(traj30))
                print('ION', len(traj17))
                print('MoS', len(traj28))
                print('MoN', len(traj03))
                print('Mas', len(traj21))
                print('Sey', len(traj13))
 
            else:
                print('Oth', 100*len(traj19)/percent)
                print('IOE', 100*len(traj24)/percent)
                print('IOS', 100*len(traj30)/percent)
                print('ION', 100*len(traj17)/percent)
                print('MoS', 100*len(traj28)/percent)
                print('MoN', 100*len(traj03)/percent)
                print('Mas', 100*len(traj21)/percent)
                print('Sey', 100*len(traj13)/percent)


        else:
            # STARTING PLOT #
            # opening histogram
            plt.figure(figsize=(12, 8), dpi=100, facecolor='w', edgecolor='w')
            plt.grid()
            plt.title(title)
            plt.xlabel('Zones')
            plt.xticks(np.arange(1.005,1.085,0.01), ('IOE', 'IOS', 'ION', 'MoS', 'MoN', 'Mas', 'Sey', 'Others'))
            plt.xlim(0.995,1.095)

            width=0.01
            ind=1

            if (percent==False):
                print('percent=False')
                p24=plt.bar(ind+0*width,len(traj24),width,color=colormany[0])
                p30=plt.bar(ind+1*width,len(traj30),width,color=colormany[14])
                p17=plt.bar(ind+2*width,len(traj17),width,color=colormany[15])
                p28=plt.bar(ind+3*width,len(traj28),width,color=colormany[1])
                p03=plt.bar(ind+4*width,len(traj03),width,color=colormany[2])
                p21=plt.bar(ind+5*width,len(traj21),width,color=colormany[3])
                p13=plt.bar(ind+6*width,len(traj13),width,color=colormany[16])
                p19=plt.bar(ind+7*width,len(traj19),width,color=colormany[17])


                plt.ylabel('Number of turtles')
                if month=='all':
                    plt.ylim(0,10200)

            else:
                print('percent=True')
                p24=plt.bar(ind+0*width,100*len(traj24)/percent,width,color=colormany[0])
                p30=plt.bar(ind+1*width,100*len(traj30)/percent,width,color=colormany[14])
                p17=plt.bar(ind+2*width,100*len(traj17)/percent,width,color=colormany[15])
                p28=plt.bar(ind+3*width,100*len(traj28)/percent,width,color=colormany[1])
                p03=plt.bar(ind+4*width,100*len(traj03)/percent,width,color=colormany[2])
                p21=plt.bar(ind+5*width,100*len(traj21)/percent,width,color=colormany[3])
                p13=plt.bar(ind+6*width,100*len(traj13)/percent,width,color=colormany[16])
                p19=plt.bar(ind+7*width,100*len(traj19)/percent,width,color=colormany[17])


                plt.ylabel('% of turtles over '+str(int(percent)))
                plt.ylim(0,110)

            self.autolabel(p24)
            self.autolabel(p30)
            self.autolabel(p17)
            self.autolabel(p28)
            self.autolabel(p03)
            self.autolabel(p21)
            self.autolabel(p13)
            self.autolabel(p19)
          

            # legend only if all particles
            if month=='all':
                plt.legend((p24[0],p30[0],p17[0],p28[0],p03[0],p21[0],p13[0],p19[0]),('Indian Ocean East', 'Indian Ocean South', 'Indian Ocean North', 'Mozambique South', 'Mozambique North', 'Mascareignes', 'Seychelles', 'Other'),loc=loca)

             
    def hist_6zones_ocean(self,title='Blabla',percent=False,loca='best', time=365, hist=True):
        """
        Display histogram of particles that have been in different oceanic zones (used for the Indian Ocean ONLY)

        title: histogram title
        percent: if False, y-axis in number of Turtles. Else, number of total turtles, gives percentage
        loca: localization of legend on histogram
        month: if given number in [1:12], specific month and all years
        time: number of days considered (ex: time=61 => first 2 months of simulation)
        
        Oceanic zones (1):
        1 = China sea
        2 = Pacific North
        3 = Pacific South
        4 = Indian Ocean
        5 = Indonesia
        6 = Celebes

        calls self.traj_month, self.list_traj_through_zone

        """

        # testing if particles in month went through zone
        traj1 = self.list_traj_through_zone(zone=1,time=time,last=True)
        traj2 = self.list_traj_through_zone(zone=2,time=time,last=True)
        traj3 = self.list_traj_through_zone(zone=3,time=time,last=True)
        traj4 = self.list_traj_through_zone(zone=4,time=time,last=True)
        traj5 = self.list_traj_through_zone(zone=5,time=time,last=True)
        traj6 = self.list_traj_through_zone(zone=6,time=time,last=True)


        if hist==False:
            if percent==False:
                print('Chi', len(traj1))
                print('NP', len(traj2))
                print('SP', len(traj3))
                print('IO', len(traj4))
                print('Ind', len(traj5))
                print('Cel', len(traj6))

 
            else:
                print('Chi', 100*len(traj1)/percent)
                print('NP', 100*len(traj2)/percent)
                print('SP', 100*len(traj3)/percent)
                print('IO', 100*len(traj4)/percent)
                print('Ind', 100*len(traj5)/percent)
                print('Cel', 100*len(traj6)/percent   )

        else:
            # STARTING PLOT #
            # opening histogram
            plt.figure(figsize=(12, 8), dpi=100, facecolor='w', edgecolor='w')
            plt.grid()
            plt.title(title)
            plt.xlabel('Zones')
            plt.xticks(np.arange(1.005,1.065,0.01), ('Chi', 'NP', 'SP', 'IO','Ind', 'Cel'))
            plt.xlim(0.995,1.075)

            width=0.01
            ind=1

            if (percent==False):
                print('percent=False')
                p1=plt.bar(ind+0*width,len(traj1),width,color=colormany[0])
                p2=plt.bar(ind+1*width,len(traj2),width,color=colormany[14])
                p3=plt.bar(ind+2*width,len(traj3),width,color=colormany[15])
                p4=plt.bar(ind+3*width,len(traj4),width,color=colormany[1])
                p5=plt.bar(ind+4*width,len(traj5),width,color=colormany[2])
                p6=plt.bar(ind+5*width,len(traj6),width,color=colormany[3])


                plt.ylabel('Number of turtles')
                #if month=='all':
                    #plt.ylim(0,22400)
                #plt.ylim(0,450)

            else:
                print('percent=True')
                p1=plt.bar(ind+0*width,100*len(traj1)/percent,width,color=colormany[0])
                p2=plt.bar(ind+1*width,100*len(traj2)/percent,width,color=colormany[14])
                p3=plt.bar(ind+2*width,100*len(traj3)/percent,width,color=colormany[15])
                p4=plt.bar(ind+3*width,100*len(traj4)/percent,width,color=colormany[1])
                p5=plt.bar(ind+4*width,100*len(traj5)/percent,width,color=colormany[2])
                p6=plt.bar(ind+5*width,100*len(traj6)/percent,width,color=colormany[3])


                plt.ylabel('% of turtles over '+str(int(percent)))
                plt.ylim(0,110)

            self.autolabel(p1)
            self.autolabel(p2)
            self.autolabel(p3)
            self.autolabel(p4)
            self.autolabel(p5)
            self.autolabel(p6)

            plt.legend((p1[0],p2[0],p3[0],p4[0],p5[0],p6[0]),('(Chi) China','(NP) North Pacific','(SP) South Pacific','(IO) Indian Ocean', '(Ind) Indonesia', '(Cel) Celebes'),loc='best')


    def noevol_temperature(self, traj='all', titre='Number of living turtles over the time', Tmin=18.0, tolerance=False, plot=False, xy=False,many=False):
        """
        Plots all trajectories except those which are dead due to cold waters
        (according to law SCL = f(age) and min_temp = f(SCL))

        traj
        titre
        Tmin = fixed minimal temperature
        tolerance = if number, number of days Tmin tolerated
        plot

        returns vector cold_death containing time of death for each particle
        """

        if traj=='all':
            traj=range(self.ntraj)

        ntraj = len(traj)
        ind = np.zeros(ntraj)

        #Positions and temperature loading
        t = self.temp.table[:, traj]

        cold_death = len(t[:,0])*np.ones(len(traj)) # dead particles
        counter = 0

        print('  - Parcours des particules pour voir leur mortalite')

        if tolerance==False:
            print('tolerance = False')
            for part in range(len(traj)):
                test_temp = np.where(t[:,part] < Tmin)
                if len(test_temp)>0:
                    cold_death[part] = test_temp[0]
        else:
            print('tolerance =', tolerance)
            """
            for part in range(len(traj)):
                for time_ind in range(len(t[:,part])):
                    if (Tmin >= t[time_ind,part]):
                        counter = counter + 1
                    else:
                        counter = 0
                    if counter>tolerance:
                        cold_death[part] = time_ind
                        break
            """
            for part in range(len(traj)):
                ind = np.array(np.where(t[:,part] < Tmin))
                ind = ind[0]
                
                if len(ind)>=tolerance:
                    for i in range(len(ind)-tolerance):
                        if (ind[i+tolerance-1]==ind[i]+tolerance-1):
                            cold_death[part]=ind[i]
                            break
            #print('cold_death', cold_death)
        
        #set_printoptions(threshold=nan)

        if plot==True: 
            print('   - Plot' )
            ### Plots the time of death of each particle
            #plt.plot(traj, cold_death, '-', lw=1)
            #xlabel('Particle number')
            #ylabel('Time of death (day)')
            #title('T = f(SCL)')
            #grid(True)

            ### Courbes du nombre de tortues en vie en fonction de l'age
            compt = np.zeros(len(t[:,0]))
            for time_ind in range(len(t[:,0])):
                for turt in range(len(traj)):
                    if (cold_death[turt] > time_ind):
                        compt[time_ind] = compt[time_ind]+1   

            ### Pourcentage de survie
            #plt.plot(range(len(temp)), compt/len(traj), '-', lw=2)
            #ylabel('Turtles left (ratio)')          
            if many==False:
                plt.plot(range(len(t[:,0])), compt, '-', lw=2, color='#FF6600')
            else:
                if xy=='5eur' or xy=='5tro' or xy=='5glo' or xy=='5may':
                    plt.plot(range(len(t[:,0])), 100*compt/len(traj), '-', lw=2, color=colormany[many])
                    legend('JFMAMJJASOND',loc='best')
                    plt.ylim(0,103)
                    print('% living: ', 100*compt[-1]/len(traj))
                else:            
                    plt.plot(range(len(t[:,0])), compt, '-', lw=2, color=colormany[many])
                    legend('JFMAMJJASOND',loc='best')
            xlabel('Time (days)')
            ylabel('Turtles left')
            title(titre)
            grid(True)
            plt.xlim(0,370)

            if xy=='1':
                1
            elif xy=='2eur':
                plt.ylim(400,1900)
            elif xy=='2tro':
                plt.ylim(1500,1900)    
            elif xy=='2glo':
                plt.ylim(1200,1900)
            elif xy=='2may':
                plt.ylim(1200,1900) 

            elif xy=='3eur':
                plt.ylim(400,1900)
            elif xy=='3tro':
                plt.ylim(1500,1900)    
            elif xy=='3glo':
                plt.ylim(1200,1900)
            elif xy=='3may':
                plt.ylim(1200,1900) 

            elif xy=='4eur':
                plt.ylim(50,320)
            elif xy=='4tro':
                plt.ylim(200,320)    
            elif xy=='4glo':
                plt.ylim(170,320)
            elif xy=='4may':
                plt.ylim(110,320)  
       
        return cold_death


    ##Tmin_hatch et Tmin_adult ==> Valeurs limites 

    def evol_temperature(self, ms=0.05, traj='all', plot=True, titre='Number of living turtles over the time', Tmin_hatch = 24.0, Tmin_adult = 15.0, tolerance=False):
        """
        IMPORTANT : Set parameters for LUTH or GREEN!
        Plots all trajectories except those which are dead due to cold waters
        (according to law SCL = f(age) and min_temp = f(SCL))

        ms: markersize
        """

        if traj=='all':
            traj=range(self.ntraj) ## Indice des particules.

        ntraj = len(traj) ## Nombre de particules.
        ind = np.zeros(ntraj) ##Défini comme un vecteur de 0 de la dimension : nombre de particules.

        #Positions and temperature loading
        t = self.temp.table[:, traj] ## Température loading pour toutes les particules si traj=all.

        #Variables LUTH : converting age into SCL
        ###if tortle_type=='luth':
            ###Linf = 1.43 
            ###k = 0.226
            ###age0 = -0.17
            ###L0 = 0.053

        #Temperature en fonction de SCL - necessary variables
        ####a = (Tmin_hatch - Tmin_adult)/(L0-Linf)
        ####b = Tmin_hatch - a*L0

        
        ####A METTRE EN PARAMETRES DE CLASSE.
        #Variables VERTES : converting age into SCL
        ###if tortle_type=='verte':
            ###Linf = 0.983 
            ###k = 0.074
            ###age0 = -0.707
        
        #Local variables
        temp = 0.0*np.zeros(t.shape[0]) # critical temperatures depending on age
        ##0.0 pour que temp soit en float, equivalent à temp[:]=float(temp[:])???
        ## t.shape[0]===dimension de temps.        
        cold_death = len(temp)*np.ones(len(traj)) # dead particles
        ##Particules mortes à chaque instant, cold_death est initialisé comme un vecteur de dimension : nombre de particules qui contient pour chaque   composante le nombre de jours de la simulation.
        counter = 0

        #Critical temperatures for different ages

        print('  - Calcul des temperatures tolerees')
        
        for time in range(len(temp)): #nb de lignes => nb de jours
            age = float(time) / 365.0 ## age en années.
            #print('age:', age)
            ####SCL = Linf*(1-exp(-k*(age-age0)))
            #print('SCL:', SCL)
            temp[time] = 13.0 + 11.4*exp(-k*(age-age0)) ## Vecteur avec les températures tolérées pour chaque insant.

        print('  - Parcours des particules pour voir leur mortalite')
        if tolerance==False:
            for part in range(len(traj)): ##PARTICULES
                for time_ind in range(len(temp)): ##TIME
                    if (temp[time_ind] > t[time_ind,part]):
                        cold_death[part] = time_ind
                        counter = counter + 1
                        break
        else:
            for part in range(len(traj)):
                for time_ind in range(len(temp)):
                    if (temp[time_ind] > t[time_ind,part]):
                        counter = counter + 1
                    else:
                        counter = 0
                    if counter>tolerance: # tolerance 1 week 7 days
                        cold_death[part] = time_ind
                        break
        
        ##TOLERANCE = Il existe un nombre de particules mortes maximum ==> FAUX
        ##TOLERANCE = Tolerance sur le nombre de jours en zone à température critique avant de mourir.



        #set_printoptions(threshold=nan)
        if plot==True: 
            print('   - Plot' )
            ### Plots the time of death of each particle
            #plt.plot(traj, cold_death, '-', lw=1)
            #xlabel('Particle number')
            #ylabel('Time of death (day)')
            #title('T = f(SCL)')
            #grid(True)

            ### Courbes du nombre de tortues en vie en fonction de l'age
            compt = np.zeros(len(temp))
            print(len(temp))
            for time_ind in range(len(temp)):
                for turt in range(len(traj)):
                    if (cold_death[turt] > time_ind):
                        compt[time_ind] = compt[time_ind]+1   
 
            ### Pourcentage de survie
            #plt.plot(range(len(temp)), compt/len(traj), '-', lw=2)
            #ylabel('Turtles left (ratio)')          
            
            plt.plot(range(len(temp)), compt, '-', lw=2, color='r')
            xlabel('Age (days)')
            ylabel('Turtles left')
            #ylim(0,5000)
            yticks(np.arange(0,1500,100))
            xticks(np.arange(0,2400,200))
            xlim(0,2200)
            title(titre)
            grid(True)
       
        return cold_death
###################################################################################

    def evol_temperature_zones(self, nb_zones = 6, ms=0.05, traject='all', plot=True, titre='Number of living turtles over the time', Tmin_hatch = 24.0, Tmin_adult = 15.0,year=6):
        """
        Plots all trajectories except those which are dead due to cold waters
        (according to law SCL = f(age) and min_temp = f(SCL))

        ms: markersize
        """

        if traject=='all':
            traject=range(self.ntraj)

        ########## ZONE OCEAN INDIEN + MER CHINOISE/INDONESIENNE + BAY TOMINI
        traj = []
        for zone in [1,4,5,6]:
            traj_zone = self.list_traj_in_zone_at(zone=zone,year=year,traj=traject)
            print('Zone ', zone, 'Turtles', len(traj))
            traj = np.hstack((traj,traj_zone))
        traj = np.array(traj,dtype=int)
        traj = np.unique(traj)
        print('Zone concatenee ; Turtles', len(traj))
        print(traj)

        #Positions and temperature loading
        x = self.x.table[:, traj]
        y = self.y.table[:, traj]
        t = self.temp.table[:, traj]

        #Variables : converting age into SCL
        Linf = 1.43 
        k = 0.226
        age0 = -0.17
        L0 = 0.053

        #Temperature en fonction de SCL - necessary variables
        a = (Tmin_hatch - Tmin_adult)/(L0-Linf)
        b = Tmin_hatch - a*L0


        #Local variables
        temp = 0.0*np.zeros(t.shape[0]) # critical temperatures depending on age
        cold_death = len(temp)*np.ones(len(traj)) # dead particles
        counter = 0

        #Critical temperatures for different ages
        #temp[10] => age = 9 days
        print('  - Calcul des temperatures tolerees')
        for time in range(len(temp)): #nb de lignes => nb de jours
            age = float(time) / 365.0
            #print('age:', age)
            SCL = Linf*(1-exp(-k*(age-age0)))
            #print('SCL:', SCL)
            temp[time] = a*SCL + b

        print('  - Parcours des particules pour voir leur mortalite')
        for part in range(len(traj)):
            #print(part)
            for time_ind in range(len(temp)):
                if (temp[time_ind] > t[time_ind,part]):
                    #print(ntraj, 'temp[time_ind]', temp[time_ind], 't[time_ind,part]', t[time_ind,part])
                    cold_death[part] = time_ind
                    counter = counter + 1
                    #print(counter)
                    break
        # OK !
        
        #set_printoptions(threshold=nan)
        if plot==True: 
            print('   - Plot' )
            ### Plots the time of death of each particle
            #plt.plot(traj, cold_death, '-', lw=1)
            #xlabel('Particle number')
            #ylabel('Time of death (day)')
            #title('T = f(SCL)')
            #grid(True)

            ### Courbes du nombre de tortues en vie en fonction de l'age
            compt = np.zeros(len(temp))
            for time_ind in range(len(temp)):
                for turt in range(len(traj)):
                    if (cold_death[turt] > time_ind):
                        compt[time_ind] = compt[time_ind]+1   
 
            ### Pourcentage de survie
            #plt.plot(range(len(temp)), compt/len(traj), '-', lw=2)
            #ylabel('Turtles left (ratio)')          
            
            ### Pourcentage de survie
            plt.plot(range(len(temp)), compt/len(traj), '-', lw=2)
                     
            ### Nombre de survie
            #plt.plot(range(len(temp)), compt, '-', lw=2)
        
        for zone in [2,3]:
            traj = self.list_traj_in_zone_at(zone=zone,year=year,traj=traject)
            ntraj = len(traj)
            ind = np.zeros(ntraj)
            print('Zone ', zone, 'Turtles', len(traj))

            #Positions and temperature loading
            x = self.x.table[:, traj]
            y = self.y.table[:, traj]
            t = self.temp.table[:, traj]

            #Variables : converting age into SCL
            Linf = 1.43 
            k = 0.226
            age0 = -0.17
            L0 = 0.053

            #Temperature en fonction de SCL - necessary variables
            a = (Tmin_hatch - Tmin_adult)/(L0-Linf)
            b = Tmin_hatch - a*L0


            #Local variables
            temp = 0.0*np.zeros(t.shape[0]) # critical temperatures depending on age
            cold_death = len(temp)*np.ones(len(traj)) # dead particles
            counter = 0

            #Critical temperatures for different ages
            #temp[10] => age = 9 days
            print('  - Calcul des temperatures tolerees')
            for time in range(len(temp)): #nb de lignes => nb de jours
                age = float(time) / 365.0
                #print('age:', age)
                SCL = Linf*(1-exp(-k*(age-age0)))
                #print('SCL:', SCL)
                temp[time] = a*SCL + b

            print('  - Parcours des particules pour voir leur mortalite')
            for part in range(len(traj)):
                #print(part)
                for time_ind in range(len(temp)):
                    if (temp[time_ind] > t[time_ind,part]):
                        #print(ntraj, 'temp[time_ind]', temp[time_ind], 't[time_ind,part]', t[time_ind,part])
                        cold_death[part] = time_ind
                        counter = counter + 1
                        #print(counter)
                        break
            # OK !
            
            #set_printoptions(threshold=nan)
            if plot==True: 
                print('   - Plot' )
                ### Plots the time of death of each particle
                #plt.plot(traj, cold_death, '-', lw=1)
                #xlabel('Particle number')
                #ylabel('Time of death (day)')
                #title('T = f(SCL)')
                #grid(True)

                ### Courbes du nombre de tortues en vie en fonction de l'age
                compt = np.zeros(len(temp))
                for time_ind in range(len(temp)):
                    for turt in range(len(traj)):
                        if (cold_death[turt] > time_ind):
                            compt[time_ind] = compt[time_ind]+1   
     
                ### Pourcentage de survie
                plt.plot(range(len(temp)), compt/len(traj), '-', lw=2)
                         
                ### Nombre de survie
                #plt.plot(range(len(temp)), compt, '-', lw=2)
            

        xlabel('Time (days)')
        #ylabel('Turtles left')
        ylabel('Turtles left (ratio)') 
        title(titre)
        #legend('123456')
        legend('123')
        grid(True)



    def autolabel(self,rects):
        """
        Display routine, puts height of an histogram bar

        called by: hist_zones
        should not be called directly
        """
        # attach some text labels
        for rect in rects:
            height = rect.get_height()
            plt.text(rect.get_x()+rect.get_width()/2., 5+height, '%d'%int(height),
                    ha='center', va='bottom',fontweight='bold',size='small')


    def mean_var_by_zone(self,ax,nom_var, year=1):
        """
        Plots primary production mean
        datas are separeted with the geographic zone of the particles at age=year

        ax: figure axes where to plot data
        year: float, age of the particles for the test, default is 1
        """

        if (nom_var.lower() == 'sst') or (nom_var.lower() == 'temp'):
            var = self.temp.table
            nb_filtre=2*year
        elif (nom_var.lower() == 'pp'):
            var = self.pp.table
            nb_filtre=10*year
        else:
            print('variable '+str(nom_var)+' non prise en charge\n,' +\)
            exit()

        # important indexes
        imin=0
        imax=year*365
        #nb_filtre=10*year

        # table initializations
        tt=np.array(np.arange(imin,imax),dtype=np.float32)
        mean_var=np.zeros(tt.shape,dtype=np.float32)
        mean_var_f=np.zeros(tt.shape,dtype=np.float32)
        mask=np.ones(tt.shape,dtype=np.bool)

        for zone in range(1,7):
            # masked position array by zone
            jlist=self.list_traj_in_zone_at(zone,year)
            sub_var=var[imin:imax,jlist]
            # mean tables
            mean_var=np.ma.mean(sub_var,1)
            for j in range(imin+nb_filtre,imax-nb_filtre):
                mean_var_f[j]=np.mean(mean_var[j-nb_filtre:j+nb_filtre+1])
                mask[j]=0
            mean_var_f=np.ma.MaskedArray(mean_var_f,mask=mask)
            # plots
            ax.plot(tt,mean_var_f,color[zone-1],lw=2)

        ax.legend(('(NC) North China sea','(NP) North Pacific','(SP) South Pacific','(IO) Indien Ocean', '(SC) South China sea', '(SS) Sulu sea'),loc='best')


    def list_traj_through_poly(self, lats, lons, closed=False):
        """
        Calculates for each time step which trajectories go through any faces
        of the polygon determined by lats/lons arrays. Distinction is made 
        between trajectories which cross a segment upwards and downwards.

        The function returns 3 1D arrays:
        ltime contains time indexes for each crossing
        lind contains trajectory indexes
        lsens contains 1, 2, 3,... when the trajectory crosses the 1st, 2nd, 3rd
        etc. face downwards ; -1, -2, -3, etc. when upwards
        
        These 3 arrays can be used for accurate statistics about a specific
        region, for example data_map.poly function prints statistics for
        each segment and plots all the trajectories in lind.
        """

        lats = np.array(lats).tolist()
        lons = np.array(lons).tolist()


        ltime = []
        lind = []
        lsens = []

        if len(lats) != len(lons):
            print('Les vecteurs de longitude/latitude doivent etre de meme taille')
            exit()

        if closed==True and not (lats[0] == lats[-1] and lons[0] == lons[-1]):
            lats.append(lats[0])
            lons.append(lons[0])

        for seg in range(len(lats)-1): # iteration sur les segments du polygone
            latA = lats[seg]
            latC = lats[seg+1]
            lonA = lons[seg]
            lonC = lons[seg+1]


            if(latA == latC and lonA == lonC):
                print('Deux points consecutifs doivent etre distincts !')
                exit()

            # Equation de la droite :
            #   (yC - yA)*x + (xA - xC)*y + xC*yA - xA*yC = 0
            #   (lonC - lonA)*x + (latA - latC)*y + latC*lonA - latA*lonC = 0
            #   A*x + B*y + C = 0
            # x <=> lat, y <=> lon


            # precalcul
            diffLat = latC-latA
            diffLon = lonC-lonA
            C = latC*lonA - latA*lonC

            # iteration sur les indices temporels
            for i in range(self.nb_output-1):
                # BD = segment de trajectoire entre les temps i et i+1
                latB = self.lat.table[i,:] # self.lat.table.shape = (nb_output, ntraj)
                latD = self.lat.table[i+1,:]
                lonB = self.lon.table[i,:]
                lonD = self.lon.table[i+1,:]

                det1=diffLat*(lonB-lonA)-diffLon*(latB-latA)
                det2=diffLat*(lonD-lonA)-diffLon*(latD-latA)
                det3=(latD-latB)*(lonA-lonB)-(lonD-lonB)*(latA-latB)
                det4=(latD-latB)*(lonC-lonB)-(lonD-lonB)*(latC-latB)
                # Les points A, B, C et D forment un quadrilatere et la trajectoire
                # (le segment) BD coupe le segment AC si ABCD est convexe.
                # ABCD est convexe ssi det1*det2 < 0 et det3*det4 < 0

                # calcul de la position de B par rapport a la droite (AC)
                pos = diffLon*latB - diffLat*lonB + C # > 0 ssi latB au dessus de AC

                # trajectoires dans le sens "positif" +1
                # Test sur la longitude pour regler les cas ou elle passe de -180 a
                # + 180 : dans ce cas, l'ecart est d'environ 360 degres.
                ind = np.ma.where((pos >= 0) & (det1*det2 < 0.0) & (det3*det4 < 0.0) \
                               & (abs(lonD - lonB) < 180.0))[0].tolist()
                ltime=ltime+len(ind)*[i]
                lind=lind+ind
                lsens = lsens+len(ind)*[(seg+1)]

                if(len(ind)>0):
                    print('Age : ' + str(i) + ' jour(s)')
                    print(' particule(s) (sens +1) :')
                    print(' '+str(ind))

                # trajectoires dans le sens "negatif" -1
                ind = np.ma.where((pos < 0) & (det1*det2 < 0.0) & (det3*det4 < 0.0) \
                               & (abs(lonD - lonB) < 180.0))[0].tolist()
                ltime=ltime+len(ind)*[i]
                lind=lind+ind
                lsens = lsens+len(ind)*[-(seg+1)]

                if(len(ind)>0):
                    print('Age : ' + str(i) + ' jour(s)')
                    print(' particule(s) (sens -1) :')
                    print(' '+str(ind))

            # difference entre le nombre de trajectoires qui traversent dans un
            # sens ou dans l'autre. Cela donne une idee de la tendance qu'on les
            # trajectoires a faire "demi-tour" une fois la ligne traversee pour 
            # revenir du cote initial, ou plutot a ne pas revenir sur leur pas
            # une fois la ligne franchie
            print('Differentiel = '+str(np.sum(np.array(lsens))))

        return np.array(ltime),np.array(lind),np.array(lsens)


############################ CLASSE FILLE DATA_MAP #################################


class data_map(data):
    """
    class map adds to class data the positions on a defined map
    contains also all map display routines
    """
    def __init__(self,filename,m,xoffset=0,yoffset=0):
        """
        Called with the path to the netcdf Ariane output file as an argument
        load all data and compute the positions on map

        filename: string, path to netcdf file
        m: Basemap reference
        xoffset: float, offfset for parallels name, default is 0
        yoffset: float, offfset for meridians name, default is 0
        """
        #load classical data
        data.__init__(self,filename)
        self.map=m

        # allow cylindrical maps with higher longitude than 180, e.g. lon1=140, lon2=240 for pacific view
        lon_c=np.array(self.lon.table)

        if (self.map.projection=='cyl'):
            lon_c[np.where(lon_c<self.map.llcrnrlon)]+=360


        # calculate and save positions on map
        x,y=m(lon_c,self.lat.table)
        self.x0,self.y0=m(self.lon0,self.lat0)
        self.xoffset=xoffset
        self.yoffset=yoffset
        self.x=var(x,self.mask,np.float32,'x - position on map')
        self.y=var(y,self.mask,np.float32,'y - position on map')
        print('   => positions on map loaded')

        # Masque les donnees situees en dehors de la zone de dessin
        if(self.map.projection == 'cyl'):
            ind = np.where((lon_c > self.map.urcrnrlon) |
                           (lon_c < self.map.llcrnrlon))
            self.x.table.mask[ind] = True
            self.y.table.mask[ind] = True
            self.lon.table.mask[ind] = True
            self.lat.table.mask[ind] = True

    def fill_map(self,lat_space=10, lon_space=20):
        """
        Fills the map on the current figure
        Draws continents, countries, coastlines, parallels and merdidians
        Also shows the nesting beach
        lat_space, lon_space : spacings between latitudes and longitudes

        NB : function create_map already calls fill_map
        """

        self.map.fillcontinents(color='0.4',zorder=-1)

        # draw parallels and meridians.
        if (self.map.projection=='cyl'):        
            # Calcul auto des etiquettes lon/lat, pour qu'elles ne sortent pas de la figure
            def getTicks(lmin, lmax, step):
                lmaxabs = (int(max(abs(lmax), abs(lmin)))/step+1)*step
                return np.intersect1d(np.arange(lmin+1, lmax), np.arange(-lmaxabs, lmaxabs, step))

            # [0,1,0,0] ; [0,0,0,1]
            self.map.drawparallels(getTicks(self.map.llcrnrlat, self.map.urcrnrlat, lat_space), labels=[0,1,0,0],xoffset=self.xoffset, fontsize=12)
            self.map.drawmeridians(getTicks(self.map.llcrnrlon, self.map.urcrnrlon, lon_space), labels=[0,0,0,1],yoffset=self.yoffset, fontsize=12)
        else: # si la projection n'est pas cylindrique, les variables llcrnr* et
              # urcrnr* ne sont pas forcement definies, donc dans le doute on 
              # affiche tous les meridiens et paralleles
            self.map.drawparallels(np.arange(-80.,81.,10.),labels=[0,1,0,0],xoffset=self.xoffset, fontsize=12)
            self.map.drawmeridians(np.arange(-180.,360.,20.),labels=[0,0,0,1],yoffset=self.yoffset, fontsize=12)

        # draw coast&countries
        self.map.drawcoastlines(color='0.6',linewidth=0.5,zorder=2)
        self.map.drawcountries(color='0.6',linewidth=0.3,zorder=2)
        self.map.plot([self.x0,self.x0],[self.y0,self.y0],markerfacecolor='w',markeredgecolor='k',marker='o',ms=6,mew=1.,zorder=999)


    def create_map(self, titre = "", lat_space=10, lon_space=20):
        """
        Creates a new figure and calls fill_map
        """
        fig=plt.figure(figsize=(11,7), dpi=100, facecolor='w', edgecolor='w')
        fig.subplots_adjust(bottom=0.02,top=0.98,left=0.01,right=0.99,hspace=0.2,wspace=0.08)
        #fig=plt.figure(figsize=(13.5,9), dpi=100, facecolor='w', edgecolor='w')        
        #fig.subplots_adjust(bottom=0.01,top=0.96,left=0.01,right=0.99,hspace=0.1,wspace=0.3)
        if titre!="":
            plt.title(titre)
        self.fill_map(lat_space=lat_space,lon_space=lon_space)
        return fig


    def display_trajectories(self,\
                          lw=0.005, ms=0.0, col='b', alpha=0.5,\
                          traj=None, age=None,\
                          trajmin=None, trajmax=None,\
                          agemin=None, agemax=None,\
                          var=None, show_cb=True, pad=+0.03,\
                          vmin=None, vmax=None,\
                          cold=0, Tmin_adult=15.0, Tmin_hatch=24.0,\
                          tolerance=False,\
                          points_death=False,\
                          pp_routes=False, temp_routes=False, traj_mortes=False\
                          ):
        """
        Highly parametrized display function, some settings however are still
        hard-coded (colorbar).

        lw: linewidth
        ms: markersize
        col: color
        alpha: ?

        traj: trajectories to be displayed ; [trajmin trajmax]
        age: restriction on the age of the particles ; [agemin agemax]

        var: variable (<=> color), {'time','sst','pp'}
        show_cb: if True a colorbar is added under the map
        pad: fraction of original axes between colorbar and new image axes = 0.035
        vmin, vmax: used for pp only

        cold: 0 if no cold_death, 1 if cold_death non evol, 2 if cold_death evol
        agemax: in years, age wanted
        Tmin_adult: if cold=1 Tmin ; if cold=2 Tmin_adult
        tolerance: if cold, number of days tolerated

        points_death: if True, displays the points where particles die

        calls self.evol_temperature, self.noevol_temperature
        returns number of living turtles at age
        """
        
        
        # Selection des trajectoires
        traj, age, cold_death = self.cold_selection(traj=traj,age=age,agemin=agemin,agemax=agemax,cold=cold,Tmin_adult=Tmin_adult,Tmin_hatch=Tmin_hatch, tolerance=tolerance, traj_mortes=traj_mortes)

        lturt = len(traj)
        print('Number of survivors :', lturt)

        traj = np.array(traj).tolist() ##tolist() convert to a standard python list.
        
        if traj == []:
            print('\n traj==[] : Aucune trajectoire a afficher')
            return lturt

        age = np.array(age).tolist()
        
        if age == []:
            print('\n age==[] : Aucune trajectoire a afficher')
            exit()


        # On map, displays points of death (red crosses) and returns number of death particles
        
        if points_death==True:


            print('points_death=True')
            print('cold_death', cold_death)
            print('traj', traj)

            if cold_death == []:
                print('\n cold_death==[] : Aucune particule morte')
                exit()

            # Recup des part mortes avant agemax
            if traj_mortes==True:
                trajm = np.array(traj)
	    else :
                indm = np.where(cold_death < agemax*365-1)[0] ## Pour quoi le -1 si c'est strictement inférieur.
                ## indm = indm[0] 
                trajm = np.array(traj)
                trajm = trajm[indm]
            ##indm = indm[0]
            ##trajm = np.array(trajex)
           
            ##trajm = trajex[indm]
            
            print('trajm', trajm)
            print(len(trajm))
            ##print(len(indm))
            

            nb_traj = np.size(trajm)
            nb_output = np.size(age)
            final_point_x = np.zeros(nb_traj)
            final_point_y = np.zeros(nb_traj)
    
            print('FINAL POINT'                )
            for t in range(nb_traj):
                indicec=int(cold_death[t]) ##indicec correspond à l'age entier de chaque particule.
                indicouz=np.where(cold_death>=365*agemax-1)
                self.x.table[indicouz,t]=365*agemax

                if indicec!=nb_output+1:
                    final_point_x[t] = self.x.table[indicec,t] ##RECUPERATION DES LONGITUDE ET LATITUDE DU JOUR DE MORT.
                    final_point_y[t] = self.y.table[indicec,t]

            final_point_x = np.array(final_point_x)
            final_point_y = np.array(final_point_y)

            print('x', final_point_x)
            print('y', final_point_y)

            print('IND')
            indx = np.where(final_point_x != 0)
            indy = np.where(final_point_y != 0)
            print('indx', indx[0])
            print('indy', indy[0])

            indx = [int(i) for i in indx[0]] ## Indices des particules pour lesquelles les positions ne sont pas nulles
            indy = [int(i) for i in indy[0]]

            print('final_point_x')
            final_point_x = final_point_x[indx]
            final_point_y = final_point_y[indy]

            print('plot')
            p=self.map.plot(final_point_x, final_point_y, lw=0.0, markersize=5, marker='x', c='r')
            
            return nb_traj

        # On map, displays best routes for pp (criteria to be chosen)
        if pp_routes==True:
            print('pp=True')
            print('cold_death', cold_death)
            print('traj', traj # survivantes, ok )

            nb_traj = np.size(traj)
            nb_output = np.size(age)

            compt_pp = np.zeros(nb_traj)

            print('COMPT PP et SEUIL')
            moy = 0                
            for p in range(nb_traj):
                for t in range(nb_output):
                    compt_pp[p] = compt_pp[p] + self.pp.table[t,traj[p]]
                moy = moy + compt_pp[p]
            print('compt_pp', compt_pp)
            
            print('max', max(compt_pp))

            # Criteria
            seuil = max(compt_pp) - 15000
            #seuil = moy/nb_traj
            #seuil = 11239.17 #moyenne des 4 trajectoires
            print('seuil', seuil)
        
             
            indinf = np.where(compt_pp <= seuil)
            indsup = np.where(compt_pp > seuil)

            indinf = [int(i) for i in indinf[0]]
            indsup = [int(i) for i in indsup[0]]

            trajinf = traj[indinf]
            trajsup = traj[indsup]

            print('trajinf', len(traj[indinf]), traj[indinf])
            print('trajsup', len(traj[indsup]), traj[indsup])

            print('self.x.table.shape', self.x.table.shape)
            print('self.y.table.shape', self.y.table.shape)

            print('age', age)

            xinf = self.x.table[:,trajinf][age,:]
            yinf = self.y.table[:,trajinf][age,:]
            xsup = self.x.table[:,trajsup][age,:]
            ysup = self.y.table[:,trajsup][age,:]

            print('xinf',xinf.shape, xinf)
            print('yinf',yinf.shape, yinf)
            print('xsup',xsup.shape, xsup)
            print('ysup',ysup.shape, ysup)

            print('plot')
            p1=self.map.plot(xinf, yinf, ls='-', ms='.', c='b', lw=lw, markersize=ms)
            print('p1 ok')
            p2=self.map.plot(xsup, ysup, ls='-', ms='.', c='r', lw=0.8, markersize=ms)
            print('p2 ok')

            return nb_traj

        # On map, displays best routes for temperature (criteria to be chosen)
        if temp_routes==True:
            print('temp=True')
            print('cold_death', cold_death)
            print('traj', traj # survivantes, ok )

            nb_traj = np.size(traj)
            nb_output = np.size(age)

            compt_temp = np.zeros(nb_traj)

            print('COMPT temp et SEUIL')
            moy = 0                
            for p in range(nb_traj):
                for t in range(nb_output):
                    compt_temp[p] = compt_temp[p] + self.temp.table[t,traj[p]]
                moy = moy + compt_temp[p]
            print('compt_temp', compt_temp)
            compt_temp = compt_temp/360
            print('compt_temp_norm', compt_temp)
            
            print('max', max(compt_temp))

            # Criteria
            seuil = max(compt_temp) - 0.25
            #seuil = moy/nb_traj
            #seuil = 11239.17 #moyenne des 4 trajectoires
            print('seuil', seuil)
        
             
            indinf = np.where(compt_temp <= seuil)
            indsup = np.where(compt_temp > seuil)

            indinf = [int(i) for i in indinf[0]]
            indsup = [int(i) for i in indsup[0]]

            trajinf = traj[indinf]
            trajsup = traj[indsup]

            print('trajinf', len(traj[indinf]), traj[indinf])
            print('trajsup', len(traj[indsup]), traj[indsup])

            print('self.x.table.shape', self.x.table.shape)
            print('self.y.table.shape', self.y.table.shape)

            print('age', age)

            xinf = self.x.table[:,trajinf][age,:]
            yinf = self.y.table[:,trajinf][age,:]
            xsup = self.x.table[:,trajsup][age,:]
            ysup = self.y.table[:,trajsup][age,:]

            print('xinf',xinf.shape, xinf)
            print('yinf',yinf.shape, yinf)
            print('xsup',xsup.shape, xsup)
            print('ysup',ysup.shape, ysup)

            print('plot')
            p1=self.map.plot(xinf, yinf, ls='-', ms='.', c='b', lw=lw, markersize=ms)
            print('p1 ok')
            p2=self.map.plot(xsup, ysup, ls='-', ms='.', c='r', lw=0.8, markersize=ms)
            print('p2 ok')

            return nb_traj                           
        
        #########################
        # DISPLAYS TRAJECTORIES #
        #########################

        x = self.x.table[:, traj][age, :]
        y = self.y.table[:, traj][age,:]

        nb_traj = np.size(traj)
        nb_output = np.size(age)

        # Choix de la variable a afficher
        str_var = str(var).lower()
        cb_title=None
        if str_var == 'time':
            col = np.meshgrid(np.arange(self.ntraj), np.arange(self.nb_output))[1][:, traj][age,:]
            cb_title = 'Age (in years)'

        elif str_var == 'traj':
            col = np.meshgrid(np.arange(nb_traj), np.arange(nb_output))[0]

        elif (str_var == 'sst') or (str_var == 'temp'):
            col = self.temp.table[:, traj][age, :]
            cb_title = "Sea Surface Temperature ($^{\circ}$C)"

        elif str_var == 'pp':
            col = self.pp.table[:, traj][age, :]
            ind = np.where(np.array(col) < 1e10)
            x = x[ind]
            y = y[ind]
            
            # unite : mgC
            col = 12.*col[ind]
            if vmax is not None:
                col[np.where(np.array(col) > vmax)] = vmax
            cb_title = "Net Primary Production (mg C.m$^{-2}$.day$^{-1}$)"

        elif str_var == 'zone':
            col = self.geotracer.table[:, traj][age, :]
            cb_title = 'Zones'

        else:
            var = None
        
        # Affichage
        if var is None:
            if lw == 0.0:
                print('x', x)
                print('y', y)
                p=self.map.scatter(x, y, s=ms, c=col, edgecolor='none', alpha=alpha, vmin=vmin, vmax=vmax)
            else:
                p=self.map.plot(x, y, ls='-', ms='.', c=col, lw=lw, markersize=ms)
        else:
           
            ##
            print('NB_TRAJ= ',nb_traj)
            print('NB_JOURS= ',nb_output)
            flag=np.zeros((2191,5000))
	    print('Dimensions du flag : ', np.shape(flag))
            ##PARTIE MODIFIEE PAR ROMAIN
            ##On donne une valeur au drapeau pour chaque jours et pour toutes les particules==> 1 lorsque la particule est morte
            if traj_mortes==True:
                for indice in range(nb_traj):  
	            if cold_death[indice] < agemax*365+1:
                        print('COLD_DEATH[INDICE] NUMERO  ',indice,' = ', cold_death[indice])
                        flag[int(cold_death[indice]):,indice]=1      
                        print(flag[:,indice])
                        print('SHAPE = ',np.shape(flag))
            indd = np.where((self.map.llcrnrlon < np.array(x)) & (np.array(x) < self.map.urcrnrlon) & (self.map.llcrnrlat < np.array(y)) & (np.array(y) < self.map.urcrnrlat) & (np.array(flag)==1))

            print('LEN(INDD[0]) = ',len(indd[0]), 'LEN(X) = ', len(x)           )
            print('SHAPE INDD[0] = ',np.shape(indd[0]))
            print('FORME DE X', np.shape(x))
            print('FORMAT DU FLAG = ',np.shape(flag))
            print('FORME DE INDD = ',np.shape(indd))

            p=self.map.scatter(x[indd], y[indd], c=col[indd], s=ms, edgecolor='none', alpha=alpha)
            
            if (cb_title is not None) and (show_cb==True):
                if var=='pp':
                    if vmax is not None:
                        cb=plt.colorbar(orientation='horizontal', shrink=0.8, pad=pad, ticks=range(0, vmax+1, 100))
                        cb.ax.set_xticklabels(np.arange(0, vmax-1, 100).astype('|S3').tolist()+['>'+str(vmax)])
                    else:
                        cb=plt.colorbar(orientation='horizontal', shrink=0.8, pad=pad)
                elif var=='time':
                    cb=plt.colorbar(orientation='horizontal', shrink=0.8, pad=pad, ticks=range(0, 2191, 365))
                    cb.ax.set_xticklabels(range(0,7,1))
                else:
                    cb=plt.colorbar(orientation='horizontal', shrink=0.9, pad=pad)
                cb.set_label(cb_title)
            
        return lturt

    def zones(self,year=1,lw=0.05, ms=0.0, showcrash=True, traj='all', time='all', nb_zones=4):
        """
        Plots the trajectories with different colors for the zone on the current figure
        test on the zone is done at age = year
        trajectories are shown from age=0 to age = year

        year: float, age of the particles for the zone test, default is 1
        lw: float, allows to fix the line width, usually low for a better view
        traj: trajectories to display
        time: portion of trajectories to be displayed

        """
        ind=year*365
        if time == 'all':
            time = range(self.nb_output)

        if nb_zones==4:

            for zone in [2, 3, 4, 1]: # bleu au dessus
                jlist=self.list_traj_in_zone_at(zone,year)
                print(len(jlist))
                if len(jlist) == 0:
                    continue

                if traj != 'all':
                    #jlist = np.round(np.intersect1d(jlist, traj)).astype(int)
                    jlist = np.intersect1d(jlist, traj).astype(int)

                if showcrash==False:
                    tempx = self.x.table[-5:,:][:,jlist]
                    tempy = self.y.table[-5:,:][:,jlist]
                    # indices des particules qui ne se crashent pas
                    sublist = np.where((tempx.min(axis=0) != tempx.max(axis=0)) |
                                     (tempy.min(axis=0) != tempy.max(axis=0)))[0]
                    jlist = jlist[sublist]

                if lw != 0.0:
                    self.map.plot(self.x.table[:,jlist][time, :],self.y.table[:,jlist][time, :],
                                  ls='-', linewidth=lw, ms='.', markersize=ms,
                                  c=color[zone-1], zorder=0)
                else:
                    self.map.scatter(self.x.table[:,jlist][time, :],self.y.table[:,jlist][time, :],
                                     s = ms, edgecolor='none', c = color[zone-1])

                print('   => zone '+str(zone)+' completed')
                print(str(len(jlist)) + ' trajectoires tracees')

        elif nb_zones==6:
            #for zone in [2, 3, 4, 5, 6, 1]: # bleu au dessus
            jlist=[]

            for zone in [4,5,6,1]:
                jlist=self.list_traj_in_zone_at(zone,year)
                print(len(jlist))
                if len(jlist) == 0:
                    continue

                if traj != 'all':
                    #jlist = np.round(np.intersect1d(jlist, traj)).astype(int)
                    jlist = np.intersect1d(jlist, traj).astype(int)

                if showcrash==False:
                    tempx = self.x.table[-5:,:][:,jlist]
                    tempy = self.y.table[-5:,:][:,jlist]
                    # indices des particules qui ne se crashent pas
                    sublist = np.where((tempx.min(axis=0) != tempx.max(axis=0)) |
                                     (tempy.min(axis=0) != tempy.max(axis=0)))[0]
                    jlist = jlist[sublist]

                if lw != 0.0:
                    self.map.plot(self.x.table[:,jlist][time, :],self.y.table[:,jlist][time, :],
                                  ls='-', linewidth=lw, ms='.', markersize=ms,
                                  c=color[zone-1], zorder=0)
                else:
                    self.map.scatter(self.x.table[:,jlist][time, :],self.y.table[:,jlist][time, :],
                                     s = ms, edgecolor='none', c = color[0])

                print('   => zone '+str(zone)+' completed')
                print(str(len(jlist)) + ' trajectoires tracees')

            for zone in [2, 3]: # bleu au dessus
                jlist=self.list_traj_in_zone_at(zone,year)
                print(len(jlist))
                if len(jlist) == 0:
                    continue

                if traj != 'all':
                    #jlist = np.round(np.intersect1d(jlist, traj)).astype(int)
                    jlist = np.intersect1d(jlist, traj).astype(int)

                if showcrash==False:
                    tempx = self.x.table[-5:,:][:,jlist]
                    tempy = self.y.table[-5:,:][:,jlist]
                    # indices des particules qui ne se crashent pas
                    sublist = np.where((tempx.min(axis=0) != tempx.max(axis=0)) |
                                     (tempy.min(axis=0) != tempy.max(axis=0)))[0]
                    jlist = jlist[sublist]

                if lw != 0.0:
                    self.map.plot(self.x.table[:,jlist][time, :],self.y.table[:,jlist][time, :],
                                  ls='-', linewidth=lw, ms='.', markersize=ms,
                                  c=color[zone-1], zorder=0)
                else:
                    self.map.scatter(self.x.table[:,jlist][time, :],self.y.table[:,jlist][time, :],
                                     s = ms, edgecolor='none', c = color[zone-1])

                print('   => zone '+str(zone)+' completed')
                print(str(len(jlist)) + ' trajectoires tracees')


    def zones_season(self,year=1,season=1):
        """
        Same as zones, but computes only the particles of one given season

        year: float, age of the particles for the test, default is 1
        season: integer between 1 and the number of seasons, default is 1
        """

        imin=0
        imax=year*365
        jmin=(season-1)*self.ntraj/self.nb_season
        jmax=(season)*self.ntraj/self.nb_season


        for zone in range(1,7):
            list=self.list_traj_in_zone_at(zone,year)
            jlist=list[np.where((list>=jmin)&(list<jmax))]
            self.map.plot(self.x.table[imin:imax,jlist],self.y.table[imin:imax,jlist],color[zone-1]+'-',linewidth=0.2,zorder=0)
            print('   => zone '+str(zone)+' completed')


    def poly(self, lats, lons, closed=False, lw=0.005, debug=True):
        """
        Plots only the trajectories that enter a given polygon
        Also prints statistics about how many trajectories cross the segment for
        each time step.

        closed: if True, first and last points are linked so lats and lons 
        define a closed path.

        lw: linewidth

        debug: if set to 'Full', detailed (daily) statistics are printed. 
        if True, only a summary for each segment is shown.
        if False (or anything that is not True or 'Full'), no information is printed.

        calls: list_traj_through_poly
        """
        ltime,lind,lsens = self.list_traj_through_poly(lats, lons, closed)

        for i in range(1, np.max(lsens)+1):
            # indice des trajectoires qui passent dans un sens
            ind = lind[np.where(lsens==i)]
            ind_traj=np.unique(ind)
            # et dans l'autre
            nind = lind[np.where(lsens==-i)]
            nind_traj=np.unique(lind[np.where(lsens==-i)])

            # print(statistics)
            if debug==True or debug == 'full':
                print('')
                print(' segment ' + chr(65+i-1))
                print('============================================')
                print(str(len(ind))+' traversee(s) de haut en bas ('+str(len(np.unique(ind)))+' trajectoires distinctes)')
                ltime_downwards = ltime[np.where(lsens==i)] # indices des particules qui traversent de haut en bas
                if debug == 'full':
                    for age in np.unique(ltime_downwards):
                        print(' age : '+str(age)+' jour(s), '+str(np.size(np.where(ltime_downwards==age)))+' particule(s)')
                print('============================================')
                print(str(len(nind))+' traversee(s) de bas en haut ('+str(len(np.unique(nind)))+' trajectoires distinctes)')
                ltime_upwards = ltime[np.where(lsens==-i)] # indices des particules qui traversent de bas en haut
                if debug == 'full':
                    for age in np.unique(ltime_upwards):
                        print(' age : '+str(age)+' jour(s), '+str(np.size(np.where(ltime_upwards==age)))+' particule(s)')
                print('============================================')
                print('  Soit '+str(len(ind)+len(nind))+' traversee(s) de particules non distinctes au total')
                print('       '+str(len(np.unique(ind.tolist()+nind.tolist())))+' trajectoires distinctes')
                print('')

            # display
            if len(ind_traj) > 0:
                self.map.plot(self.x.table[:,ind_traj],self.y.table[:,ind_traj],ls='-',lw=lw,c='b',alpha=1,zorder=0)
            if len(nind_traj) > 0:
                self.map.plot(self.x.table[:,nind_traj],self.y.table[:,nind_traj],ls='--',lw=lw,c='b',alpha=1,zorder=0)   

        self.plotPolygon(lats, lons, closed)


    def propagation_month(self, month=1, lw=0.05, traj='all'):
        """
        Displays propagation of particles for a given month, each season is given
        a different color.
        month: has to be chosen between 1 (january) and 12 (december)
        
        ind: indexes array of trajectories to focus on. Leave it to 'all' to consider all trajectories

        Note: only tested with 1-year long seasons (ie trajectories start 
        regularly throughout the year)
        """

        print('')

        if traj=='all':
            traj=range(self.ntraj)

        print(str(len(traj)) + ' trajectoires')
        init_t = self.init_t[...,traj]
        xpos = self.x.table[...,traj]
        ypos = self.y.table[...,traj]

        for s in range(self.nb_season):
            if month=='all':
                ind = np.where((init_t>=(365*s))&(init_t<(365*(s+1))))[0]
            else:
                ind = np.where((init_t>=(365*s+30*(month-1)))&(init_t<(365*s+30*month)))[0]
            if ind!=[]:
                print('  season '+str(s+1)+', trajectories '+str(ind[0])+' to '+str(ind[-1]))
                self.map.plot(xpos[:,ind],ypos[:,ind],ls='-',lw=lw,c=color[s],alpha=1,zorder=0)
            self.map.plot(-1000,-1000,ls='-', c=color[s], label='year '+str(2002+s) + ' (' +str(len(ind)) + ')') # dummy (for legend)

        plt.legend(loc=0)

    def propagation_temperature(self, ms=0.05, temp=0, traj='all', mintime = 1, plot=False):
        """
        Plots all trajectories except those which spend more than mintime days
        in waters where temperature is below temp.

        ms: markersize
        temp: temperature threshold
        """
        if traj=='all':
            traj=range(self.ntraj)

        ntraj = len(traj)
        ind = np.zeros(ntraj)

        x = self.x.table[:, traj]
        y = self.y.table[:, traj]
        t = self.temp.table[:, traj]

        for tnum in range(ntraj):
            # indices des positions >= temp
            ind_zone = np.ma.where(t[:,tnum] >= temp)[0] 
            # ecart entre ces indices = 1+temps passe hors de la zone
            hors_zone = np.where(ind_zone[1:] - ind_zone[:-1] - 1 > mintime)[0] 
            if hors_zone.shape[0]==0: # si la trajectoire ne sort pas, retourne le temps max
                ind[tnum] = ind_zone[-1]+1
            else:
                # on ne tue les particules qu'apres avoir passe 15j < 19deg, donc on rajoute les mintime jours
                ind[tnum] = min(np.min(ind_zone[hors_zone])+1+mintime, ind_zone[-1]+1)

            if plot==True:
                self.map.scatter(x[:ind[tnum], tnum], y[:ind[tnum], tnum], c=t[:ind[tnum], tnum],
                                 s=ms, edgecolor='none', alpha=0.5)

        if plot==True:
            cb=plt.colorbar(orientation='horizontal', shrink=0.8)
            cb.set_label('SST')

        return ind

    def propagation_month_poly(self, month, lats, lons, closed=False, lw=0.05):
        """
        Same behaviour as propagation_month but only plots trajectories which 
        enter a given polygon

        month: between 1 and 12
        lats: list of latitude coordinates
        lons: list of longitude coordinates

        calls: list_traj_through_poly, propagation_month
        """
        # computes indexes of crossing trajectories
        ltime,lind,lsens = self.list_traj_through_poly(lats, lons, closed)
        
        ind = np.unique(lind) # indices des trajectoires qui franchissent la ligne

        # displays line segment
        self.plotPolygon(lats, lons, closed)

        # calls progatation_month with crossing trajectories indexes only
        self.propagation_month(month, lw, ind)

        return ind

    def plotPolygon(self, lats, lons, closed=False, lw=2):
        """
        Displays a polygon determined by two arrays of lat/lon coordinates
        """

        x, y = self.map(lons, lats)

        if closed==True:
            x.append(x[0])
            y.append(y[0])

        self.map.plot(x, y, lw=lw, ls='-', c='g', alpha=1)
        
        for i in range(len(x)-1):
            pass
            #plt.text((x[i]+x[i+1])/2, (y[i]+y[i+1])/2, chr(65+i), fontsize=20)


    def ind_distance_to_beach(self, 
                              lw=0.1, ms=1,                               
                              traj=None, trajmax=None, 
                              age=None, agemax=None,
                              seuil=None, distmax=None,
                              pos=None,
                              zones=True, showcrash=False, smooth=None):
        """
        Plots distance to beach (or a given position) vs time for each particle
        
        lw, ms: linewidth, markersize
        traj: indexes of trajectories to work with
            trajmax: same as traj=range(trajmax)
        age: restriction on the age of the particles
            agemax: same as age=range(agemax)
        seuil: if specified, only trajectories which get closer than the 
               specified distance to their spawning point after a 100-day 
               drift will be displayed
        distmax: trajectories going further than this distance from their 
                 spawning point won't be displayed
        pos: if not None, distances are computed between trajectories and 
             the given [lat, lon] position (in degrees)
        zones: if =True, the curves are colored depending on the current 
               geographical area of the particles
        showcrash: if =False trajectories that get stuck on the shores 
                   (and remain at a same position forever) are discarded
        smooth: smoothing factor, each point i in the distance sequence is 
                the mean of dist[i-smooth, i+smooth+1]
        """

        R = 6378.137 # rayon de la Terre en km (sphere GRS80)
        ind_traj = np.arange(self.ntraj)

        # Restriction des trajectoires a traiter
        if traj is None:
            if trajmax is not None:
                traj = range(trajmax)
            else:
                traj = range(self.ntraj)
        ind_traj = ind_traj[traj]

        if age is None:
            if agemax is not None:
                age=range(agemax)
            else:
                age=range(self.nb_output)
        
        lon = self.lon.table[age, :][:, traj] * np.pi/180. # Conversion degrad
        lat = self.lat.table[age, :][:, traj] * np.pi/180.
        tracer = self.geotracer.table[age, :][:, traj]
        ntraj = len(traj)
        noutput = len(age)
        nzones = int(tracer.max())

        if pos is not None:
            lat0 = pos[0] * np.pi/180
            lon0 = pos[1] * np.pi/180
        else:
            lat0 = lat[0, :]
            lon0 = lon[0, :]

        # Calcul de la distance http://fr.wikipedia.org/wiki/Distance_du_grand_cercle
        tab_dist = 2*R * np.arcsin(np.sqrt(np.sin((lat - lat0)/2)**2 +\
                                           np.cos(lat)*\
                                           np.cos(lat0)*\
                                           np.sin((lon - lon0)/2)**2))
        """# moins precis pour de petites valeurs
        tab_dist = R * np.arccos(np.sin(lat)*np.sin(lat[0,:]) +\
                                 np.cos(lat)*np.cos(lat[0,:])*np.cos(lon-lon[0,:]))
        """

        # Filtrage des particules echouees
        if showcrash == False:
            ind = np.ma.where(np.ma.min(tab_dist[-5:], axis=0) != np.ma.max(tab_dist[-5:], axis=0))[0]
            #print(str(len(ind_traj) - len(ind)) + ' particules echouees')
            tab_dist = tab_dist[:, ind]
            tracer = tracer[:, ind]
            ind_traj = ind_traj[ind]
            ntraj = len(ind)

        # Selection des trajectoires qui bouclent
        if seuil is not None:
            if pos is not None:
                timethresh=0
            else:
                timethresh=100
            ind = np.ma.where(np.ma.min(tab_dist[timethresh:, :], axis = 0) < seuil)[0]
            #ind = ind[np.where(tab_dist[50, ind] > 1000)[0]]#5traj mayotte
            tab_dist = tab_dist[:, ind]
            tracer = tracer[:, ind]
            ind_traj = ind_traj[ind]
            ntraj=len(ind)

        # Filtrage des particules qui s'eloignent trop
        if distmax is not None:
            ind = np.ma.where(np.ma.max(tab_dist, axis = 0) < distmax)[0]
            tab_dist = tab_dist[:, ind]
            tracer = tracer[:, ind]
            ind_traj = ind_traj[ind]
            ntraj=len(ind)          


        ## Affichage
        if ntraj == 0:
            print("\n Pas de trajectoire a afficher, sortie")
            exit()


        fig = plt.figure()

        # Lissage des donnees
        # pas terrible, les bords ne sont pas traites. et pas tres utile...
        if smooth is not None:
            aux = np.zeros((noutput-(2*smooth+1), ntraj))
            for i in np.arange(-smooth, smooth+1, 1):
                aux += tab_dist[smooth+i:-smooth+i-1,:]
            tab_dist[smooth:-smooth-1,:] = aux/(2*smooth+1)

        if zones == False:
            plt.plot(tab_dist, lw=1)
        else:
            time = range(noutput)
            cmap = mpl.colors.ListedColormap(color[:nzones])
            bounds = np.arange(0.5, nzones+1+0.5, 1.0)
            norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

            for tnum in range(ntraj):
                plt.scatter(time, tab_dist[:, tnum], c=tracer[:, tnum], \
                            s=ms, edgecolor='none', alpha=0.5, cmap=cmap, norm=norm)
                plt.plot(time, tab_dist[:, tnum], '-k', lw=lw)

            plt.axis([0, noutput, 0, tab_dist.max()])

            cb=plt.colorbar(orientation='horizontal', ticks=np.arange(1, nzones+1, 1), shrink=0.9)
            cb.set_label('zones traversees')
            plt.xlabel('temps en jours')
            plt.ylabel('distance (km)')

            print("\n"+str(ntraj) + " trajectoires tracees sur " + str(self.ntraj))

        return fig, ind_traj

    def plot_init_positions(self,lon_island=None, lat_island=None, plot=True):
        """
        Displays initial positions on map
        """
        vmin = None
        vmax = None
        ms = 0.02      
        alpha = 0.5
        lw = 0.02
        traj = np.arange(0, self.ntraj)
        x = self.x.table[:, traj][0, :]
        y = self.y.table[:, traj][0, :]
        print('init_x', x)
        print('init_y', y)
        if plot:
            p=self.map.plot(x, y, lw=0.0, markersize=2, marker='o', c='b')
            if (lon_island is not None) and (lat_island is not None):
                print('island plotted')
                p=self.map.plot(lon_island,lat_island, marker='p', markersize=15, color='#FFFFFF', markeredgecolor='k', markeredgewidth=4)
        return x,y


class var:
    """
    class var contains all information about one variable

    should not be directly called
    called by data and data_map
    """
    def __init__(self,tab,mask,type,name,vmin=0,vmax=0):
        # Seules les valeurs vraies sont prises en compte pour toutes opérations faites dans ce tableau.
        self.table=np.ma.masked_array(tab,mask=mask,fill_value=-1E10,dtype=type)
        self.name=name
        self.min=vmin
        self.max=vmax

    def show(self,data,numfig=1,sub='111'):
        fig=plt.figure(num=numfig, figsize=(8, 4), dpi=100, facecolor='w', edgecolor='w')
        fig.add_subplot(sub)
        plt.matshow(self.table.filled(),cmap=cm.jet,fignum=False,origin='upper',vmin=self.min,vmax=self.max)
        xx=np.arange(data.nb_season*12)*data.ntraj/data.nb_season/12
        yy=np.arange(data.nb_year)*365
        plt.xticks(xx,[])
        plt.yticks(yy,np.arange(data.nb_year),fontsize=16)
        plt.grid(lw=2,ls='--')
        CB=plt.colorbar()
        CB.ax.set_position([0.90,0.1,1,0.8])
        plt.title(self.name,fontsize=18)
        #fig.subplots_adjust(bottom=0.02,top=0.98,left=0.05,right=0.89,hspace=0.2,wspace=0.08)
        fig.subplots_adjust(bottom=0.01,top=0.96,left=0.01,right=0.99,hspace=0.2,wspace=0.08)
        plt.show()

