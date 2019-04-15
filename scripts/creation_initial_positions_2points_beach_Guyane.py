#!/usr/bin/python2.7
#-*- coding:utf-8 -*-




## creation du fichier initial_positions.txt

#import Scientific.IO.NetCDF as nc
from mpl_toolkits.basemap import Basemap
import numpy as np
import os, sys
import datetime
import random as rd
import matplotlib.pyplot as plt
import time


#personal libraries
import IOlib as IO
import librumeau as brum




#### declaration des classes
class nesting_period :
    """
    class nesting_period used as a structure
    contains the begin, end and peak of the nesting season
    """
    def __init__(self,first=0.,last=0.,peak=0.):
        self.begin = np.float32(first)
        self.end = np.float32(last)
        self.max = np.float32(peak)
        

class beach :
    """
    class beach is used as a structure
    contains all information needed about the nesting beach
    => launching area
    => nesting season
    => geographic coordinates of the beach
    => name of the beach
    => eventually the reference of those information
    """

    
    
    def __init__(self,
                time=nesting_period(),
                north_pt=[1,1],
                south_pt=[0,0],
                d_min = 26.125,
                d_max = 53.875,
                ref="no reference",
                beach_name="no name",
                file_name="no name"):
        self.time=time
        self.north_pt=north_pt
        self.south_pt=south_pt
        self.d_min=d_min
        self.d_min=d_min
        self.d_max=d_max
        

        self.ref=str(ref)
        self.beach_name=str(beach_name)
        self.correction_incubation()
    
    def correction_incubation(self):
        """
        add 60 days to the nesting season for the incubation time
        """
	print "\n=> correcting nesting periods with incubation times"
        self.time.begin+=0.0
        self.time.end+=0.0
        self.time.max+=0.0
        
class position_4D:
    """ 
    class position_4D contains one particle initial position
    x, y, z, t
    """
    def __init__(self,x=0.,y=0.,t=0.,z=-1.0):
        self.x=np.float32(x)
        self.y=np.float32(y)
        self.z=np.float32(z)
        self.t=np.float32(t)
    def write_line(self,f):
        """
        write as a proper line readable for Ariane the position in a given file
        """
        line=' '+"%.3f" % self.x+' '+\
             "%.3f" % self.y+'    '+\
             "%.1f" % self.z+'  '+\
             "%.2f" % self.t+'     '+\
             '1.0\n'

        line=line.replace(',','.')
        f.write(line)

class echantillon :
    """
    class echantillon contains an sample of initial positions and its computation
    contains also the reference to the beach structure used for computation
    """
    ## initialisation
    def __init__(self,beach):
        """
        initialization
        return the beach structure and an empty list of positions
        """
        self.positions=[]
        self.beach=beach

    ## public - to be used
    def gen_positions(self,nb_part,nb_year,lat_mat,lon_mat,coord_mode):
        """
        generate the sample
        the number of particles and number of seasons are given as arguments
        """
	print "\n=> generating initial positions"
        self.nb_part=nb_year*int(nb_part/nb_year)
        self.nb_year=nb_year     
        self.random_time()
        for j in range(self.nb_part):
            #print self.vect_t[j]
            self.random_position(self.vect_t[j],lat_mat,lon_mat,coord_mode)



    def gen_positions_all_year(self,nb_part,nb_year,lat_mat,lon_mat):
        """
        generate the sample
        the number of particles and number of seasons are given as arguments
        """
        self.nb_year=nb_year    
        
        h=61.*np.float32(nb_year)/np.float32(nb_part)

        self.vect_t=np.arange(0.5,61.*nb_year+0.5,h)
        self.nb_part=self.vect_t.shape[0]

        for j in range(self.nb_part):
            self.random_position(self.vect_t[j],lat_mat,lon_mat)


    def write_all(self,file):
        """
        write the sample in a given file
        """
        for j in range(E.nb_part):
            E.positions[j].write_line(file)


    def histogram(self):
        """
        display in a new figure the histogram of the sample
        as time index, one unit correspond to sampling time of the ORCA025 model used
        for exemple, for 6day_mean fields, one year is approximated to 61units=366j
        """
        data=self.vect_t
        plt.figure(num=1, figsize=(6,3), dpi=100, facecolor='w', edgecolor='w')
        n, bins, patches = plt.hist(data, 120, normed=1, facecolor='green', alpha=0.75)        
        plt.show()

    ## private 
    # time randomisation
    def random_time(self):
        """
        compute the set of intial times        
        create a sorted list of times self.vect_t
        
        called by : gen_position
        call to : echantillon_normal_cut
        should not be called directly
        """
        print "\n=> Computing the set of departure times"
        tt=[]

        counter=0
        for year in range(self.nb_year):

            #tt=tt+self.echantillon_normal_cut(year)
            tt=tt+self.echantillon_normal_asym(year)
        print("\n=> "+str(self.nb_part)+' initial positions computed')
        tt=np.array(tt,dtype=np.float32)
        tt.sort()
        self.vect_t=tt
        print 'fin departure time'
        #print tt


    def variance95(self):
        """
        compute the variance for the time-dimension normal distribution
        troncate approximatively 5% of the distribution (true if symetric values are used)
        
        called by : echantillon_normal_cut
        should not be called directly
        """
        M=self.beach.time.end
        m=self.beach.time.begin
        s=((M-m)/3.92)
        return s
    def variance98(self):
        """
        compute the variance for the time-dimension normal distribution
        troncate approximatively 2% of the distribution (true if symetric values are used)
        
        called by : echantillon_normal_cut
        should not be called directly
        """
        M=self.beach.time.end
        m=self.beach.time.begin
        s=((M-m)/(2*2.33))
        return s

    def echantillon_normal_asym(self,year):
        """
        compute the set of intial times for one year    
        random a sample following a troncated normal distribution
        the distribution is replicated nb_season times, pushed over 61units (366j) each year
        
        called by : random_time
        call to : variance95 or variance98
        should not be called directly
        """
        peak=self.beach.time.max
        m=self.beach.time.begin
        M=self.beach.time.end
        #print 'MMMMM',m,M,peak
        h=[]
        Ntot=self.nb_part/self.nb_year
        for day in np.arange(m+0.5,M+1.5):
            #print day
            if day<=peak:
                h.append(np.exp(-2*((day-peak)/(peak-m))**2))
            elif day>peak:
                h.append(np.exp(-2*((day-peak)/(M-peak))**2))
        h = np.array(h)
        S = sum(h)
        a = Ntot/S
        n_days = np.round(a*h)
        Ntot_bis= np.sum(np.round(a*h))
        err_round = Ntot-Ntot_bis
        #print err_round
        n_days[np.where(n_days==np.max(n_days))[0][np.shape(np.where(n_days==np.max(n_days))[0])[0]/2]] += err_round
        #ind_max = np.where(n_days==np.max(n_days))[0][np.shape(np.where(n_days==np.max(n_days))[0])[0]/2]
        #print ind_max
        #n_days[ind_max] += 2
        #n_days[ind_max-1] += 1
        #n_days[ind_max-2] += 1
        #n_days[ind_max-3] += 1
        #n_days[ind_max+1] += 1
        #n_days[ind_max+2] += 1
        #n_days[ind_max+3] += 1
        time=np.zeros([Ntot])
        part_hatched=0
        for i in range(np.shape(n_days)[0]):
           #print "CP0", time
            time[int(part_hatched):int(part_hatched+n_days[i])]=(np.float(m+i+0.5))+365.25*year
            part_hatched = part_hatched+n_days[i]

        #print list(time)
        plt.hist(time,bins=np.arange(m,M,1))
        plt.show()
        return list(time)


    def echantillon_normal_cut(self,year):
        """
        compute the set of intial times for one year    
        random a sample following a troncated normal distribution
        the distribution is replicated nb_season times, pushed over 61units (366j) each year
        
        called by : random_time
        call to : variance95 or variance98
        should not be called directly
        """
        peak=self.beach.time.max
        m=self.beach.time.begin
        M=self.beach.time.end
        #print 'MMMMM',m,M
        taille=self.nb_part/self.nb_year
        sigma=self.variance95() #can be change to variance98
        
        time=range(taille)
        for j in range(taille):
            ####################
            #Normal distribution
            ####################
            a=rd.normalvariate(peak,sigma)
            while a<m or a>M :
                a=rd.normalvariate(peak,sigma)
                #print a
            ####################
            #Uniform ditrib
            ####################
            #a = rd.uniform(m,M)

            time[j]=np.floor(a+365.25*year)+0.5
        #plt.hist(time,bins=np.arange(m,M,1))
        #plt.show()
        return time

    # space randomisation
  
    def random_position(self,tt,lat_mat,lon_mat,coord_mode):
        """
           beach and dropping zone in (X,Y) coordinates ((X,Y) is a rotation of the geographical coordinates (x,y) of angle theta)
    
            A ----------------C--------D
            |                 |dropping|
            |B                |  zone  |
            |E                |        | 
            |A                |        |
            |C                |        |
            |H                |        |
            |                 |        |
            B ----------------F--------E
             <--------------->
                    dmin
             <------------------------->
                         dmax
        """
        # Génère les coordonnées des points extrèmes de la zone de largage dans le repère (X,Y)
        # obenu par rotation du repère géographique d'angle theta (angle entre le segment
        # représentant la plage et le Nord).
        
        #lecture des coordonnées dans le repère géographique et conversion en radians :
        xA=self.beach.north_pt[1]*np.pi/180
        yA=self.beach.north_pt[0]*np.pi/180
        xB=self.beach.south_pt[1]*np.pi/180
        yB=self.beach.south_pt[0]*np.pi/180
        Ox=xB
        Oy=yB
        theta=-np.arctan((xB-xA)/(yB-yA))
    
        #conversion de la translation dans le repère géographique
        xT_min,yT_min=brum.rotation(self.beach.d_min,0,Ox,Oy,-theta)
        xT_max,yT_max=brum.rotation(self.beach.d_max,0,Ox,Oy,-theta)
        
        #conversion des distances en angles (approx. du plan tangent)
        lat_moy=(yA+yB)/2
        R_terre=6371
        xT_min=xT_min/R_terre/abs(np.cos(lat_moy))
        xT_max=xT_max/R_terre/abs(np.cos(lat_moy))
        yT_min=yT_min/R_terre
        yT_max=yT_max/R_terre
        XT_min=np.sqrt(xT_min**2+yT_min**2)
        XT_max=np.sqrt(xT_max**2+yT_max**2)
        YT_min=0
        YT_max=0
    
        #conversion vers le repère rectifié :
        XB,YB=brum.rotation(xB,yB,Ox,Oy,theta)
        XA,YA=brum.rotation(xA,yA,Ox,Oy,theta) 
        XC=XA+XT_min
        XD=XA+XT_max
        XF=XB+XT_min
        XE=XB+XT_max
        YC=YA+YT_min
        YD=YA+YT_max
        YF=YB+YT_min
        YE=YB+YT_max
        
        #Tirage aléatoire de coordonnées (X,Y) dans le repère rectifié :

        #Test pour ne pas tirer de positions sur la terre

        x_grid = 0
        y_grid = 0
        k = 0
        #print temperature
        
        while np.isclose(temperature[0,int(y_grid),int(x_grid)],1e34) or k ==0:
            
            X=rd.uniform(XC,XD)
            Y=rd.uniform(YF,YC)
    
    
            # Conversion vers le repère géographique:
            xC,yC=brum.rotation(XC,YC,Ox,Oy,-theta)
            xF,yF=brum.rotation(XF,YF,Ox,Oy,-theta)
            x,y=brum.rotation(X,Y,Ox,Oy,-theta)
        
            #conversion en degrés
            x=x*180/np.pi
            y=y*180/np.pi

            #conversion en coordonnées de grille
            #print x,y,lon_mat,lat_mat
            x_grid,y_grid=brum.geo_to_grid(x,y,lon_mat,lat_mat)
            k = k+1

        if coord_mode == 'grid':
            self.positions=self.positions+[position_4D(x_grid,y_grid,tt)]
        elif coord_mode == 'latlon':
            self.positions=self.positions+[position_4D(x,y,tt)]






def YALIMAPO_HATTES_SUMMER(nesting_year,d_min,d_max):

    """                                                                                                   
    Nesting from 15 march to 15 august 
    Hatching from 15 may to 15 october
    REF = Study of a bimodal nesting season for leatherback turtles (Dermochelys coriacea) in French Guiana
    Chevalier et al, 1999
    """

    start_day = 135+(365*(nesting_year-2002))
    end_day = 288+(365*(nesting_year-2002))
    peak_day = (start_day+end_day)/2.
    time=nesting_period(start_day,end_day,peak_day) #ajouter test pour vérifier first < peak < last
    ref="no ref"
    beach_name="Yalimapo_hattes"
    file_name="Yalimapo_hattes"
    #north_pt=(-5.748738, -53.909114)  #coordonnées  lat,lon
    north_pt=(5.75001, 306.05)
    south_pt=(5.75, 306.0825)

    return beach(time,
                 north_pt,south_pt,
                 d_min,d_max,
                 ref,
                 beach_name,
                 file_name)  

def CAYENNE_SUMMER(nesting_year,d_min,d_max):

    """                                                                                                   
    Nesting from 15 march to 15 august 
    Hatching from 15 may to 15 october
    REF = Study of a bimodal nesting season for leatherback turtles (Dermochelys coriacea) in French Guiana
    Chevalier et al, 1999
    """

    start_day = 135+(365*(nesting_year-2002))
    end_day = 288+(365*(nesting_year-2002))
    peak_day = (start_day+end_day)/2.
    time=nesting_period(start_day,end_day,peak_day) #ajouter test pour vérifier first < peak < last
    ref="no ref"
    beach_name="Cayennes"
    file_name="Cayennes"

    north_pt=(4.95, 307.7)
    south_pt=(4.90, 307.77)

    return beach(time,
                 north_pt,south_pt,
                 d_min,d_max,
                 ref,
                 beach_name,
                 file_name)  

def YALIMAPO_HATTES_SUMMER_25km(nesting_year,d_min,d_max):

    """                                                                                                   
    Nesting from 15 march to 15 august 
    Hatching from 15 may to 15 october
    REF = Study of a bimodal nesting season for leatherback turtles (Dermochelys coriacea) in French Guiana
    Chevalier et al, 1999
    """

    start_day = 135+(365*(nesting_year-2002))
    end_day = 288+(365*(nesting_year-2002))
    peak_day = (start_day+end_day)/2.
    time=nesting_period(start_day,end_day,peak_day) #ajouter test pour vérifier first < peak < last
    ref="no ref"
    beach_name="Yalimapo_hattes"
    file_name="Yalimapo_hattes"
    #north_pt=(-5.748738, -53.909114)  #coordonnées  lat,lon

    north_pt=(5.75003976421, 305.953266328)
    south_pt=(5.74997023579, 306.179233672)

    return beach(time,
                 north_pt,south_pt,
                 d_min,d_max,
                 ref,
                 beach_name,
                 file_name) 

def CAYENNE_SUMMER_25km(nesting_year,d_min,d_max):

    """                                                                                                   
    Nesting from 15 march to 15 august 
    Hatching from 15 may to 15 october
    REF = Study of a bimodal nesting season for leatherback turtles (Dermochelys coriacea) in French Guiana
    Chevalier et al, 1999
    """

    start_day = 135+(365*(nesting_year-2002))
    end_day = 288+(365*(nesting_year-2002))
    peak_day = (start_day+end_day)/2.
    time=nesting_period(start_day,end_day,peak_day) #ajouter test pour vérifier first < peak < last
    ref="no ref"
    beach_name="Cayennes"
    file_name="Cayennes"

    north_pt=(4.9905, 307.643353594)
    south_pt=(4.8595, 307.826646406)

    return beach(time,
                 north_pt,south_pt,
                 d_min,d_max,
                 ref,
                 beach_name,
                 file_name)  
                 
                 
def YALIMAPO_HATTES_SUMMER_25km_025sq(nesting_year,d_min,d_max):

    """                                                                                                   
    Nesting from 15 march to 15 august 
    Hatching from 15 may to 15 october
    REF = Study of a bimodal nesting season for leatherback turtles (Dermochelys coriacea) in French Guiana
    Chevalier et al, 1999
    """

    start_day = 135+(365*(nesting_year-2002))
    end_day = 288+(365*(nesting_year-2002))
    peak_day = (start_day+end_day)/2.
    time=nesting_period(start_day,end_day,peak_day) #ajouter test pour vérifier first < peak < last
    ref="no ref"
    beach_name="Yalimapo_hattes"
    file_name="Yalimapo_hattes"
    #north_pt=(-5.748738, -53.909114)  #coordonnées  lat,lon

    #north_pt=(5.841683, 305.960736)
    #south_pt=(5.706046, 306.256929)

    north_pt=(5.81, 305.99)
    south_pt=(5.72, 306.23)


    return beach(time,
                 north_pt,south_pt,
                 d_min,d_max,
                 ref,
                 beach_name,
                 file_name) 

def CAYENNE_SUMMER_25km_025sq(nesting_year,d_min,d_max):

    """                                                                                                   
    Nesting from 15 march to 15 august 
    Hatching from 15 may to 15 october
    REF = Study of a bimodal nesting season for leatherback turtles (Dermochelys coriacea) in French Guiana
    Chevalier et al, 1999
    """

    start_day = 135+(365*(nesting_year-2002))
    end_day = 288+(365*(nesting_year-2002))
    peak_day = (start_day+end_day)/2.
    time=nesting_period(start_day,end_day,peak_day) #ajouter test pour vérifier first < peak < last
    ref="no ref"
    beach_name="Cayennes"
    file_name="Cayennes"

    #north_pt=(5.050857, 360-52.386974)
    #south_pt=(4.838257, 360-52.174374)
    
    north_pt=(5.0325, 307.6075)
    south_pt=(4.8675, 307.7925)

    return beach(time,
                 north_pt,south_pt,
                 d_min,d_max,
                 ref,
                 beach_name,
                 file_name)  

                                                 
##execution_example=\
##""" 
#nombre de particules
nb_part_total = 2500
nb_part=2500

nb_year=1
nb_season=1
nesting_year = 2002

#La représentation de la zone de départ dans le papier luth NATL est légèrement différente de la zone de départ réellement utilisée:
#La longueur de la plage était un peu plus grande que 0.25°. Nous l'avons donc raccourcie pour que les zones de départ soient bien des carrée de 0.25°x0.25°
# centrés 40km au large des plages de ponte
d_min = 26.125
d_max = 53.875
coord_mode = 'grid'
#fichier a ecrire
dirout='/homelocal/otitaud/gitlab/lmtl-wp4/config/'
#file = open ( dirout+'Jamursba_summer_'+str(nesting_year)+'_2000turt'+'.txt', 'w' )
outfile = dirout+'/YalimapoCayenne_summer_25km_'+str(d_min)+'_'+str(d_max)+'_'+str(nesting_year)+'_'+str(nb_part_total)+'turt_version_gaussienne_glorys2v4'+'.v2.txt'
file = open (outfile, 'w' )


#lecture du fichier contenant les informations sur la grille pour la conversion (lon,lat) -> (x,y)
#infile='/data/FSHML/Tortues/LMTL-WP4/mercatorglorys2v4_1998_2015/sad/mesh_reg.glorys2v4.nc'
infile='/data/FSHML/Tortues/LMTL-WP4/mercatorglorys2v4_1998_2015/forcings/yyyymmdd/phys/mercatorglorys2v4_T_19980101.nc'
nc_dico=brum.read_nc(infile,['longitude','latitude','temperature'])
lon_mat = nc_dico['longitude']
lat_mat = nc_dico['latitude']
temperature = nc_dico['temperature']


E=echantillon(YALIMAPO_HATTES_SUMMER_25km_025sq(nesting_year,d_min,d_max))
E.gen_positions(nb_part,nb_year,lat_mat,lon_mat,coord_mode)
E.write_all(file)
E=echantillon(CAYENNE_SUMMER_25km_025sq(nesting_year,d_min,d_max))
E.gen_positions(nb_part,nb_year,lat_mat,lon_mat,coord_mode)
E.write_all(file)

#E.histogram



print "wrote ",outfile
#print "\n\n Termine"

#fermeture du fichier
file.close()

