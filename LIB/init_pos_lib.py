#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 17 11:31:26 2019

@author: tcandela
"""

"""
Librairies to read nesting beach features and to compute initial positions
"""
# =============================================================================
# IMPORTS
# =============================================================================
from mpl_toolkits.basemap import Basemap
import numpy as np
import datetime
import random as rd
import matplotlib.pyplot as plt
import json
from shapely.geometry import Point, Polygon
import math

#personal libraries
import librumeau as brum
# =============================================================================
# FUNCTIONS
# =============================================================================
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
                release_zone_size = 0.25,
                d_min = 40 - (111*0.25)/2,
                d_max = 40 + (111*0.25)/2,
                ref="no reference",
                beach_name="no name",
                season_type="nesting",
		nb_turtles= 1000,
        release_mode = 'square'):

        self.time=time
        self.north_pt=north_pt
        self.south_pt=south_pt
        self.d_min=d_min
        self.d_max=d_max
        self.nb_turtles=nb_turtles
        self.ref=str(ref)
        self.beach_name=str(beach_name)
        self.season_type=str(season_type)
        self.release_mode = str(release_mode)

        if self.season_type=="nesting":
            self.correction_incubation()
    
    def correction_incubation(self):
        """
        add 60 days to the nesting season for the incubation time
        """
        self.time.begin+=60.0
        self.time.end+=60.0
        self.time.max+=60.0
        
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
    def gen_positions(self, nb_part, nb_year, lat_mat, lon_mat, mask, coord_mode, grid_type, beach_orientation, disk_radius, width):
        """
        generate the sample
        the number of particles and number of seasons are given as arguments
        """
        self.nb_part=nb_year*int(nb_part/nb_year)
        self.nb_year=nb_year 
        self.beach_orientation = beach_orientation
        self.random_time()
        for j in range(self.nb_part):
            self.random_position(self.vect_t[j], lat_mat, lon_mat, mask, coord_mode, grid_type, beach_orientation, disk_radius, width)
            print("Turtle %s of %d"%(j+1,nb_part), end='\r')



    def gen_positions_all_year(self,nb_part,nb_year,lat_mat,lon_mat):
        """
        generate the sample
        the number of particles and number of seasons are given as arguments
        """
        self.nb_year=nb_year    
        
        h=61.*np.float32(nb_year)/np.float32(nb_part) #anna 61 ???

        self.vect_t=np.arange(0.5,61.*nb_year+0.5,h)
        self.nb_part=self.vect_t.shape[0]

        for j in range(self.nb_part):
            self.random_position(self.vect_t[j], lat_mat, lon_mat, grid)


    def write_all(self,file):
        """
        write the sample in a given file (already open)
        """
        for j in range(self.nb_part):
            self.positions[j].write_line(file)

    def histogram(self, filename):
        """
        display in a new figure the histogram of the sample
        as time index, one unit correspond to sampling time of the ORCA025 model used
        for exemple, for 6day_mean fields, one year is approximated to 61units=366j
        """
        for j in range(self.nb_part):
            plt.style.use('ggplot')
            fig = plt.figure(num=1, figsize=(12,6), dpi=100, facecolor='w', edgecolor='w')
            ax = fig.add_subplot(111)
            plt.subplots_adjust(left=0.1, bottom=0.2, right=0.95, top=0.95, wspace=None, hspace=None)
            n, bins, patches = ax.hist(self.vect_t, bins=np.arange(self.beach.time.begin, self.beach.time.end, 1))
            ax.set_xlabel('Hatchling period (julian day)')
            ax.set_ylabel('Number of turtles')
            plt.show()
            #plt.savefig(filename)
        
    def release_map(self,filename,lon_mat,lat_mat):
        fig=plt.figure(facecolor='w')
        fig.set_size_inches(12,6)
        plt.title('Release Map',Fontsize=18,Fontweight='bold')
        x=np.zeros(np.size(self.positions))
        y=np.zeros(np.size(self.positions))
        for k in range(np.size(self.positions)):
            x[k] = self.positions[k].x
            y[k] = self.positions[k].y
            
        lon, lat = brum.grid_to_geo(x, y, lon_mat, lat_mat)
        distance = 2.5 #°from beach
        map = Basemap(llcrnrlon=self.beach.north_pt[1] - distance,llcrnrlat=self.beach.north_pt[0]-distance,urcrnrlon=self.beach.north_pt[1]+distance,urcrnrlat=self.beach.north_pt[0]+distance,projection='cyl',resolution='h')
        map.drawcountries()
        map.fillcontinents(color='grey')
        map.drawcoastlines()
        parallels = np.arange(-90,90,2)
        map.drawparallels(parallels,labels=[True,False,True,False],linewidth=0.5,fontsize=16,dashes=[1, 4])
        meridians = np.arange(-180,180,2)
        map.drawmeridians(meridians,labels=[True,False,False,True],linewidth=0.5,fontsize=16,dashes=[1, 4])
        map.scatter(lon, lat, c='b', marker='.', linewidth=0)
        plt.show()
        #plt.savefig(filename)

    ## private 
    # time randomisation
    def random_time(self):
        """matplotlib
        ://en.wikipedia.org/wiki/Pan_(genus)
        compute the set of intial times        
        create a sorted list of times self.vect_t
        
        called by : gen_position
        call to : echantillon_normal_cut
        should not be called directly
        """
        #print("\n=> Computing the set of departure times")
        tt=[]

        counter=0
        for year in range(self.nb_year):

            #tt=tt+self.echantillon_normal_cut(year)
            tt=tt+self.random_time_uniform(year)
#            tt=tt+self.echantillon_normal_asym(year)
        print('random_time_uniform')
        print("\n=> "+str(self.nb_part)+' initial positions computed\n')
        tt=np.array(tt,dtype=np.float32)
        tt.sort()
        self.vect_t=tt

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
    
    def random_time_uniform(self,year):
        m=self.beach.time.begin
        M=self.beach.time.end
        h=[]
        Ntot=self.nb_part/self.nb_year
        tur_per_day = np.round(Ntot/((M+0.5) - (m + 1.5)))
        for day in np.arange(m+0.5,M+1.5):
            h.append(tur_per_day)
        h = np.array(h)
        S = sum(h)
        a = Ntot/S
        n_days = np.round(a*h)
        Ntot_bis= np.sum(np.round(a*h))
        err_round = Ntot-Ntot_bis
        add = err_round/len(n_days)
        if add < 1:
            add = 1
        else :
            add = int(add)
        for i in np.arange(len(n_days)):
            if np.sum(n_days) != Ntot:
                n_days[i] += add
        time=np.ones([int(Ntot)])*(m+0.5)
        part_hatched=0
        for i in range(np.shape(n_days)[0]):
            time[int(part_hatched):int(part_hatched+n_days[i])]=(np.float(m+i+0.5))+365.25*year
            part_hatched = part_hatched+n_days[i]

        #plt.hist(time,bins=np.arange(m,M,1))
        #plt.show()


        return list(time) 
    
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
        h=[]
        Ntot=self.nb_part/self.nb_year
        for day in np.arange(m+0.5,M+1.5):
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
        n_days[np.where(n_days==np.max(n_days))[0][0]] += err_round
        time=np.ones([int(Ntot)])*(m+0.5)
        part_hatched=0
        for i in range(np.shape(n_days)[0]):
            time[int(part_hatched):int(part_hatched+n_days[i])]=(np.float(m+i+0.5))+365.25*year
            part_hatched = part_hatched+n_days[i]

        #plt.hist(time,bins=np.arange(m,M,1))
        #plt.show()


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
            ####################
            #Uniform ditrib
            ####################
            #a = rd.uniform(m,M)

            time[j]=np.floor(a+365.25*year)+0.5
        #plt.hist(time,bins=np.arange(m,M,1))
        #plt.show()
        return time

    # space randomisation
  
    def random_position(self, tt, lat_mat, lon_mat, mask, coord_mode,grid_type, beach_orientation, disk_radius, width):
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
#        xA=self.beach.north_pt[1]*np.pi/180
#        yA=self.beach.north_pt[0]*np.pi/180
#        xB=self.beach.south_pt[1]*np.pi/180
#        yB=self.beach.south_pt[0]*np.pi/180
#        Ox=xB
#        Oy=yB
#        theta=-np.arctan((xB-xA)/(yB-yA))
#        
#        #conversion de la translation dans le repère géographique
#        xT_min,yT_min=brum.rotation(self.beach.d_min,0,Ox,Oy,-theta)
#        xT_max,yT_max=brum.rotation(self.beach.d_max,0,Ox,Oy,-theta)
#        
#        #conversion des distances en angles (approx. du plan tangent)
#        lat_moy=(yA+yB)/2
#        R_terre=6371
#        xT_min=xT_min/R_terre/abs(np.cos(lat_moy))
#        xT_max=xT_max/R_terre/abs(np.cos(lat_moy))
#        yT_min=yT_min/R_terre
#        yT_max=yT_max/R_terre
#        XT_min=np.sqrt(xT_min**2+yT_min**2)
#        XT_max=np.sqrt(xT_max**2+yT_max**2)
#        YT_min=0
#        YT_max=0
#    
#        #conversion vers le repère rectifié :
#        XB,YB=brum.rotation(xB,yB,Ox,Oy,theta)
#        XA,YA=brum.rotation(xA,yA,Ox,Oy,theta)
#        XC=XA+XT_min
#        XD=XA+XT_max
#        XF=XB+XT_min
#        XE=XB+XT_max
#        YC=YA+YT_min
#        YD=YA+YT_max
#        YF=YB+YT_min
#        YE=YB+YT_max
        if self.beach.release_mode == 'disk':
            #In disk mode, the center of the disk is (lon_0;lat0).
            x0=self.beach.north_pt[1]
            y0=self.beach.north_pt[0]
            #
            r = disk_radius * math.sqrt(rd.random())
            theta = rd.random() * 2 * math.pi
            #
            X = x0 + r * math.cos(theta) #/ (111195 * math.cos(y0*math.pi/180)) uncomment if disk_radius in meters
            Y = y0 + r * math.sin(theta) #/ 111195
            #
            if grid_type == 'orca':
                x_grid,y_grid=brum.fx_inv(X,Y,lon_mat,lat_mat)
            else:
                x_grid,y_grid=brum.geo_to_grid(X,Y,lon_mat,lat_mat)
       
        # Beach release
        if self.beach.release_mode == 'beach':
            # suppose que lon_mat et lat_mat sont 1D
            # attention il peut y avoir des pb des points sur la terre car on ne regarde pas le masque ici...
            xA=self.beach.north_pt[1]
            yA=self.beach.north_pt[0]
            xB=self.beach.south_pt[1]
            yB=self.beach.south_pt[0]
            
            
            # angle entre (AB) et l'est (positif de 0 à pi)
            theta = np.arctan((yA-yB)/(xB-xA)) if xB > xA else np.pi - np.arctan((yA-yB)/(xA-xB))

            # find horizontal edge towards east or west (beach otientation)
            if theta < np.pi/4 or theta > 3 * np.pi / 4:
                if beach_orientation == 'E':
                    a=2 # A implémenter...
             
            # find vertical edge towards east or west (beach orientation)
            else:
                #point de grille extérieur à A
                XA = list(lon_mat > xA).index(True)
                XB = list(lon_mat > xB).index(True)
                
                # point de grille extérieur à B
                YA = list(lat_mat > yA).index(True)
                YB = list(lat_mat > yB).index(True) - 1
                
                if beach_orientation == 'W':
                    # on prend l'indice précédent
                    XA -= 1
                    XB -= 1
                
                Y = rd.uniform(yB, yA) 
                if XA == XB:
                    X = lon_mat[XA]
                else:
                    a=2# à impléménter...

            if grid_type == 'orca':
                x_grid,y_grid=brum.fx_inv(X,Y,lon_mat,lat_mat)
            else:
                x_grid,y_grid=brum.geo_to_grid(X,Y,lon_mat,lat_mat)


        if self.beach.release_mode == 'rectangle':
            xB=self.beach.north_pt[1]#*np.pi/180
            yB=self.beach.north_pt[0]#*np.pi/180
            xA=self.beach.south_pt[1]#*np.pi/180
            yA=self.beach.south_pt[0]#*np.pi/180
            
            theta=np.arctan((yB-yA)/(xB-xA))
            lat_moy=(yA+yB)/2 * np.pi/180
            AB = np.sqrt(((xB-xA)*111.195*np.cos(lat_moy))**2 + ((yB-yA)*111.195)**2) #km
            l = rd.uniform(-1, 1) * width/2
            d = rd.random() * AB 
            xC = xB + d * np.cos(theta) / (111.195*np.cos(lat_moy))
            yC = yB + d * np.sin(theta) / 111.195
            
            Y = yC + l * np.sin(theta) / 111.195

            if xB > xA:
                if l <= 0:
                    X = xC - l * np.sin(theta) / (111.195*np.cos(lat_moy))
                elif l > 0:
                    X = xC + l * np.sin(theta) / (111.195*np.cos(lat_moy))

            elif xB < xA:
                if l <= 0:
                    X = xC + l * np.sin(theta) / (111.195*np.cos(lat_moy))
                elif l > 0:
                    X = xC - l * np.sin(theta) / (111.195*np.cos(lat_moy))

            if grid_type == 'orca':
                x_grid,y_grid=brum.fx_inv(X,Y,lon_mat,lat_mat)
            else:
                x_grid,y_grid=brum.geo_to_grid(X,Y,lon_mat,lat_mat)           
            
        if self.beach.release_mode == 'square':
            xA=self.beach.north_pt[1]*np.pi/180
            yA=self.beach.north_pt[0]*np.pi/180
            xB=self.beach.south_pt[1]*np.pi/180
            yB=self.beach.south_pt[0]*np.pi/180
            theta=-np.arctan((xB-xA)/(yB-yA))
            
            xT_min = np.cos(theta) * self.beach.d_min
            yT_min = np.sin(theta) * self.beach.d_min
            xT_max = np.cos(theta) * self.beach.d_max
            yT_max = np.sin(theta) * self.beach.d_max
            
            lat_moy=(yA+yB)/2
            xT_min = xT_min/(111.195*np.cos(lat_moy))
            yT_min = yT_min/111.195
            xT_max = xT_max/(111.195*np.cos(lat_moy))
            yT_max = yT_max/111.195
            
            if beach_orientation == 'E':
                XC=(xA*180/np.pi)+xT_min
                XD=(xA*180/np.pi)+xT_max
                XF=(xB*180/np.pi)+xT_min
                XE=(xB*180/np.pi)+xT_max
                YC=(yA*180/np.pi)+yT_min
                YD=(yA*180/np.pi)+yT_max
                YF=(yB*180/np.pi)+yT_min
                YE=(yB*180/np.pi)+yT_max
            if beach_orientation == 'W':
                XC=(xA*180/np.pi)-xT_min
                XD=(xA*180/np.pi)-xT_max
                XF=(xB*180/np.pi)-xT_min
                XE=(xB*180/np.pi)-xT_max
                YC=(yA*180/np.pi)-yT_min
                YD=(yA*180/np.pi)-yT_max
                YF=(yB*180/np.pi)-yT_min
                YE=(yB*180/np.pi)-yT_max        
            
            
            
            poly_coord = [(XC, YC), (XD, YD), (XE, YE), (XF, YF)]
            poly = Polygon(poly_coord) 
#        print(XC,YC,XD,YD,XF,YF, XE,YE)        
        
        #Tirage aléatoire de coordonnées (X,Y) dans le repère rectifié :

        #Test pour ne pas tirer de positions sur la terre
   
            x_grid = 0
            y_grid = 0
            k = 0
    
            while mask[int(y_grid),int(x_grid)]==0 or k ==0:
                X = rd.uniform(min([XC,XD,XE,XF]), max([XC,XD,XE,XF]))
                Y = rd.uniform(min([YC,YD,YE,YF]), max([YC,YD,YE,YF]))
                p1 = Point(X, Y)
                while poly.contains(p1) == False:
                    X = rd.uniform(min([XC,XD,XE,XF]), max([XC,XD,XE,XF]))
                    Y = rd.uniform(min([YC,YD,YE,YF]), max([YC,YD,YE,YF]))
                    p1 = Point(X, Y)
                                           
                # Conversion vers le repère géographique:
    #            xC,yC=brum.rotation(XC,YC,Ox,Oy,-theta)
    #            xF,yF=brum.rotation(XF,YF,Ox,Oy,-theta)
    #            x,y=brum.rotation(X,Y,Ox,Oy,-theta)
            
                #conversion en degrés
    #            x=x*180/np.pi
    #            y=y*180/np.pi
    
                #conversion en coordonnées de grille
                if grid_type == 'orca':
                    x_grid,y_grid=brum.fx_inv(X,Y,lon_mat,lat_mat)
                else:
                    x_grid,y_grid=brum.geo_to_grid(X,Y,lon_mat,lat_mat)
                #print(x,y,x_grid,y_grid)
                k = k+1
            
        if coord_mode == 'grid':
            self.positions=self.positions+[position_4D(x_grid,y_grid,tt)]
        elif coord_mode == 'latlon':
            self.positions=self.positions+[position_4D(X,Y,tt)]


def beach_json(release_zone_size, lon_name, lat_name, beach_carac, nesting_year, d_min, d_max, date_ref):

    #function generalized by anna (with beach input read on a .json file stored in config/beaches/
    #to do add d_min/d_max not default

    start = datetime.datetime.strptime(str(nesting_year)+' '+beach_carac['start_month']+' '+beach_carac['start_day'], '%Y %B %d') 
    print(start, nesting_year)
    start_jday = (start-date_ref).days
    end = datetime.datetime.strptime(str(nesting_year)+'/'+beach_carac['end_month']+'/'+beach_carac['end_day'], '%Y/%B/%d')
    end_jday = (end-date_ref).days
    if end_jday < start_jday:
        end_jday += 365
        
    peak_jday = (start_jday + end_jday)/2. #gaussian
    north_pt=(beach_carac["lat_0"], beach_carac["lon_0"])
    south_pt=(beach_carac["lat_1"], beach_carac["lon_1"])
    try:
        release_mode = beach_carac["release_mode"]
    except:
        release_mode = 'square'
    if release_mode not in ['square', 'disk', 'rectangle', 'beach']:
        raise ValueError('Please set release_mode to square, rectangle, beach or disk')
    
    lon_bc = (north_pt[1]+south_pt[1])/2 #longitude of beach center
    lat_bc = (north_pt[0]+south_pt[0])/2 #latitude of beach center
    if north_pt[1] == south_pt[1]:
        theta = 0
    else:
        theta = np.arctan((north_pt[0]-south_pt[0])/abs(north_pt[1]-south_pt[1]))
    x = (((release_zone_size*np.pi)/180)/2)*(np.cos(theta))
    y = (((release_zone_size*np.pi)/180)/2)*(np.sin(theta))
    x = x*180/np.pi
    y = y*180/np.pi
    
    if release_mode == 'square': 
        if north_pt[1] < lon_bc:
            north_pt = (lat_bc + y, lon_bc - x)
            south_pt = (lat_bc - y, lon_bc + x)
        else:
            north_pt = (lat_bc + y, lon_bc + x)
            south_pt = (lat_bc - y, lon_bc - x)

    return beach(time=nesting_period(start_jday,end_jday,peak_jday),
                 north_pt=north_pt,south_pt=south_pt,
                 ref=beach_carac['ref'],
                 beach_name=beach_carac['beach_name'],
		         nb_turtles=np.int(beach_carac['nb_turtles']),
                 season_type=beach_carac['season_type'],
                 d_min=d_min,d_max=d_max, release_mode = release_mode)

def read_beach_json(file, grid, lon_name, lat_name):
    file_c = open(file)#.read()
    beach_carac = json.loads(file_c.read())

    for k in ["lon_0", "lat_0", "lon_1", "lat_1"]:
        try:
            beach_carac[k] = np.float(beach_carac[k])
        except:
            print('%s not found.'%k)


    for k in ["lon_0", "lon_1"]: # anna attention ne marche qu'en global
        if beach_carac[k] < np.min(grid[lon_name]) :
            beach_carac[k] += 360
        if beach_carac[k] > np.max(grid[lon_name]) :
            beach_carac[k] -= 360

    beach_carac["nb_turtles"] = np.int(beach_carac["nb_turtles"])
    
    file_c.close()

    return beach_carac
