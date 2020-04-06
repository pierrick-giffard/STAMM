#!/usr/bin/env python
#-*- coding:utf-8 -*-
"""
    Author:             Benjamin Rumeau (intern at CLS in 2014 under supervision
                        of Philippe Gaspard and Anne-Cecile Dragon.
    Module description: Library containing all the geographical classes and fun-
                        ctions.
"""
import numpy as np


# ||==========================================
# ||                                        || 
# ||           Objets géographiques         ||
# ||             et biologiques             ||
# ||                                        ||
# ==========================================||


class NestingPeriod :
    """
    class nesting_period used as a structure
    contains the begin, end and peak of the nesting season
    """
    def __init__(self,first=0.,last=0.,peak=0.,incubation_time=60):
        self.begin = np.float32(first)
        self.end = np.float32(last)
        self.peak = np.float32(peak)
        self.incubation_time=incubation_time

    def incubation(self):
        """                                                                  
        Corrects nesting period to take into account incubation time.
        Incubation time must be expressed in days.
        """
        print("\n=> correcting nesting periods with incubation times")
        self.begin+=self.incubation_time
        self.end+=self.incubation_time
        self.peak+=self.incubation_time
        

class Beach :
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
                 nesting_period=NestingPeriod(),
                 north_pt=[1,1],
                 south_pt=[0,0],
                 orientation="east",
                 d_min=30,
                 d_max=50,
                 ref="no reference",
                 beach_name="no name",
                 file_name="no name"):
        self.nesting_period=nesting_period
        self.north_pt=north_pt
        self.south_pt=south_pt
        if orientation=="east":
            self.orientation=1
        elif orientation=="west":
            self.orientation=-1
        else:
            print("Error: orientation must be east or west")
        self.d_min=d_min
        self.d_max=d_max
        self.ref=str(ref)
        self.beach_name=str(beach_name)
        self.nesting_period.incubation()


class Geolocation:
    """ 
    class position_4D contains the a 4 uplet x,y,z,t giving the time ans space
    localisation of a particle, as well as an information about the coordinate
    system.
    """
    def __init__(self,x=0.,y=0.,t=0.,z=-1.0,coord_sys="geo"):
        self.x=np.float32(x)
        self.y=np.float32(y)
        self.z=np.float32(z)
        self.t=np.float32(t)
        self.coord_sys=coord_sys
    
    def geo_to_grid(self,lon_mat,lat_mat):
        """
        converts (lon,lat) geographical coordinate to grid coordinates (conti-
        nuous conversion as coordinates of the closest neighbour are computed
        in Ariane.
        - lon,lat are the geographical coordinate of the point.
        - lon_mat, lat are matrices containing the longitudes (resp. latitudes) of the centers of the cells of the Arakawa C grid.
        NB : this function is only appropriate for regular, square grids.
        """
        if self.coord_sys=="grid":
            print("\n Warning: already in grid coordinates")
        else:
            # Read extreme coordinates
            lon_max=np.max(lon_mat)
            lon_min=np.min(lon_mat)
            lat_max=np.max(lat_mat)
            lat_min=np.min(lat_mat)
            # Converts to grid coordinate
            if ((self.x> lon_max) or (self.x< lon_min) or (self.y>lat_max) or (self.y<lat_min)):
                print("\n Longitude or latitude out of range, must be between ", lon_min, " & ", lon_max, " for longitude and ", lat_min, " & ", lat_max, " for latitude \n")
            else:
                h_x=(lon_max-lon_min)/len(lon_mat)
                h_y=(lat_max-lat_min)/len(lat_mat)
                self.x=(self.x-lon_min)/h_x
                self.y=(self.y-lat_min)/h_y
                self.coord_sys="grid"

    def grid_to_geo(self, lon_mat, lat_mat) : 
        if self.coord_sys=="grid":
            print("\n Warning: already in geographical coordinates")
        # Read extreme coordinates
        lon_max=np.max(lon_mat)
        lon_min=np.min(lon_mat)
        lat_max=np.max(lat_mat)
        lat_min=np.min(lat_mat)
        # Converts to grid coordinate
        h_x=(lon_max-lon_min)/len(lon_mat)
        h_y=(lat_max-lat_min)/len(lat_mat)
        self.x=h_x * self.x + lon_min
        self.y=h_y * self.y + lat_min
        self.coord_sys="grid"


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



# ||==========================================
# ||                                        || 
# ||       Transformations géométriques     ||
# ||                                        ||
# ==========================================||


def rotation(x,y,Ox,Oy,theta):
    """
        Computes coordinates X,Y in referential R' from a poin (x,y) in a referential R
        where R' is obtained by a rotation of R of angle theta and center O(Ox,Oy)
    """      
    X=Ox+(x-Ox)*np.cos(theta)+(y-Oy)*np.sin(theta)
    Y=Oy-(x-Ox)*np.sin(theta)+(y-Oy)*np.cos(theta)
    return(X,Y)







def geo_to_grid(lon,lat,lon_mat,lat_mat):
    """
    converts (lon,lat) geographical coordinate to grid coordinates (conti-
    nuous conversion as coordinates of the closest neighbour are computed
    in Ariane.
    - lon,lat are the geographical coordinate of the point.
    - lon_mat, lat are matrices containing the longitudes (resp. latitudes) of the centers of the cells of the Arakawa C grid.
    NB : this function is only appropriate for regular, square grids.
    """
    # Read extreme coordinates
    lon_max=np.max(lon_mat)
    lon_min=np.min(lon_mat)
    lat_max=np.max(lat_mat)
    lat_min=np.min(lat_mat)
    # Converts to grid coordinate
    if ((lon> lon_max) or (lon< lon_min) or (lat>lat_max) or (lat<lat_min)):
        print("\n Longitude or latitude out of range, must be between ", lon_min, " & ", lon_max, " for longitude and ", lat_min, " & ", lat_max, " for latitude \n")
    else:
        h_x=(lon_max-lon_min)/len(lon_mat)
        h_y=(lat_max-lat_min)/len(lat_mat)
        i1 = int((lon-lon_min)/h_x)
        i2 = i1 + 1
        I = [i1,i2]
        i0 = np.argmin([abs(lon - lon_mat[i1]), abs(lon - lon_mat[i2])])
        i = I[i0]
        j1 = int((lat-lat_min)/h_y)
        j2 = j1 + 1
        J = [j1,j2]
        j0 = np.argmin([abs(lat - lat_mat[j1]), abs(lat - lat_mat[j2])])
        j = J[j0]
    return i,j











































