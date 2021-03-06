#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import time
import netCDF4 as nc
import matplotlib.pylab as pl
from scipy.interpolate import interp2d
from mpl_toolkits.basemap import Basemap
import operator
import datetime as dt
from matplotlib.patches import Rectangle


# Personal libraries 
import biogeolib as bg

# Declaration of functionsi


# ||==========================================
# ||                                        || 
# ||               Tout venant              ||
# ||                                        ||
# ||                                        ||
# ==========================================||

def height_trapezoid(A, h1, delta):
    """
    Compute height of trapezoid of area A, small side distance h1 
    from circle center, and angle delta.
    Solve a degree 2 polynomial : L=2htan(delta/2) ; h=2A/(l+L)
    """ 
    l = 2 * h1 * np.tan(delta/2)
    d = l**2 + 4 * np.tan(delta/2) * (h1 * l + 4 * A + h1**2 * np.tan(delta/2))
    h = (2 * h1 * np.tan(delta/2) - l + np.sqrt(d)) / (4 * np.tan(delta/2))
    return h - h1

def compute_directions(speed, theta, nb_directions, bins_classes):
    nb_classes = len(bins_classes) - 1
    #flatten data and remove nan
    speed = speed[~np.isnan(speed)]
    theta = theta[~np.isnan(theta)]
    
    #Compute histogram
    hist = np.histogram(theta, bins=nb_directions, range=(-np.pi, np.pi))   

    #Compute mean and std for each bin
    rose = np.zeros((nb_directions, nb_classes))
    for k in range(nb_directions):
        idx = np.where((theta > hist[1][k]) & (theta < hist[1][k+1]))
        sp = speed[idx]
        h = np.histogram(sp, bins=bins_classes)
        rose[k, :] = h[0]
    
    directions = (hist[1][:-1] + hist[1][1:]) / 2
    return rose, directions


def mean_direction(ltheta):
    """
    calculate mean direction of a list (or array) of angles expressed in radians
    """
    #remove nan
    ltheta = ltheta[~np.isnan(ltheta)]
    
    #calculate sum of sin and cos of angles
    S = np.sum(np.sin(ltheta))
    C = np.sum(np.cos(ltheta))

    return np.arctan2(S,C)


def compute_speed(us, vs):
    """
    us and vs are 2 arrays from batch
    """
    vector = np.array([us, vs])
    speed = np.linalg.norm(vector, axis = 0)
    return speed

def shift_coord(lat,lon) :
    """ Shift coordinates so that 36.125<lon<396.375 and -60<lat<66.5 """
    for i in range(len(lat)) :
        if lon[i] < 36.125 :
            lon[i] += 360.
        if lat[i] > 66.5 :
            lat[i] -= 180.
        elif lat[i] < -60 :
            lat[i] += 180.
    return lat,lon
      
def read_nc_index(infile,dict_keys,i0,j0) :
    """ Similar to read_nc, but only loads specified index """
    nc_file=nc.Dataset(infile, 'r')
    data = {}.fromkeys(dict_keys)
    for key in data.keys():
        data = (np.array(nc_file.variables[key][0,j0,i0]))
    nc_file.close()
    return data

def geo_to_grid(lon, lat, lon_mat, lat_mat):
    """computes grid coordinates (x,y) of the closest cell center from a (lon,lat) geographical point.
    - lon,lat are the geographical coordinate of the point.
    - lon_mat, lat are matrices containing the longitudes (resp. latitudes) of the centers of the cells of the Arakawa C grid.
    NB : this function is only appropriate for regular, square grids.
    """
    # Read extreme coordinates
    lon_max = np.max(lon_mat)
    lon_min = np.min(lon_mat)
    lat_max = np.max(lat_mat)
    lat_min = np.min(lat_mat)
    h_x = (lon_max-lon_min)/len(lon_mat)
    h_y = (lat_max-lat_min)/len(lat_mat)
    # Converts to grid coordinate
    if ((lon> lon_max) or (lon< lon_min) or (lat>lat_max) or (lat<lat_min)):
        print("\n Longitude or latitude out of range, must be between ", lon_min, " & ", lon_max, " for longitude and ", lat_min, " & ", lat_max, " for latitude \n")
    elif lon < lon_mat[0]:
        x = (360 + lon - lon_mat[0])/h_x      
        y = (lat - lat_min)/h_y 
    else:
        x=(lon-lon_min)/h_x
        y=(lat-lat_min)/h_y

    return x,y
 



def fx_inv(lon,lat,lon_mat,lat_mat, integer = False):
    """
    This function gets the x,y position on grid corresponding to longitude and latitude
        lon : f, longitude
        lat : f, latitude
        lon_mat : matrix containing the value of the longitude at each grid point
        lat_mat : matrix containing the value of the latitude at each grid point
        integer : True to return integers
    """

    lon_min = np.min(lon_mat)
    lon_max = np.max(lon_mat)
    lat_min = np.min(lat_mat)
    lat_max = np.max(lat_mat)

    # Compute position of nearest point on grid
    if ((lon> lon_max) or (lon< lon_min) or (lat>lat_max) or (lat<lat_min)):
        print("Longitude or latitude out of range, must be between ", np.min(lon_mat), " & ", np.max(lon_mat), " for longitude and ", np.min(lat_mat), " & ", np.max(lat_mat), " for latitude")
        return 
    else:
        distance = (lon - lon_mat)**2  + (lat - lat_mat)**2
        inds = np.argmin(distance)
        #argmin compute index of the flattened array (array.flat()), unravel_index gives the corresponding indexes of the 2D array
        i2, i1 = np.unravel_index(inds, distance.shape)
        if integer:
            return i1, i2

    # Compute position of 4 neighbours on grid
        #   nw ---- ne  i2+1
        #   |        |
        #   |        |
        #   |        |
        #   sw ---- se  i2
        #   i1    i1+1

    if (lon>=lon_mat[i2,i1]):
        if (lat>=lat_mat[i2,i1]):
            sw = i2,i1
            nw = i2+1,i1
            se = i2,i1+1
            ne = i2+1,i1+1
        else:
            nw = i2,i1
            sw = i2-1,i1
            se = i2-1,i1+1
            ne = i2,i1+1
    else:
        if (lat>=lat_mat[i2,i1]):
            sw = i2,i1-1
            nw = i2+1,i1-1
            se = i2,i1
            ne = i2+1,i1
        else:
            nw = i2,i1-1
            sw = i2-1,i1-1
            se = i2-1,i1
            ne = i2,i1

    a = lon - lon_mat[sw]
    b = lat - lat_mat[sw]
    xsize = lon_mat[se] - lon_mat[sw]
    ysize = lat_mat[nw] - lat_mat[sw]
        #Weighted barycentre of the 4 neighbour points (cell can be trapezoidal in Orca)
    i = (xsize-a)/xsize*(ysize-b)/ysize*sw[1] + (xsize-a)/xsize*b/ysize*nw[1] + a/xsize*(ysize-b)/ysize*se[1] + a/xsize*b/ysize*ne[1]+1
    j = (xsize-a)/xsize*(ysize-b)/ysize*sw[0] + (xsize-a)/xsize*b/ysize*nw[0] + a/xsize*(ysize-b)/ysize*se[0] + a/xsize*b/ysize*ne[0]+1
    return i,j



def grid_to_geo(x, y, lon_mat, lat_mat) : 
    # Read extreme coordinates
    lon_max=np.max(lon_mat)
    lon_min=np.min(lon_mat)
    lat_max=np.max(lat_mat)
    lat_min=np.min(lat_mat)
    # Converts to grid coordinate
    h_x=(lon_max-lon_min)/len(lon_mat)
    h_y=(lat_max-lat_min)/len(lat_mat)    
    lon=h_x * x + lon_mat[0]
    lat=h_y * y + lat_min
    if lon < lon_mat[0]:
        print('recalibrage')
        lon = (h_x * x) - 360 + lon_mat[0]
        lat = (h_y * y) + lat_min
    
    return lon, lat
    
    
def rotation(x,y,Ox,Oy,theta):
    """ """
    X=Ox+(x-Ox)*np.cos(theta)+(y-Oy)*np.sin(theta)
    Y=Oy-(x-Ox)*np.sin(theta)+(y-Oy)*np.cos(theta)
    return(X,Y)

def Lagrange(x, x_values, y_values):
    """ """
    def _basis(j):
        p = [(x - x_values[m])/(x_values[j] - x_values[m]) for m in xrange(k) if m != j]
        return reduce(operator.mul, p)
    assert len(x_values) != 0 and (len(x_values) == len(y_values)), 'x and y cannot be empty and must have the same length'
    k = len(x_values)
    return sum(_basis(j)*y_values[j] for j in xrange(k))

def plot_Lagrange(X, x_values, y_values) :
    """ """
    Y = np.zeros(len(X))
    i = 0
    for x in X :
        Y[i] = Lagrange(x, x_values, y_values)
        i+=1
    pl.plot(X,Y)
    return Y



def read_tracking_file(file_name) :
    """ Convert tracking data in .dat format to python array. """    
    data = {}.fromkeys(['date','day','month','year','lat', 'lon'])
    for k in data.keys() :
          data[k] = []
    File = open(file_name,'rb')
    for row in File:
        values = list(row.strip().split(' '))
        # data['date'].append(float(values[0]))
        data['day'].append(int(values[2]))
        data['month'].append(int(values[1]))
        data['year'].append(int(values[0]))
        data['lon'].append(-float(values[3]))
        data['lat'].append(float(values[4]))
    return data


def read_csv_file(file_name,variables,types,separator,header) :
    """ Convert CSV file to Python array. """
    data = {}.fromkeys(variables)
    for k in data.keys() :
          data[k] = []
    File = open(file_name,'rb')
    i = 0
    for row in File :
        if i != 0 or header == False:
            values = list(row.strip().split(separator))
            j = 0
            test_NA = 1
            k = 0
            for variable1 in variables :
                test_NA *= (values[k] != 'NA')
                k+=1
            for variable in variables :
                if test_NA == 1 :
                    data[variable].append(types[j](values[j]))
                j+=1
        i += 1
    return data


    
def data_to_same_calendar(vars, dates) :
    """ Sample tracking data to daily values and set them to the same calendar
    that each line of lat/lon matrices correspond to the same date."""
    nturtle = len(vars.values()[0])
    new_dates = []

    # Resample date to 1 day time step.
    for turtle in range(nturtle) :
        new_dates.append(dates[turtle][0::8])
    dates = new_dates[:]

    # Build common daily calendar.
    start_date = min([min(dates[turtle]) for turtle in range(len(dates))])
    end_date  = max([max(dates[turtle][:]) for turtle in range(len(dates))])
    numdays = (end_date - start_date).days
    newdates = [start_date + dt.timedelta(days = x) for x in range(0,numdays)]

    # Compute daily mean positions and build position matrices.
    mean = {}  
    new = {}
    for turtle in range(nturtle) :
        for var in vars.keys() :
            if mean.has_key(var) == False :
                mean[var] = []
            # Compute daily mean position.
            mean[var].append([np.mean(vars[var][turtle][i:i+8]) for i in range(len(vars[var][turtle])/8)])
            # Build position matrix.
            start_date = dates[turtle][0]
            delay = (start_date - newdates[0]).days
            before = [float('nan') for i in range(delay)]
            after = [float('nan') for i in range(numdays-(delay + len(mean[var][turtle][:])))]
            if new.has_key(var) == False :
                new[var] = []
            new[var].append(before + mean[var][turtle] + after)
            new[var].append(before + mean[var][turtle] + after)
    return new, dates
           


def geod_dist(A,B) :
    """Computes geodesic distance (in m) between two points of coordinates 
    A = (latA,lonA) and B = (latB,lonB)""" 
    # Convert angles to radian
    A = np.asarray(A)
    B = np.asarray(B)
    A=A*np.pi/180.
    B=B*np.pi/180.
    # Assign explicit names and load earth radius R
    latA = A[0]
    lonA = A[1]
    latB = B[0]
    lonB = B[1]  
    R = 6378*1000
    # Compute distance
    a = np.arccos(np.sin(latA)*np.sin(latB) + np.cos(latA)*np.cos(latB)*np.cos(lonB-lonA))
    return R*a

def plan_tangeant(A,B) :
    """ """
    R = 6371 # Rayon terrestre
    # conversion en radian
    A=np.asarray(A) 
    A=A*np.pi/180
    B=np.asarray(B)
    B=B*np.pi/180
    # Correction 
    lat_moy = (B[0]+A[0])/2
    X = (A[1]-B[1])*R*np.cos(lat_moy)
    Y = (A[0]-B[0])*R
    return np.sqrt(X**2+Y**2)*1000.

def linear_regression(x,y) :
    """ """
    A = np.array([ x, np.ones(len(x))])
    # linearly generated sequence
    w = np.linalg.lstsq(A.T,y)[0] # obtaining the parameters
    # plotting the line
    line = w[0]*x+w[1] # regression line
    err=np.sqrt(sum((line-y)**2)/len(x))
    return w[1],w[0],err

def read_current(data_path,daynum) :
    """ """
    t = time.time()
    U_file = data_path + "U/GLORYS2_"+str(daynum)+"_gridU.nc"
    V_file = data_path + "V/GLORYS2_"+str(daynum)+"_gridV.nc"
    U_data = read_nc(U_file,['x','y','vozocrtx'])
    V_data = read_nc(V_file,['x','y','vomecrty'])
    xu = np.asarray(U_data['x'])
    yu = np.asarray(U_data['y'])
    xv = np.asarray(V_data['x'])
    yv = np.asarray(V_data['y'])
    U_cur = U_data['vozocrtx'][0,:,:]
    V_cur = V_data['vomecrty'][0,:,:]
    return U_cur,V_cur,xu,xv,yu,yv

def interpolate_var(var, xgrid , ygrid, lon, lat) :
    """ """
    t = time.time()
    i0,j0 = bg.geo_to_grid(lon,lat,xgrid,ygrid)
    if lon > xgrid[i0] :
        i1 = i0 + 1
    else :
        i1 = i0 -1
    if lat > ygrid[j0] :
        j1 = j0 + 1
    else :
        j1 = j0 - 1
    var0 = var[j0,i0]
    var1 = var[j0,i1]
    var2 = var[j1,i1]
    var3 = var[j1,i0]
    X = [xgrid[i0],xgrid[i1],xgrid[i1],xgrid[i0]]
    Y = [ygrid[j0],ygrid[j0],ygrid[j1],ygrid[j1]]
    Var = [var0,var1,var2,var3]
    var_interp = interp2d(X,Y,Var)
    var = var_interp(lon,lat)[0]
    elapsed = time.time() - t
    return var

#----------------------------------------------------------------------------- 
   
# ||==========================================
# ||                                        || 
# ||           Fonctions de tracé           ||
# ||                                        ||
# ||                                        ||
# ==========================================||




def show_control_zone(latlim, lonlim, color = '#8dd3c7') :
    """Plot the rectangle corresponding to latlim and lonlim """
    # Show control zone
    latmin = min(latlim)
    latmax = max(latlim)
    lonmin = min(lonlim)
    lonmax = max(lonlim)
    currentAxis = pl.gca()
    currentAxis.add_patch(Rectangle((lonmin,latmin), lonmax-lonmin, 
                          latmax-latmin, facecolor=color, linestyle="dashed",
                          alpha=0.5))

def plot_current(GLORYS_path,numday) :
    """ """
    # Read data.
    U_cur,V_cur,xu,xv,yu,yv = read_current(GLORYS_path,numday)
    # Interpolate V over U grid :
    V_interp = np.zeros(np.shape(V_cur))
    for i in range(np.shape(V_interp)[0]) :
        for j in range(np.shape(V_interp)[1]) :
            V_interp[i,j] = interpolate_var(V_cur, xv , yv, xu[i], yu[j])
    xu,U_cur = shift_lon(xu,U_cur)  
    xv,V_interp = shift_lon(xv,V_interp)
    pl.quiver(xu,yu,U_cur,V_interp)
    pl.show() 
    
def shift_lon(lon,M) :
    """ """
    lon1 = lon.copy()
    i = 0
    while lon[i]<360.25 :
        i+=1
    i0 = i
    lon1[i0:] -= 360.
    lon1 = np.sort(lon1)
    inf = M[:,i0:]
    M1 = np.concatenate((inf,M[:,0:i0]), axis=1)
    return lon1,M1

#----------------------------------------------------------------------------- 
   
# ||==========================================
# ||                                        || 
# ||          Fonctions de tortues          ||
# ||                                        ||
# ||                                        ||
# ==========================================||



def t_hab(T,SCL,To) :
    """ Compute temperature habitat """
    T_hab = np.ones([np.shape(T)[0],np.shape(T)[1]])
    ########################
    #Tasym
    ########################
    Tmin = 24.4-7.9*SCL 
    sigma = (To-Tmin)/2
    #inf = np.where(T < To)
    #sup = np.where(T>=To)
    #T_hab = np.ones([np.shape(T)[0],np.shape(T)[1]])
    #T_hab[inf] =  np.exp(-(T[inf]-To)**2/(2*sigma**2))
    T_hab =  np.exp(-(T-To)**2/(2*sigma**2))
    #T_hab[sup] = 1.
    #########################
    #VarTasym
    #########################
    """Mass = 0.000214*(SCL*100)**2.86
    
    #Topt_newborn = 24.4
    #Topt_grownup = 17.
    #Tmin_newborn = 23.
    #Tmin_grownup = 5.
    #Topt = ((Topt_grownup - Topt_newborn)/(312.))*Mass + Topt_newborn
    #Tmin = ((Tmin_grownup - Tmin_newborn)/(312.))*Mass + Tmin_newborn
    
    #Topt = 24.4-(0.464*np.sqrt(Mass))
    #Tmin = Topt - (0.35*((Mass)**(2./3.)))
    
    Topt = 24.4 - 0.21*np.sqrt(Mass)
    Tmin = 24.4 - 0.84*np.sqrt(Mass)
    sigma = (Topt-Tmin)/2

    #inf = np.where(T < Tmin)
    #mid = np.where((Tmin<=T) & (T<=Topt))
    mid = np.where(T<=Topt)
    sup = np.where(T > Topt)
    #T_hab[inf]=0.0
    #T_hab[mid]=0.5 * (1-np.sin((np.pi*(Tmin+Topt-2*T[mid]))/(-2*(Tmin-Topt))))
    T_hab[mid]=np.exp(-(T[mid]-Topt)**2/(2*sigma**2))
    T_hab[sup]=1.0"""
    
    '''
    Tbo = 24.4	#Optimal temperature (nesting beach) from Gaspar et al 2012	
    Tmin = Tbo-7.9*SCL     #Critical minimum temperature from Gaspar et al 2012
    mu = Tbo
    sigma = (Tbo-Tmin)/2
    if T <= Tbo:
        T_hab =  np.exp(-(T-Tbo)**2/(2*sigma**2))
    else :    
        T_hab = 1.
    '''
    #T_hab = np.exp((-(T-To)**2)/2/sigma**2)
    #Tmin = 14 - 2/(1.088 - 0.038)*SCL
    return T_hab

def food_hab(mnk,mnk_max) :
    """ Compute food habitat"""
    Food_hab = mnk/mnk_max
    Food_hab[Food_hab>1] = 1
    return Food_hab

def age_to_SCL(age) :
    """ Compute SCL for a given age using a model proposed in F. Perham et
    al., Age and growth of Loggerhead sea turtle of coastal Georgia, 1997"""
    #A = 1.088
    #B = 0.9649
    #k = 0.0739
    #SCL = A*(1 - B*np.exp(-k*age/365.))
    SCL = 1.43*(1-np.exp(-0.226*(age/365.+0.17)))
    return SCL

def compute_M(SCL):
    """
    Compute all turtles mass(kg)

    Ref : Jones, T.T., Hastings, M.D., Bostrom, B.L., Pauly, D., Jones, D.R., 2011. Growth of captive leatherback turtles, Dermochelys coriacea, with inferences on growth in the wild: Implications for population decline and recovery. Journal of Experimental Marine Biology and Ecology 399, 84–92.
    """
    M = 0.000214*(SCL*100)**2.86
    return M
    
def compute_PPmax(age):
    """
    Compute food threshold
    Ref : Jones, T.T., Bostrom, B.L.,  Hastings, M.D., Van Houtan K.S., Pauly, D., 2012. Resource requirements of the Pacific Leatherback Turtle Population
    """
    #Besoin annuels en valeur absolue
    PPmax = 312*2.86*0.299*(((1-np.exp(-0.299*(age/365.+0.17)))**(2.86-1))*(np.exp(-0.299*(age/365.+0.17))))/(1-(1-np.exp(-0.299*(age/365.+0.17)))**(2.86*0.0328))
    #Ramenés à une valeur de production primaire maximale
    PPmax = PPmax/56.8
    return PPmax
def vmax(SCL) :
    """ Compute maximum sustainable speed"""
    #A = 0.26
    #B = -1.12
    #Vmax = SCL * A * SCL**B # A et B sont obtenus par reeression lineaire à 
                            # partir des données de M. Abecassis.
    
    Vmax = 3.5e4*2.7*SCL/86400
    return Vmax

def bathy_hab(bathy,dmin,dmax) :
    """ """
    Pred_hab = np.zeros(np.shape(bathy))
    bathy = np.asarray(-bathy)
    alpha = 1/(dmax-dmin)
    beta = - dmin/(dmax-dmin)
    depth_pref = alpha * bathy + beta
    for i in range(np.shape(Pred_hab)[0]) :
        for j in range(np.shape(Pred_hab)[1]) :
            Pred_hab[i,j] = max(min(depth_pref[i,j],1),0)
    return Pred_hab



















