#!/usr/bin/env python
#-*- coding:utf-8 -*-

"""
        Author : Amaury Dehecq (intern at CLS Summer 2011 working with Phillipe Gaspar)

        Module description : All input/output functions to read :
        - input parameters in namelist file
        - initial positions
        - external variables files (U,V,T...)
        - grid parameters
        and write output file.

        Last modified : 10/08/2012
"""

import numpy as np
import netCDF4 as nc
import csv
import os
import sys


##########
# Inputs #
##########

def read_namelist(filename):
    """
    Read the namelist file and return a dictionnary containing the name of the items (str), and its corresponding value (str)
    read_namelist(filename) -> items
    - filename : str, name of the namelist file
    - items : dictionnary

    For more information on the items, read 'IBM-tutorial' file 
    """

    items = {'init_file':'',
             'type':'',
             'nturtles':'',
             'tstep':'',
             'nsteps_simu':'',
             'species':'',
             'mode':'',
             'alpha':'',
             'vscale':'',
             'Fa':'',
             'key_alltracers':'',
             'nsteps_max':'',
             'key_periodic':'',
             'overlap':'',
             'key_jfold':'',
             'pivot':'',
             'c_dir_zo':'',
             'c_prefix_zo':'',
             'ind0_zo':'',
             'indn_zo':'',
             'maxsize_zo':'',
             'c_suffix_zo':'',
             'nc_var_zo':'',
             'nc_lon_zo':'',
             'nc_lat_zo':'',
             'nc_att_mask_zo':'',
             'c_dir_me':'',
             'c_prefix_me':'',
             'ind0_me':'',
             'indn_me':'',
             'maxsize_me':'',
             'c_suffix_me':'',
             'nc_var_me':'',
             'nc_lon_me':'',
             'nc_lat_me':'',
             'nc_att_mask_me':'',
             'c_dir_te':'',
             'c_prefix_te':'',
             'ind0_te':'',
             'indn_te':'',
             'maxsize_te':'',
             'c_suffix_te':'',
             'nc_var_te':'',
             'nc_lon_te':'',
             'nc_lat_te':'',
             'nc_att_mask_te':'',
             'c_dir_pp':'',
             'c_prefix_pp':'',
             'ind0_pp':'',
             'indn_pp':'',
             'maxsize_pp':'',
             'c_suffix_pp':'',
             'nc_var_pp':'',
             'nc_lon_pp':'',
             'nc_lat_pp':'',
             'nc_att_mask_pp':'',
             'dir_mesh':'',
             'fn_mesh':'',
             'nc_var_xx_tt':'',
             'nc_var_xx_uu':'',
             'nc_var_yy_tt':'',
             'nc_var_yy_vv':'',
             'nc_var_e2u':'',
             'nc_var_e1v':'',
             'nc_var_e1t':'',
             'nc_var_e2t':'',
             'nc_var_tmask':'',
             }

    namelist = open(filename,'r')
    text = csv.reader(namelist,delimiter='=')
    
    print("    ************")
    print("    * NAMELIST *")
    print("    ************")
    for line in text:
        print(str(line).replace('[','').replace(']','').replace('\'','').replace(',',':'))
        for key in items.keys():
            if line[0].__contains__(key):
                items[key]=line[1].replace('\'','').replace(' ','').replace(',','')

    """
    All items are read as strings, they must be converted to correct type
    """
    #Convert integers
    for key in ['nturtles','nsteps_simu','nsteps_max','overlap','ind0_zo','indn_zo','maxsize_zo','ind0_me','indn_me','maxsize_me']:
        try:
            items[key] = int(items[key])
        except ValueError:
            sys.exit("ERROR : %s must be integer" %(key))

    #convert floats
    for key in ['tstep']:
        try:
            items[key] = float(items[key])
        except ValueError:
            sys.exit("ERROR : %s must be float" %(key))

    #convert booleans
    for key in ['key_alltracers','key_periodic','key_jfold']:
        try:
            items[key] = items[key].__contains__('T')
        except ValueError:
            sys.exit("ERROR : %s must be boolean" %(key))

    #Species name to lower caracters
    items['species'] = items['species'].lower()

    #Check mode
    items['mode'] = items['mode'].lower()
    if items['mode'] not in ['passive','active','diffusion']:
        sys.exit("ERROR : mode must be 'passive', 'active' or 'diffusion'")
        
    
    #Check pivot
    if items['pivot']!='T' and items['pivot']!='F':
        sys.exit("ERROR : pivot must be 'T' or 'F'")

    #optional items
    if items['key_alltracers']==True:
        for key in ['ind0_te','indn_te','maxsize_te','ind0_pp','indn_pp','maxsize_pp']:
            try:
                items[key] = int(items[key])
            except ValueError:
                sys.exit("ERROR : %s must be integer" %(key))

    for key in ['alpha','vscale','Fa']:
        if items['mode']=='active':
            try:
                items[key] = float(items[key])
            except ValueError:
                sys.exit("ERROR : %s must be float" %(key))
        else:
            items[key] = 0.


    print("    ****************")
    print("    * END NAMELIST *")
    print("    ****************")

    return items


def fx_inv(lon,lat,lon_mat,lat_mat):
    """
    This function gets the x,y position on grid corresponding to longitude and latitude
        lon : f, longitude
        lat : f, latitude
        lon_mat : matrix containing the value of the longitude at each grid point
        lat_mat : matrix containing the value of the latitude at each grid point
    """
        #Dimensions of the grid
    xdim = lon_mat.shape[1]
    ydim = lon_mat.shape[0]

    lon_min = np.min(lon_mat)
    lon_max = np.max(lon_mat)
    lat_min = np.min(lat_mat)
    lat_max = np.max(lat_mat)

    # Compute position of nearest point on grid
    if ((lon> lon_max) or (lon< lon_min) or (lat>lat_max) or (lat<lat_min)):
        print("Longitude or latitude out of range, must be between ", np.min(lon_mat), " & ", np.max(lon_mat), " for longitude and ", np.min(lat_mat), " & ", np.max(lat_mat), " for latitude")
        return 
    else:
        distance = (lon-lon_mat)**2+(lat-lat_mat)**2
        inds = np.argmin(distance)     
        #argmin compute index of the flattened array (array.flat()), unravel_index gives the corresponding indexes of the 2D array
        i2,i1 = np.unravel_index(inds,distance.shape) 

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


def read_positions(param,glamu,gphiv,mask):
    """
    Function to read initial positions in file.
    - param : dictionnary containing all items of the namelist, output of read_namelist
    - glamu :  matrix containing the value of the longitude at each grid point
    - gphiv :  matrix containing the value of the latitude at each grid point
    Return : np arrays, x_init, y_init, t_init initial position in grid indices and initial times
    """

    nturtles = param['nturtles']

    #Count number of lines and check if there are enough positions 
    try:
        init = open(param['init_file'],'r')
    except IOError:
        sys.exit("Initial positions file does not exist")
    lines = 0
    for line in init:
        lines += 1
    init.close()

    if nturtles>lines:
        print("ERROR - There are not enough initial positions for intended simulation. Add more lines to initial positions file or choose a lower nturtles")
        sys.exit(1)

    x_init = np.zeros(nturtles,dtype='float64')
    y_init = np.zeros(nturtles,dtype='float64')
    t_init = np.zeros(nturtles,dtype='float64')

    #if initial positions are saved as grid indices (x,y)
    if param['type']=='x/y':
        print('****************************************************')
        print("Read initial positions in file : \n ", param['init_file'])
        print('****************************************************')
        init = open(param['init_file'],'r')
        x_init, y_init, t_init = np.loadtxt(init,usecols=(0,1,3),unpack=True)
        x_init, y_init, t_init = x_init[:nturtles], y_init[:nturtles], t_init[:nturtles] 
        init.close()
        
    #if initial positions are saved as lon/lat
    elif param['type']=='lon/lat':
        print('****************************************************')
        print("Read initial positions in file : \n"+ param['init_file'] + "\n and convert to grid points")
        print('****************************************************')
        init = open(param['init_file'],'r')
        #Attention sur l'agencement des colonnes respectivement avec lon/lat
        lat_init, lon_init, t_init = np.loadtxt(init,usecols=(1,0,3),unpack=True)
        lon_init, lat_init, t_init = lon_init[:nturtles], lat_init[:nturtles], t_init[:nturtles]
        print(lon_init[0], lat_init[0])
        for i in range(nturtles):
            x_init[i],y_init[i] = fx_inv(lon_init[i],lat_init[i],glamu,gphiv)
        init.close()

    else:
        print("ERROR - Type must be 'lon/lat' or 'x/y'")
        sys.exit(1)

    #Check that number of days to be simulated does not exceed max number of input files
    if param['nsteps_simu']+np.max(t_init)>param['nsteps_max']:
        print('WARNING - There are not enough input files : Loop on time when last one is reached')
        print('WARNING - max(t_init) = %d ' % np.max(t_init))
        print('WARNING - nsteps_max  = %d ' % (param['nsteps_max']))
        print('WARNING - nsteps_simu = %d, should be less than %d ' % (param['nsteps_simu'],param['nsteps_max']-np.max(t_init)))
    
    #For compatibility with Ariane, we read grid indices as Fortran arrays (start at 1) but then convert to Python (start at 0)
    x_init = x_init-1
    y_init = y_init-1

    #Initial cell cannot be on land
    i0 = np.int32(x_init) + 1
    j0 = np.int32(y_init) + 1

    position_on_mask = mask[j0,i0]
    
    land = np.where(position_on_mask==0)[0]
    
    if len(land)!=0:
        # import pandas as pda
        # dfpos = pda.DataFrame()
        # # i      j           mand t          mand
        # #1081.244 265.187    -1.0  135.50     1.0
        # dfpos["x_init"]=pda.Series()
        # dfpos["x_init"]=x_init
        # dfpos["y_init"]=pda.Series()
        # dfpos["y_init"]=y_init
        # dfpos["mand_1"]=pda.Series()
        # dfpos["mand_1"]="-1.0"
        # dfpos["t_init"]=pda.Series()
        # dfpos["t_init"]=t_init
        # dfpos["mand_2"]=pda.Series()
        # dfpos["mand_2"]="1.0"
        
        # dfpos["i0"]=pda.Series()
        # dfpos["i0"]=i0
        # dfpos["j0"]=pda.Series()
        # dfpos["j0"]=j0
        # dfpos["OnLand"] = dfpos[["i0","j0"]].apply( lambda row: False if mask[row[1],row[0]]==1 else np.nan,axis=1)
        # dfpos.dropna(inplace=True)
        # for col in ["i0","j0","OnLand"]:
        #     del dfpos[col]

        # outfile="./init_pos_on_ocean.txt"
        # dfpos.to_csv(outfile,sep=" ",header=False,index=False)

        print("E R R O R : found positions on land:")
        print(str(land))
        # print("File without position on land (%d positions): %s" %(len(dfpos),outfile))
        # print("Exiting. You can run the code with this init file %s" %(outfile))
        quit()

    return x_init, y_init, t_init
 

##############
# Grid input #       
##############

class Grid(object):
    
    def __init__(self,param):
        """
        Class for grid objects, containing all parameters of the grid :
        - glamt, gphit are matrices containing latitudes and longitudes at each point of the T grid
        - glamu, gphiv : the same for grid U & V
        - e1u, e2v... : matrices containing size in km of grid cell along direction 1 (x), resp. 2 (y) for grid U, V or T
        - mask : matrix containing continental mask (0 for continents, 1 for oceans)
        """
        filename = param['dir_mesh']+param['fn_mesh']

        infile = nc.Dataset(filename)
        e2u = np.squeeze(infile.variables[param['nc_var_e2u']])
        e1v = np.squeeze(infile.variables[param['nc_var_e1v']])
        mask = np.squeeze(infile.variables[param['nc_var_tmask']])
        e1t = np.squeeze(infile.variables[param['nc_var_e1t']])
        e2t = np.squeeze(infile.variables[param['nc_var_e2t']])
        glamt = np.squeeze(infile.variables[param['nc_var_xx_tt']])
        gphit = np.squeeze(infile.variables[param['nc_var_yy_tt']])
        glamu = np.squeeze(infile.variables[param['nc_var_xx_uu']])
        gphiv = np.squeeze(infile.variables[param['nc_var_yy_vv']])
        infile.close()

        self.e2u = np.float64(e2u,order='Fortran')
        self.e1v = np.float64(e1v,order='Fortran')
        self.mask = np.int32(mask,order='Fortran')
        self.tfac = np.float64(e1t*e2t,order='Fortran')
        self.glamt = np.float64(glamt,order='Fortran')
        self.gphit = np.float64(gphit,order='Fortran')
        self.glamu = np.float64(glamu,order='Fortran')
        self.gphiv = np.float64(gphiv,order='Fortran')
        
        self.key_periodic = param['key_periodic']
        self.overlap = param['overlap']
        self.key_jfold = param['key_jfold']
        self.pivot = param['pivot']

        #Mask in the U grid -> all points to the left of a masked point must be masked for U, to avoid particles to enter the mask
        umask = np.zeros_like(self.mask)
        xdim = np.shape(umask)[1]
        #Non overlapping points
        for i in range(0,xdim-1):
            vec1 = self.mask[:,i]
            vec2 = self.mask[:,i+1]
            umask[:,i]= np.where(vec2==0,0,vec1)
        #For the last point, the point to the right is shifted of 'overlap' points 
        vec1 = self.mask[:,-1]
        vec2 = self.mask[:,self.overlap]
        umask[:,-1]= np.where(vec2==0,0,vec1)
    
        self.umask = umask

           #Mask in the V grid -> all points below a masked point must be masked for V, to avoid particles to enter the mask
        vmask = np.zeros_like(self.mask)
        ydim = np.shape(vmask)[0]
        for j in range(0,ydim-1):
            vec1 = self.mask[j,:]
            vec2 = self.mask[j+1,:]
            vmask[j,:]= np.where(vec2==0,0,vec1)
            
        self.vmask = vmask


###################
# Variables input #
###################

class Variable(object):
    
    def __init__(self,dir,prefix,suffix,indsize,var_name,var_lon,var_lat,mask,mask_mat=0,ind=0):
        """
        Read NetCDF files of the form 'dir/prefixXXXXsuffix' with indsize the number of digits XXXX (e.g mydir/GLORYS_0001_gridU.nc, mydir/GLORYS_0002_gridU.nc ..., mydir/GLORYS_2619_gridU.nc. Here indsize=4)
        - var_name : str, name of the considered variable in NetCDF file
        - var_lon : str, name of variable longitude in NetCDF file
        - var_lat : str, name of variable latitude in NetCDF file
        - mask : 'mask' if mask is same as meshfile, 'FillValue' if masked points are determined by a FillValue
        - mask_mat : if mask='mask', matrix containing continental mask (0 for continents, 1 for oceans)
        - ind : index of the file to be read
        """
        
        self.dir = dir
        self.prefix = prefix
        self.suffix = suffix
        self.indsize= indsize
        self.var_name = var_name
        self.var_lon = var_lon
        self.var_lat = var_lat
        self.mask = mask
        self.mask_mat = mask_mat

        #index of the file
        ind_str = str("%04d" %ind)
        filename = self.dir+self.prefix+ind_str+self.suffix
            
        if os.path.exists(filename):
            infile = nc.Dataset(filename)
            xdim = infile.dimensions[self.var_lon]
            ydim = infile.dimensions[self.var_lat]
            values = np.squeeze(infile.variables[self.var_name])
            self.values = np.float64(values,order='Fortran')
            lon = np.squeeze(infile.variables[self.var_lon])
            self.lon = np.float64(lon,order='Fortran')
            lat = np.squeeze(infile.variables[self.var_lat])
            self.lat = np.float64(lat,order='Fortran')
            if self.mask=='mask':
                #Mask point using known mask
                self.values = np.where(self.mask_mat==1,self.values,np.float64(0))
            elif self.mask.lower()=='fillvalue':
                #Mask points at FillValue
                fill_val = np.float64(infile.variables[self.var_name]._FillValue)
                self.values = np.where(self.values==fill_val,np.float64(0),self.values)
            infile.close()
            
        else:
            print("ERROR - %s does not exist" %filename)
            sys.exit(1)

    
    def UpdateDay(self,ind):
        """
        Update variable values for file of index ind
        """
        #index of the file
        ind_str = str("%04d" %ind)
        filename = self.dir+self.prefix+ind_str+self.suffix
            
        if os.path.exists(filename):
            infile = nc.Dataset(filename)
            xdim = infile.dimensions[self.var_lon]
            ydim = infile.dimensions[self.var_lat]
            values = np.squeeze(infile.variables[self.var_name])
            self.values = np.float64(values,order='Fortran')

            if self.mask=='mask':
                #Mask point using known mask
                self.values = np.where(self.mask_mat==1,self.values,np.float64(0))
            elif self.mask=='FillValue':
                #Mask points at FillValue
                fill_val = np.float64(infile.variables[self.var_name]._FillValue)
                self.values = np.where(self.values==fill_val,np.float64(0),self.values)
            infile.close()
            
        else:
            print("ERROR - %s does not exist" %filename)
            sys.exit(1)


###########
# Outputs #
###########

def add_nc_var(ncfile,varname,dtype,dim,values,attr):
    """
    Create a new variable in ncfile
    
    ncfile : NetCDF file
    varname : str, name of the variable
    dtype : str, type of the data
    dim : tuple, dimensions
    values : array, values of the data
    atr : dictionnary, list of attributes and values
    """

    nc_var = ncfile.createVariable(varname,dtype,dim)
    nc_var[:] = values[:] 
    #nc_var.assignValue(values)
    
    for key in attr:
        setattr(nc_var,key,attr[key])
    


def FillNetcdfFile(param,filename,init_x,init_y,init_t,final_x,final_y,final_t,traj_lon,traj_lat,traj_time,date,\
    traj_i0,traj_j0,u_current,v_current,traj_temp=0,traj_pp=0,habT=0,habPP=0,hab=0,xgrad=0,ygrad=0,u_swim=0,v_swim=0,u_tot=0,v_tot=0):

    """
    Fill a NetCDF file called filename with output data :
        - param : dictionnary of all parameters of the simulation
        - filename : str, name of the output file
        - initial positions (grid indices) and times : init_x, init_y, init_t
        - final positions (grid indices) and times : final_x,final_y,final_t
        - positions (longitude, latitude) and age at all time : traj_lon,traj_lat,traj_time
        - positions (grid indices) at all time : traj_i0,traj_j0
        - tracers along trajectories at all time (temperature, primary production...) : traj_temp,traj_pp
        - if 'active' mode, habitat along trajectories at all time : habT, habPP
    """
    
        #Create NetCDF file and dimensions
    out_nc=nc.Dataset(filename, 'w', format='NETCDF4')     
    out_nc.createDimension('nsteps',traj_lon.shape[0])
    out_nc.createDimension('nturtles',traj_lat.shape[1])

    #Write global attributes
    setattr(out_nc, 'title', 'Sea turtle trajectories computed by IBM2D')
    setattr(out_nc, 'species', param['species'])
    setattr(out_nc, 'initial_positions_file', param['init_file'])
    setattr(out_nc, 'tstep', param['tstep'])
    setattr(out_nc, 'nsteps_simu', param['nsteps_simu'])
    setattr(out_nc, 'mode', str(param['mode']))
    if param['mode']=='active':
        print(param['alpha'],param['vscale'],param['Fa'])
        setattr(out_nc, 'alpha', param['alpha'])
        setattr(out_nc, 'vscale', param['vscale'])
        setattr(out_nc, 'Fa', param['Fa'])
    setattr(out_nc, 'meshfile', param['dir_mesh']+param['fn_mesh'])

        #Write namelist items
#        for key in param:
#            setattr(out_nc,key,param[key])

    
        #Save variables
    print("--> Saving initial positions")
    add_nc_var(out_nc,'init_x','d',('nturtles',),init_x,
               {'standard_name':'init_x',
                'long_name':'Initial x position in grid',
                'valid_min':np.min(init_x),
                'valid_max':np.max(init_x)})

    add_nc_var(out_nc,'init_y','d',('nturtles',),init_y,
               {'standard_name':'init_y',
                'long_name':'Initial y position in grid',
                'valid_min':np.min(init_y),
                'valid_max':np.max(init_y)})

    add_nc_var(out_nc,'init_t','d',('nturtles',),init_t,
               {'standard_name':'init_t',
                'long_name':'Initial time',
                'valid_min':np.min(init_t),
                'valid_max':np.max(init_t),
                'units':'index of inner calendar'})
        

    print("--> Saving final positions")
    add_nc_var(out_nc,'final_x','d',('nturtles',),final_x,
               {'standard_name':'final_x',
                'long_name':'Final x position in grid',
                'valid_min':np.min(final_x),
                'valid_max':np.max(final_x)})

    add_nc_var(out_nc,'final_y','d',('nturtles',),final_y,
               {'standard_name':'final_y',
                'long_name':'Final y position in grid',
                'valid_min':np.min(final_y),
                'valid_max':np.max(final_y)})

    add_nc_var(out_nc,'final_t','d',('nturtles',),final_t,
               {'standard_name':'final_t',
                'long_name':'Final time',
                'valid_min':np.min(final_t),
                'valid_max':np.max(final_t),
                'unit':'index of inner calendar'})

    print("--> Saving longitudes")
    add_nc_var(out_nc,'traj_lon','d',('nsteps','nturtles'),traj_lon,
               {'standard_name':'traj_lon',
                'long_name':'Longitude of all trajectories at all time',
                'valid_min':np.min(traj_lon),
                'valid_max':np.max(traj_lon)})

    print("--> Saving latitudes")
    add_nc_var(out_nc,'traj_lat','d',('nsteps','nturtles'),traj_lat,
               {'standard_name':'traj_lat',
                'long_name':'Latitude of all trajectories at all time',
                'valid_min':np.min(traj_lat),
                'valid_max':np.max(traj_lat)})

    print("--> Saving x positions")
    add_nc_var(out_nc,'traj_i0','i',('nsteps','nturtles'),traj_i0,
               {'standard_name':'traj_i0',
                'long_name':'X position of ambient cell of all trajectories at all time',
                'valid_min':np.min(traj_i0),
                'valid_max':np.max(traj_i0)})

    print("--> Saving y positions")
    add_nc_var(out_nc,'traj_j0','i',('nsteps','nturtles'),traj_j0,
               {'standard_name':'traj_j0',
                'long_name':'Y position of ambient cell of all trajectories at all time',
                'valid_min':np.min(traj_j0),
                'valid_max':np.max(traj_j0)})

    print("--> Saving age of turtles")
    add_nc_var(out_nc,'traj_time','d',('nsteps','nturtles'),traj_time,
               {'standard_name':'traj_time',
                'long_name':'Time in inner calendar of all trajectories at all time',
                'valid_min':np.min(traj_time),
                'valid_max':np.max(traj_time),
                'units':'days'})

    add_nc_var(out_nc,'date','d',('nsteps','nturtles'),date,
               {'standard_name':'date',
                'long_name':'Time of all trajectories at all time',
                'valid_min':np.min(date),
                'valid_max':np.max(date),
                'units':'days'})

    print("--> Saving zonal current speed")
    add_nc_var(out_nc,'u_current','d',('nsteps','nturtles'),u_current,
               {'standard_name':'u_current',
                'long_name':'zonal current speed',
                'valid_min':np.min(u_current),
                'valid_max':np.max(u_current),
                'units':'meters per second'})
        
    print("--> Saving meridional current speed")
    add_nc_var(out_nc,'v_current','d',('nsteps','nturtles'),v_current,
               {'standard_name':'v_current',
                'long_name':'meridional current speed',
                'valid_min':np.min(v_current),
                'valid_max':np.max(v_current),
                'units':'meters per second'})  

    if param['key_alltracers']==True:
        print("--> Saving ambiant temperature")
        add_nc_var(out_nc,'traj_temp','d',('nsteps','nturtles'),traj_temp,
                   {'standard_name':'traj_temp',
                    'long_name':'Sea surface temperature along trajectories at all time',
                    'valid_min':np.min(traj_temp),
                    'valid_max':np.max(traj_temp),
                    'units':'Celsius degrees'})

        print("--> Saving ambiant food proxy")
        #Read name and unit in file
        ind_str = str("%04d" %param['ind0_pp'])
        ppfilename = param['c_dir_pp']+param['c_prefix_pp']+ind_str+param['c_suffix_pp']

        try:
            ppfile = nc.Dataset(ppfilename)
            ppvar = ppfile.variables[param['nc_var_pp']]
        except IOError:
            print("PP file doesn't exist")
            ppvar = 0

        try:
            ppname = ppvar.long_name
        except AttributeError:
            ppname = 'Food proxy'
        try:
            ppunit = ppvar.units
        except AttributeError:
            ppunit = '??'

        try:
            ppfile.close()
        except:
            pass
        
        add_nc_var(out_nc,'traj_pp','d',('nsteps','nturtles'),traj_pp,
                   {'standard_name':'traj_pp',
                    'long_name':ppname+' along trajectories at all time',
                    'valid_min':np.min(traj_pp),
                    'valid_max':np.max(traj_pp),
                    'units':ppunit})
            
    if param['mode']=='active':
        print("--> Saving ambiant temperature habitat")
        add_nc_var(out_nc,'habT','d',('nsteps','nturtles'),habT,
                   {'standard_name':'habT',
                    'long_name':'Temperature habitat along trajectories at all time',
                    'valid_min':np.min(habT),
                    'valid_max':np.max(habT)})

        print("--> Saving ambiant foraging habitat")
        add_nc_var(out_nc,'habPP','d',('nsteps','nturtles'),habPP,
                   {'standard_name':'habPP',
                    'long_name':'Food habitat along trajectories at all time',
                    'valid_min':np.min(habPP),
                    'valid_max':np.max(habPP)})

        print("--> Saving ambiant habitat")
        add_nc_var(out_nc,'hab','d',('nsteps','nturtles'),hab,
                   {'standard_name':'hab',
                    'long_name':'Total habitat along trajectories at all time',
                    'valid_min':np.min(hab),
                    'valid_max':np.max(hab)})

        print("--> Saving zonal habitat gradient")
        add_nc_var(out_nc,'xgrad','d',('nsteps','nturtles'),xgrad,
                   {'standard_name':'xgrad',
                    'long_name':'Zonal habitat gradient along trajectories at all time',
                    'valid_min':np.min(xgrad),
                    'valid_max':np.max(xgrad)})

        print("--> Saving meridionnal habitat gradient")
        add_nc_var(out_nc,'ygrad','d',('nsteps','nturtles'),ygrad,
                   {'standard_name':'ygrad',
                    'long_name':'Meridionnal habitat gradient along trajectories at all time',
                    'valid_min':np.min(ygrad),
                    'valid_max':np.max(ygrad)})

        print("--> Saving zonal swimming speed")
        add_nc_var(out_nc,'u_swim','d',('nsteps','nturtles'),u_swim,
                   {'standard_name':'u_swim',
                    'long_name':'zonal swimming speed due to habitat',
                    'valid_min':np.min(u_swim),
                    'valid_max':np.max(u_swim),
                    'units':'meters per second'})

        print("--> Saving meridional swimming speed")
        add_nc_var(out_nc,'v_swim','d',('nsteps','nturtles'),v_swim,
                   {'standard_name':'v_swim',
                    'long_name':'meridional swimming speed due to habitat',
                    'valid_min':np.min(v_swim),
                    'valid_max':np.max(v_swim),
                    'units':'meters per second'})
                    

        print("--> Saving zonal total speed")
        add_nc_var(out_nc,'u_tot','d',('nsteps','nturtles'),u_tot,
                   {'standard_name':'u_tot',
                    'long_name':'zonal total speed',
                    'valid_min':np.min(u_tot),
                    'valid_max':np.max(u_tot),
                    'units':'meters per second'})
                 
        print("--> Saving meridional total speed")
        add_nc_var(out_nc,'v_tot','d',('nsteps','nturtles'),v_tot,
                   {'standard_name':'v_tot',
                    'long_name':'meridional total speed',
                    'valid_min':np.min(v_tot),
                    'valid_max':np.max(v_tot),
                    'units':'meters per second'})

            #~ print "--> Saving zonal habitat gradient"
        #~ add_nc_var(out_nc,'xgradh','d',('nsteps','nturtles'),xgradh,
                   #~ {'standard_name':'xgradh',
                    #~ 'long_name':'zonal habitat gradient component',
                    #~ 'valid_min':np.min(xgradh),
                    #~ 'valid_max':np.max(xgradh),
                        #~ 'units':'hab per meter'})
                    #~ 
        #~ print("--> Saving meridional habitat gradient")
        #~ add_nc_var(out_nc,'ygradh','d',('nsteps','nturtles'),ygradh,
                   #~ {'standard_name':'ygradh',
                    #~ 'long_name':'meridional habitat gradient component',
                    #~ 'valid_min':np.min(ygradh),
                    #~ 'valid_max':np.max(ygradh),
                    #~ 'units':'hab per meter'})

                     

    
    out_nc.close()

