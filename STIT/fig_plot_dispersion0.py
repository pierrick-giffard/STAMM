#!/usr/bin/env python
#-*- coding:utf-8 -*-
import os, sys
import numpy as np
import netCDF4 as nc
#from datetime import date
import matplotlib.pyplot as plt
#import matplotlib.pylab as pl
#from matplotlib.patches import Rectangle
from matplotlib import cm
from matplotlib import gridspec
import random
from mpl_toolkits.basemap import Basemap#, shiftgrid
#from mpl_toolkits.axes_grid1 import make_axes_locatable
#from pylab import rcParams

class data:
    """
    class data contains all results from simulation
    contains also all routines for data analysis
    """

    def __init__(self,filename=str()):
        """
        called with the path to the netcdf output file as an argument
        load all data

        filename : string, path to netcdf file
        """

        # file opening
        if os.path.exists(filename):
            infile = nc.Dataset(filename)
        else :
            sys.exit(str(filename)+"\n\n\n file does not exist\n\n")

        # loading data
        try:
            self.lon=infile.variables['traj_lon'][:,:]
        except:
            self.lon=np.transpose(infile.variables['lon'][:,:])
        self.lon=np.float32(self.lon)
        print('   => longitude loaded',self.lon.shape)
        #
        try:
            self.lat=infile.variables['traj_lat'][:,:]
        except:
            self.lat=np.transpose(infile.variables['lat'][:,:])
        self.lat=np.float32(self.lat)
        print('   => latitude loaded',self.lat.shape)
        #
        try:
            self.traj_time = infile.variables['traj_time'][:,:]
        except:
            self.traj_time = np.transpose(infile.variables['age'][:,:])
        self.traj_time = np.int32(self.traj_time)
        print('   => traj_time loaded')
        #
        try:
            self.init_t=np.float32(infile.variables['init_t'][:])
        except:
            self.init_t=np.float32(infile.variables['time'][:,0])/86400
        #Shape of data
        self.nb_output,self.nturtles=self.lon.shape
        infile.close()

        # assigning data
        self.filename=filename
        self.nb_output,self.nturtles=self.lat.shape
        self.nb_year=(self.nb_output-1)/365. # for daily outputs
        self.lat0=np.mean(self.lat[0,:])
        self.lon0=np.mean(self.lon[0,:])
        print('Zone de depart : ' + str(self.lat0) + ', ' + str(self.lon0))

        #Center longitudes
        self.lon = np.where(self.lon<=self.lon0-180,self.lon+360,self.lon)
        self.lon = np.where(self.lon>self.lon0+180,self.lon-360,self.lon)

        #Ariane days
        days = np.zeros(self.lon.shape,dtype='int32')
        for k in range(self.nturtles):
            days[:,k] = np.int32((self.traj_time[:,k]-self.traj_time[0,k])+self.init_t[k])
        self.days=days

def getTicks(lmin, lmax, step):
    lmaxabs = (int(max(abs(lmax), abs(lmin)))/step+1)*step
    return np.intersect1d(np.arange(lmin+1, lmax), np.arange(-lmaxabs, lmaxabs, step))


def plot_map(ars,xmin,xmax,ymin,ymax,value=0.6,res=0.25,alpha=1,lat_space=10,lon_space=20):
    #Basemap reference
    map=Basemap(ax=ars,llcrnrlon=xmin,llcrnrlat=ymin,urcrnrlon=xmax,urcrnrlat=ymax,projection='cyl',resolution='l')
    map.fillcontinents(color='0.35')
    map.drawcoastlines(color='grey',linewidth=0.2)
    map.drawcountries(color='k',linewidth=0.01)
    map.drawparallels(getTicks(ymin, ymax, lat_space), labels=[1,0,0,0], fontsize=8, linewidth=0.2)
    map.drawmeridians(getTicks(xmin, xmax, lon_space), labels=[0,0,0,1], fontsize=8, linewidth=0.2)



def display_trajectories(dataFile,f,ax,lw=0.005, ms=0.04, col='b', alpha=0.5,show_start=True):

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
    y = np.array(dataFile.lat)

    if np.mean(dataFile.lon[0,:]) < 36 :
        #Si point de départ dans l'Atlantique, on corrige l'affichage des longitudes > 180°E(cad toutes puisque 180°E est dans le Pacifique et que l'on se trouve dans l'Atlantique)
        x = np.array([l+(((1-np.sign(l))/2)*360) for l in x]) 
        # Lorsque les particules dépassent greenwich on ajoute 360 (elles passent de 359 à 361 plutot que de 359 à 1, par exemple en mediterranée)
        x[np.where(x<200)]=x[np.where(x<200)]+360
        

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
        p=ax.scatter(xbis, ybis, c=col, s=ms, edgecolor='none', vmin=colmin, vmax=colmax, alpha=alpha)


    #Shows starting point
    dataFile.lon = np.array([l+(((1-np.sign(l))/2)*360) for l in dataFile.lon])
    ax.plot((np.mean(dataFile.lon[0,:]),),(np.mean(dataFile.lat[0,:]),),markerfacecolor='w',markeredgecolor='k',marker='o',ms=6,mew=0.3,zorder=999)
    
    return p,time

def display_colorbar(f,im, ax_cb, label):
    cb=f.colorbar(im,cax=ax_cb,orientation='horizontal')
    cb.solids.set(alpha=1)
    cb.ax.tick_params(labelsize=8)
    cb.set_label(label, labelpad=1, size=8)
    cb.outline.set_linewidth(0.5)
    cb.ax.xaxis.set_tick_params(width=0.5)



def main():

    ifile = sys.argv[1]
    ofile = ifile.replace(".nc",".png")
    
    # input_filename = 'Yalimapo_Cayenne_5000indiv_NATLpaper_2002_18y_PPmax60_alpha3e6_vscale1.2_speedrecord_yesrebond'
    
    # output_path='/data/MEMMS/IBM/visu/'
    # output_filename = 'test'
    
    ## -- Si l'affichage coupe les trace à 0° de longitude, bien vérifier que x n'est pas corrigé dans 
    ## -- PlotTrajectories2(display_trajectories_all) : #x = np.array([l+(((1-np.sign(l))/2)*360) for l in x])
    ## -- il faut que cette ligne soit bien commentée

    #Entry parameters --------

    # Definition de la zone d'affichage
    xmin = 255
    xmax = 396
    ymin = -5
    ymax = 62

    # xmin = 255
    # xmax = 360
    # ymin = -5
    # ymax = 30
   

    #Afficher une colorbar/point de depart
    colorbar = True
    start = True

    # Plot figure ------
    c = 0.89
    f = plt.figure(figsize = (12*c/2.54,8*c/2.54))
    gs = gridspec.GridSpec(2,1,height_ratios=[11,1],left=0.08, right=0.98, bottom=0.07, top=0.95)

    dataFile = data(ifile)
    ax = plt.subplot(gs[0])
    
    im,time = display_trajectories(dataFile,f,ax)
    plot_map(ax,xmin,xmax,ymin,ymax)

    ax.spines['right'].set_linewidth(0.5)
    ax.spines['left'].set_linewidth(0.5)
    ax.spines['bottom'].set_linewidth(0.5)
    ax.spines['top'].set_linewidth(0.5)
        

    if np.nanmax(time)>365: 
        ax_cb = plt.subplot(gs[1])
        label = u"Age (Years)"
        display_colorbar(f,im, ax_cb, label)

    else:
        ax_cb = plt.subplot(gs[2])
        label = u"Age (Days)"
        display_colorbar(f,im, ax_cb, label)
            
        
    plt.savefig(ofile,bbox_inches='tight',dpi=800)
    print("wrote {}".format(ofile))
    plt.close()


if __name__ == '__main__':
    main()
    
