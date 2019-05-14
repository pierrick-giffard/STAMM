#!/usr/bin/env python
#-*- coding:utf-8 -*-

"""
	Author : Amaury Dehecq (intern at CLS 2011-2012 working with Phillipe Gaspar)

	Module description : contain all useful information on a turtle population and all possible operations on this population (size & mass computation, habitat computation, etc)

        Last modified : 10/08/2012
"""

#Python libraries
import numpy as np
from math import *
import random
from scipy import misc 
#Personal libraries
from mod_fx_fy import fx_fy as sf



class Turtle():
    
    def __init__(self,x_init,y_init,t_init,nsteps, key_alltracers, mode):
        """
        Contain all useful information on a population of leatherback turtles.

        x_init, y_init : array, initial position of all turtles in grid indices
        t_init : array, initial time of all turtles
        nsteps : f, number of simulation steps
        key_alltracers : bool, are optional tracers to be saved
        mode : str, 'passive', 'active' or 'diffusion'
        """

        nturtles = len(x_init)
        self.nturtles = nturtles
        self.nsteps = nsteps
        self.key_alltracers = key_alltracers
        #Current position in grid indices (x -> U grid, y -> V grid)
        self.x = np.zeros(nturtles,dtype='float64')
        self.y = np.zeros(nturtles,dtype='float64')
        #Current ambient cell (T grid)
        self.i0 = np.zeros(nturtles,dtype='int32')
        self.j0 = np.zeros(nturtles,dtype='int32')
        #Contain indices of ambient cell at any time
        self.traji0 = np.zeros((nsteps+1,nturtles),dtype='int32')
        self.trajj0 = np.zeros((nsteps+1,nturtles),dtype='int32')
        #Contain all positions at any time in lon/lat coordinates
        self.trajlon = np.zeros((nsteps+1,nturtles),dtype='float64')
        self.trajlat = np.zeros((nsteps+1,nturtles),dtype='float64')
        #Initial time of all turtles
        self.t_init = np.zeros(nturtles,dtype='float64')
        #Age (in days) of the turtle at any time and current time indice
        self.age = np.zeros((nsteps+1,nturtles),dtype='float64')
        self.tind = np.zeros((nturtles),dtype='int32')

        self.age_now = np.zeros((nturtles),dtype='float64')
        #self.age_now = np.ones((nturtles),dtype='float64')*6570
        self.u_current = np.zeros((nsteps+1,nturtles),dtype='float64')
        self.v_current = np.zeros((nsteps+1,nturtles),dtype='float64')
                
        if key_alltracers==True:
            #Temperature of the environment at any time and current time
            #Temperature is evaluated at the centre of the cell not interpolated as in Ariane
            self.Tenv = np.zeros((nsteps+1,nturtles),dtype='float64')
            self.Tenv_now = np.zeros((nturtles),dtype='float64')
            #Primary production of the environment at any time and current time
            #PP is evaluated at the centre of the cell not interpolated as in Ariane
            self.PP = np.zeros((nsteps+1,nturtles),dtype='float64')
            self.PP_now = np.zeros((nturtles),dtype='float64')
        else:
            self.Tenv = 0
            self.PP = 0
        
        if mode=='active':
            #Temperature and food habitat at any time
            self.habT = np.zeros((nsteps+1,nturtles),dtype='float64')
            self.habPP = np.zeros((nsteps+1,nturtles),dtype='float64')
            self.hab = np.zeros((nsteps+1,nturtles),dtype='float64')
            self.xgrad = np.zeros((nsteps+1,nturtles),dtype='float64')
            self.ygrad = np.zeros((nsteps+1,nturtles),dtype='float64')
            self.u_swim = np.zeros((nsteps+1,nturtles),dtype='float64')
            self.v_swim = np.zeros((nsteps+1,nturtles),dtype='float64')
            self.u_hab = np.zeros((nsteps+1,nturtles),dtype='float64')
            self.v_hab = np.zeros((nsteps+1,nturtles),dtype='float64')
            self.u_rand = np.zeros((nsteps+1,nturtles),dtype='float64')
            self.v_rand = np.zeros((nsteps+1,nturtles),dtype='float64')
            self.u_tot = np.zeros((nsteps+1,nturtles),dtype='float64')
            self.v_tot = np.zeros((nsteps+1,nturtles),dtype='float64')

        else:
            self.habT = 0
            self.habPP = 0
            self.hab = 0
            self.xgrad = 0
            self.ygrad = 0
            self.u_swim = 0
            self.v_swim = 0
            self.u_hab = 0
            self.v_hab = 0
            self.u_rand = 0
            self.v_rand = 0
            self.u_tot = 0
            self.v_tot = 0


        #State of the turtle, 0 for inactive, 1 for active, 2 for out of domain
        self.state = np.ones(nturtles,dtype='int32') 

        #Initialisation
        self.x[:] = x_init
        self.y[:] = y_init
        self.t_init = t_init
        self.traji0[0,:] = np.int32(self.x[:])+1
        self.trajj0[0,:] = np.int32(self.y[:])+1
        self.i0[:] = np.int32(self.x[:])+1
        self.j0[:] = np.int32(self.y[:])+1



    def UpDatePosition(self,grid):
        """
        Update position at current time indice
        
        grid : Object of type IOlib.Grid containing glamu, gphiv, conversion matrices from indices to lon/lat and east/west periodicity.
        Call Fortran 90/95 subroutines fx & fy 
         """
        active_turtles = np.where(self.state==1)[0]
        #print(active_turtles)
        if len(active_turtles)>0:
            x = self.x[active_turtles]+1
            y = self.y[active_turtles]+1
            lon, lat = np.zeros(self.nturtles,dtype='float64'), np.zeros(self.nturtles,dtype='float64')
            lon[active_turtles], lat[active_turtles] = sf.fx(x,y,grid.glamu,grid.key_periodic,grid.overlap), sf.fy(x,y,grid.gphiv,grid.key_periodic)
            #print(self.tind           )
            for j in active_turtles:
                i=self.tind[j]
                self.traji0[i,j] = np.int32(self.i0[j])
                self.trajj0[i,j] = np.int32(self.j0[j])
                self.trajlon[i,j] = lon[j]
                #print('lon', lon[j])
                self.trajlat[i,j] = lat[j]
                #print('lat', lat[j])
                #print('distance lon', ((self.trajlon[i,j]*np.pi/180.)-(self.trajlon[i-1,j]*np.pi/180.))*6371.*np.cos(((self.trajlat[i,j]+self.trajlat[i-1,j])/2.)*np.pi/180.))



    def UpDateTemperature(self,temp):
        """
        Update temperature at current time indice.
        
        temp : object of class IOlib.Variable
        """
        active_turtles = np.where(self.state==1)[0]
        for j in active_turtles:
            i=self.tind[j]
            self.Tenv_now[j] = temp.values[self.j0[j],self.i0[j]]
            self.Tenv[i,j] = self.Tenv_now[j]




    def UpDatePP(self,PP):
        """
        Update primary production at current time indice.
        
        PP : object of class IOlib.Variable
        """
        active_turtles = np.where(self.state==1)[0]
        for j in active_turtles:
            i=self.tind[j]
            self.PP_now[j] = PP.values[self.j0[j],self.i0[j]]
            self.PP[i,j] = self.PP_now[j]


class Leatherback(Turtle):

    def __init__(self,x_init,y_init,t_init,nsteps, key_alltracers, mode):
        
        Turtle.__init__(self,x_init,y_init,t_init,nsteps, key_alltracers, mode)

        
    def compute_SCL(self):
        """
        Compute Straight Carapace Length (m) of all turtles

        Ref : Jones, T.T., Hastings, M.D., Bostrom, B.L., Pauly, D., Jones, D.R., 2011. Growth of captive leatherback turtles, Dermochelys coriacea, with inferences on growth in the wild: Implications for population decline and recovery. Journal of Experimental Marine Biology and Ecology 399, 84–92.
        """
        
        self.SCL = 1.43*(1-np.exp(-0.226*(self.age_now/365.+0.17)))
            
    def compute_M(self):
        """
        Compute all turtles mass(kg)

        Ref : Jones, T.T., Hastings, M.D., Bostrom, B.L., Pauly, D., Jones, D.R., 2011. Growth of captive leatherback turtles, Dermochelys coriacea, with inferences on growth in the wild: Implications for population decline and recovery. Journal of Experimental Marine Biology and Ecology 399, 84–92.
        """
        self.M = 112.31*(self.SCL)**2.86
    
    def compute_PPmax(self,Fa):
        """
        Compute food threshold

        Ref : Jones, T.T., Bostrom, B.L.,  Hastings, M.D., Van Houtan K.S., Pauly, D., 2012. Resource requirements of the Pacific Leatherback Turtle Population
        """
        #Besoin annuels en valeur absolue

        #self.PPmax = 312*2.86*0.299*(((1-np.exp(-0.299*(self.age_now/365.+0.17)))**(2.86-1))*(np.exp(-0.299*(self.age_now/365.+0.17))))/(1-(1-np.exp(-0.299*(self.age_now/365.+0.17)))**(2.86*0.0328))
        self.PPmax = 266.80368*(((1-np.exp(-0.299*(self.age_now/365.+0.17)))**(2.86-1))*(np.exp(-0.299*(self.age_now/365.+0.17))))/(1-(1-np.exp(-0.299*(self.age_now/365.+0.17)))**(2.86*0.0328))
        self.PPmax = self.PPmax/(2835.24/float(Fa))
        

    def compute_vmax(self,vscale):
        
        """
        Compute maximum speed for each turtle

        Ref : Gaspar, P., Benson, S., Dutton, P., Réveillère, A., Jacob, G., Meetoo, C., Dehecq, A., Fossette, S., 2012. Oceanic dispersal of juvenile leatherback turtles: going beyond passive drift modeling. Marine Ecology Progress Series.
        """

        # Version originale avec la distance moyenne parcourue en une journée :
        # Pour une tortue de 1.43m, on a 50 km (35000*1.43) que l'on divise par
        # le nombre de secondes dans une journée (86400)
        # vmax = 3.5e4*self.SCL/86400

        # Version avec vitesse optimale = SCL**0.126
        #self.vmax = 1.5*(self.SCL**0.126)
        #print float(vscale)
        self.vmax = float(vscale)*(self.SCL**0.126)
        #print self.vmax
        #print self.SCL
        
    def habitat(self,turtle,Temp,PP,PPmax,Total=True):
        """
        Compute habitat for leatherback turtles

        Temperature habitat is a gaussian centered about an optimum temperature mu with standard deviation sigma 
        Both parameters still need to be evaluated. I used temperature gradient as estimated by Bostrom to get optimal body temperature. mu-sigma is minimum temperature used in Gaspar et al., 2012

        Food habitat is a linear function of PP value. 0 for 0, 1 for max (decile 95) of PP
        """



        #Temperature habitat
        
        #SCL = self.SCL[turtle]
        Mass = self.M[turtle]
        
        #Tbo = 24.4	#Optimal temperature (nesting beach) from Gaspar et al 2012
        ##########Tbo = 24.4-5.9*SCL	#Test pour diminuer temperature optimale avec l'age	
        #Tmin = Tbo-7.9*SCL     #Critical minimum temperature from Gaspar et al 2012
        #mu = Tbo
        ##########
        #Gaussienne
        ##########
        #sigma = (Tbo-Tmin)/2
        #T_hab=np.exp(-(Temp-Tbo)**2/(2*sigma**2))
        
        ##########
        #T_logistique
        ##########
        #sigma = (Tbo-Tmin)/7
        #T_hab = 1/(1+np.exp(-(Temp-Tmin)/sigma))
         
        ##########
        #Gaussienne tronquée
        ##########
        #sigma = (Tbo-Tmin)/2
        #if Temp <= Tbo:
            #T_hab =  np.exp(-(Temp-Tbo)**2/(2*sigma**2))
        #else :    
            #T_hab = 1.
        ##########
        #Sinus, Hab_T asymatrique. Topt et Tmin variant LINEAIREMENT avec la masse. 
        ##########
        #Topt_newborn = 24.4
        #Topt_grownup = 17.
        #Tmin_newborn = 23.
        #Tmin_grownup = 5.
        #Topt = ((Topt_grownup - Topt_newborn)/(312.))*Mass + Topt_newborn
        #Tmin = ((Tmin_grownup - Tmin_newborn)/(312.))*Mass + Tmin_newborn
        #if Temp > Topt:
            #T_hab = 1.0
        #elif Temp < Tmin:
            #T_hab = 0.0
        #else:
            #T_hab = 0.5 * (1-np.sin((np.pi*(Tmin+Topt-2*Temp))/(-2*(Tmin-Topt))))
        
        ###########
        #THERMIQUE PHYSIQUE - P.Gaspar
        #Sinus, Hab_T asymetrique. Seuils Topt et Tmin variant en fonction de la masse
        ###########
        #Topt = 24.4-(0.464*np.sqrt(Mass))
        #Tmin = Topt - (0.35*((Mass)**(2./3.)))
        #if Temp > Topt:
            #T_hab = 1.0
        #elif Temp < Tmin:
            #T_hab = 0.0
        #else:
            #T_hab = 0.5 * (1-np.sin((np.pi*(Tmin+Topt-2*Temp))/(-2*(Tmin-Topt))))
    
        ###########
        #THERMIQUE PRODUCTION DE CHALEUR - Bostrom & Jones
        #Sinus, Hab_T asymetrique.
        ###########
        #Topt = 24.4 - 0.21*np.sqrt(Mass)
        #Tmin = 24.4 - 0.84*np.sqrt(Mass)
        #if Temp > Topt:
            #T_hab = 1.0
        #elif Temp < Tmin:
            #T_hab = 0.0
        #else:
            #T_hab = 0.5 * (1-np.sin((np.pi*(Tmin+Topt-2*Temp))/(-2*(Tmin-Topt))))
        
        ###########
        #THERMIQUE PRODUCTION DE CHALEUR - Bostrom & Jones
        #Gaussienne, Hab_T asymetrique.
        ###########
        
        Topt = 24. - 0.21*np.sqrt(Mass)
        Tmin = 24. - 1.05*np.sqrt(Mass)
        #sigma = (Topt-Tmin)/2.
        if Temp >= Topt:
            T_hab = 1.0
        else:
            #T_hab = np.exp(-(Temp-Topt)**2/(2*sigma**2))
             T_hab = np.exp(-2*((Temp-Topt)/(Topt-Tmin))**2)
        
        ###########
        #THERMIQUE PRODUCTION DE CHALEUR - Bostrom & Jones
        #Gaussienne, Hab_T SYMETRIQUE.
        ###########
        """
        Topt = 24.4 - 0.21*np.sqrt(Mass)
        Tmin = 24.4  - 1.05*np.sqrt(Mass)
        #Tmin = 26. - 0.84*np.sqrt(Mass)
        sigma = (Topt-Tmin)/2

        T_hab = np.exp(-(Temp-Topt)**2/(2*sigma**2))
        
        """


        ###########
        #THERMIQUE PRODUCTION DE CHALEUR - Bostrom & Jones
        #Gaussienne, Hab_T asymetrique.
        ###########
        """
        Topt_min = 24.4 - 0.21*np.sqrt(Mass)
        Topt_max = 29
        Tmin = 24.4 - 1.05*np.sqrt(Mass)
        Tmax =32

        sigma_inf = (Topt_min-Tmin)/2
        sigma_sup = (Tmax-Topt_max)/2

        if Temp > Topt_min and Temp< Topt_max:
            T_hab = 1.0

        elif Temp <= Topt_min :
            T_hab = np.exp(-(Temp-Topt_min)**2/(2*sigma_inf**2))

        elif Temp >= Topt_max :
            T_hab = np.exp(-(Temp-Topt_max)**2/(2*sigma_sup**2))
        """
        #Food habitat
        food_hab = min(PP/PPmax,1)
        #food_hab = min(PP/56,1)

        if Total==True:
            return T_hab*food_hab
        else:
            return T_hab,food_hab


        #Complete habitat
        """
        #Temperature habitat
        SCL = self.SCL[turtle]
        mass = 0.000224*(SCL*100)**2.78  # T.T Jones 2011
        Tbo = 29       #Optimal temperature from Southwood et al 2005
        Tbc = 24.4-7.9*SCL     #Critical minimum temperature from Gaspar et al 2012
        
        deltaT = log(1./0.85)*0.22*mass**0.83/(2*pi*0.25*8*0.036*mass**(1./3))  #Bostrom 2006
        
        mu = Tbo - deltaT
        sigma = (Tbo-Tbc)/2
        
        T_hab =  np.exp(-(Temp-mu)**2/(2*sigma**2))
        
        #Food habitat
        food_hab = min(PP/PPmax,1)
        
        if Total==True:
            return T_hab*food_hab
        else:
            return T_hab, food_hab
        """

class Green(Turtle):

    def __init__(self,x_init,y_init,t_init,nsteps, key_alltracers, mode):
        
        Turtle.__init__(self,x_init,y_init,t_init,nsteps, key_alltracers, mode)

        
    def compute_SCL(self):
        """
        Compute Straight Carapace Length of all turtles

        Ref : Bjorndal K.A., Bolten A.B., Atilio L.C. Jr., Kleiber P. (1995) Estimation of green turtle (Chelonia mydas) growth rates from length-frequency analysis. Copeia, 1, 71-77.
        """

        self.SCL = 0.983*(1-np.exp(-0.949*(self.age_now/365.+0.074)))

    
    def compute_vmax(self):
        """
        Compute maximum speed for each turtle

        Ref : Gaspar, P., Benson, S., Dutton, P., Réveillère, A., Jacob, G., Meetoo, C., Dehecq, A., Fossette, S., 2012. Oceanic dispersal of juvenile leatherback turtles: going beyond passive drift modeling. Marine Ecology Progress Series.
        """

        self.vmax = 3.5e4*self.SCL/86400


    def habitat(self,turtle,Temp,PP,PPmax,Total=True):
        """
        Compute habitat for green turtles

        Temperature habitat is a gaussian centered about an optimum temperature mu with standard deviation sigma 
        Both parameters still need to be implemented.

        Food habitat is a linear function of PP value. 0 for 0, 1 for max (decile 95) of PP
        """

        #Temperature habitat
        Tbo = 18.0
        Tbc = 18.0-((1/0.983)*(18.0-15.0))*self.SCL
        deltaT = 0 #??
        mu = Tbo - deltaT
        sigma = (Tbo-Tbc)/2
        T_hab =  np.exp(-(Temp-mu)**2/(2*sigma**2))
        
        #Food habitat
        food_hab = min(PP/PPmax,1)

        if Total==True:
	        return T_hab
            #return T_hab*food_hab
        else:
            return T_hab, food_hab

class Loggerhead(Turtle):

    def __init__(self,x_init,y_init,t_init,nsteps, key_alltracers, mode):
        
        Turtle.__init__(self,x_init,y_init,t_init,nsteps, key_alltracers, mode)

        
    def compute_SCL(self):

        """
        Compute Straight Carapace Length (m) of all turtles

        Ref : Jones, T.T., Hastings, M.D., Bostrom, B.L., Pauly, D., Jones, D.R., 2011. Growth of captive leatherback turtles, Dermochelys coriacea, with inferences on growth in the wild: Implications for population decline and recovery. Journal of Experimental Marine Biology and Ecology 399, 84–92.
        """
        ##VB
        #A=1.088
        #B=0.9649
        #k=0.0739
        ## !! CCL
        A=1.088
        k=0.0981
        self.SCL=A*(1-np.exp(-k*((self.age_now/365.)+0.36)))



    def compute_M(self):

        """
        Compute all turtles mass(kg)
        """
        A= 0.00104
        b = 2.5
        self.M = A*(self.SCL*100)**b

    

    def compute_PPmax(self,Fa):
        """
        Compute food threshold

        Ref : Jones, T.T., Bostrom, B.L.,  Hastings, M.D., Van Houtan K.S., Pauly, D., 2012. Resource requirements of the Pacific Leatherback Turtle Population
        """
        #Besoin annuels en valeur absolue
        self.PPmax = 0.195*(((1-np.exp(-0.0981*((self.age_now)/365.+0.36)))**(1.5))*(np.exp(-0.0981*((self.age_now/365.)+0.36))))/(1-(1-np.exp(-0.0981*((self.age_now/365.)+0.36)))**(0.195))
        self.PPmax = self.PPmax*float(Fa)
        print(self.PPmax)

    def compute_vmax(self,vscale):
        
        """
        Compute maximum speed for each turtle
        Ref : Gaspar, P., Benson, S., Dutton, P., Réveillère, A., Jacob, G., Meetoo, C., Dehecq, A., Fossette, S., 2012. Oceanic dispersal of juvenile leatherback turtles: going beyond passive drift modeling. Marine Ecology Progress Series.
        """
        self.vmax = float(vscale)*(self.SCL**0.025)


    def habitat(self,turtle,Temp,PP,PPmax,Total=True):
        """
        Compute habitat for leatherback turtles

        Temperature habitat is a gaussian centered about an optimum temperature mu with standard deviation sigma 
        Both parameters still need to be evaluated. I used temperature gradient as estimated by Bostrom to get optimal body temperature. mu-sigma is minimum temperature used in Gaspar et al., 2012

        Food habitat is a linear function of PP value. 0 for 0, 1 for max (decile 95) of PP
        """

        ###########
        #THERMIQUE PRODUCTION DE CHALEUR - Bostrom & Jones
        #Gaussienne, Hab_T asymetrique.
        ###########
        #T1 = 14.
        #T2 = 18.
        #T3 = 29.
        #T4 = 32.

        T1 = 10.
        T2 = 18.
        T3 = 30.
        T4 = 33.
        if Temp >= T2 and Temp<=T3:
            T_hab = 1.0

        elif Temp < T2 :
            T_hab = np.exp(-2*((Temp-T2)/(T2-T1))**2)

        elif Temp > T3 :
            T_hab = np.exp(-2*((Temp-T3)/(T4-T3))**2)

        #Food habitat
        food_hab = min(PP/PPmax,1)
        #food_hab = min(PP/56,1)

        if Total==True:
            return T_hab*food_hab
        else:
            return T_hab,food_hab

