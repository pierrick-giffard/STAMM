#!/usr/bin/env python
#-*- coding:utf-8 -*-

"""
	Author : Amaury Dehecq (intern at CLS 2011-2012 working with Phillipe Gaspar)

	Module description : contain all useful information for transport of individuals. Computation of zonal and meridional transports, taxis and advection.
	Call a Fortran 90/95 subroutine advection2d
	
	Last modified : 10/08/2012
"""

#Python libraries
import numpy as np
from math import *
import random

#Personal libraries
from Advection2D_f import advection2d as adv

class Transport():
	
	def __init__(self,pop,uu,vv,grid,mode,Temp=0,PP=0):
		"""
		Objects of class Transport contain all needed information for the computation of trajectories : advection by ocean currents + taxis
		pop : object of class Turtle containing positions and caracteristics of all turtles
		uu : matrix of zonal current velocity at time of computation
		vv : matrix of meridional current velocity at time of computation
		grid : object of class Grid containing all parameters of the grid
		mode : str, 'passive', 'active' or 'diffusion'
		Temp : (optional) matrix of sea surface temperature at time of computation
		PP : (optional) matrix of primary production (or any food proxy) at time of computation	
		"""
		self.pop = pop
		self.uu = uu
		self.vv = vv
		self.grid = grid
		self.e2u = grid.e2u
		self.e1v = grid.e1v
		self.mask = grid.mask
		self.umask = grid.umask
		self.vmask = grid.vmask
		self.Temp = Temp
		self.PP = PP
		self.mode = mode


	def UpDate(self,pop,uu,vv,Temp=0,PP=0):
		"""
		Update values for new files
		"""
		self.pop = pop
		self.uu = uu
		self.vv = vv
		self.Temp = Temp
		self.PP = PP


	def compute_TotalTransport(self,alpha,traceur_pp):
		"""
		Defines functions to be used in the computation of trajectories zontr (zonal transport) and mertr (meridional transport)
		Ariane's algorithm works with transport, which is the product of the velocity times the cross section of the grid cell
		"""

		if self.mode=='passive':
			"""
			In the case of passive drift, the transport is composed of ocean transport only.
			"""
			active_turtles = np.where(self.pop.state==1)[0]
			vect_speed_ucurrent = []
			vect_speed_vcurrent = []

			def zontr(j,i,k,side,tj,sign,time_drifting,TotalTime,paire,saving):
				turtle = active_turtles[k-1]
				i0 = i-1     
				j0 = j-1
				#Fortran->Python array : -1
				if paire == 1:
				    vect_speed_ucurrent.append(self.uu[j0,i0])

				elif paire == 2 :
				    vect_speed_ucurrent[-1]=(vect_speed_ucurrent[-1]+self.uu[j0,i0])/2.

				if saving == 1:
				    self.pop.u_current[self.pop.tind[turtle],turtle] = self.pop.u_current[self.pop.tind[turtle],turtle]+((time_drifting/86400.)*vect_speed_ucurrent[-1])
                               		
				return self.uu[j0,i0]*self.e2u[j0,i0]



			def mertr(j,i,k,side,tj,sign,time_drifting,TotalTime,paire,saving):
				turtle = active_turtles[k-1]
				i0 = i-1     
				j0 = j-1
				#Fortran->Python array : -1
				if paire == 1:
				    vect_speed_vcurrent.append(self.vv[j0,i0])
				elif paire == 2 :
				    vect_speed_vcurrent[-1]=(vect_speed_vcurrent[-1]+self.vv[j0,i0])/2.
				if saving == 1:
				    self.pop.v_current[self.pop.tind[turtle],turtle] = self.pop.v_current[self.pop.tind[turtle],turtle]+((time_drifting/86400.)*vect_speed_vcurrent[-1])
                               		
				return self.vv[j0,i0]*self.e1v[j0,i0]
			
			self.zontr = zontr
			self.mertr = mertr

		elif self.mode=='diffusion':
			"""
			Diffusion only. Velocity is maximum but direction is uniformly distributed.
			"""

			active_turtles = np.where(self.pop.state==1)[0]
			self.theta = np.zeros(self.pop.nturtles,dtype='float64')
			for k in active_turtles:	
				self.theta[k] = random.vonmisesvariate(0,0)

			def zontr(j,i,k):
			        #Vectors start at 1 in Fortran, 0 in Python
				turtle = active_turtles[k-1]
				i0 = i-1     
				j0 = j-1
			
				vmax = self.pop.vmax[turtle]
				return (self.uu[j0,i0] + vmax*cos(self.theta[turtle])*self.umask[j0,i0])*self.e2u[j0,i0]
		
		
			def mertr(j,i,k):
				turtle = active_turtles[k-1]
				i0 = i-1    
				j0 = j-1

				vmax = self.pop.vmax[turtle]
				return (self.vv[j0,i0] + vmax*sin(self.theta[turtle])*self.vmask[j0,i0])*self.e1v[j0,i0]

			self.zontr = zontr
			self.mertr = mertr
		

		elif self.mode=='active':
			"""
			Compute a velocity that has both a deterministic and a stochastic component. The norm of the velocity is proportional to the habitat value, the direction is chosen randomly following the Von Mises distribution and depend on the gradient of habitat
			"""
                        active_turtles = np.where(self.pop.state==1)[0]
                        #theta is the direction taken by the turtle
                        self.theta = np.zeros(self.pop.nturtles,dtype='float64')
                        self.theta_rand = np.zeros(self.pop.nturtles,dtype='float64')
                        #grad is the habitat gradient
                        self.grad=np.zeros(self.pop.nturtles,dtype='float64')
 
                        for k in active_turtles:
                                #print k
                                #Compute maximal concentration to define forage habitat

                                PPmax = self.pop.PPmax[k]
                                #Vectors start at 1 in Fortran, 0 in Python
                                i0 = self.pop.i0[k]-1
                                j0 = self.pop.j0[k]-1

                                gradmax=2e-6
                                #Compute habitat and habitat gradient
                                h_left = self.pop.habitat(k,self.Temp[j0,i0-1],self.PP[j0,i0-1],PPmax,Total=True)
                                h_right = self.pop.habitat(k,self.Temp[j0,i0+1],self.PP[j0,i0+1],PPmax,Total=True)
                                h_bot = self.pop.habitat(k,self.Temp[j0-1,i0],self.PP[j0-1,i0],PPmax,Total=True)
                                h_top = self.pop.habitat(k,self.Temp[j0+1,i0],self.PP[j0+1,i0],PPmax,Total=True)
                                h_position = self.pop.habitat(k,self.Temp[j0,i0],self.PP[j0,i0],PPmax,Total=True)

                                xgradh = (h_right*self.umask[j0,i0] - h_left*self.umask[j0,i0-1])/(self.e1v[j0,i0]+self.e1v[j0,i0-1])
                                ygradh = (h_top*self.vmask[j0,i0] - h_bot*self.vmask[j0-1,i0])/(self.e2u[j0,i0]+self.e2u[j0-1,i0])

                                gradh = sqrt(xgradh**2+ygradh**2)
                                self.grad[k]=gradh

                                #Save habitat gradient


				if xgradh>0:
					theta0 = atan(ygradh/xgradh)
				elif xgradh<0:
					theta0 = pi + atan(ygradh/xgradh)
				else:
					theta0 = np.sign(ygradh)*pi/2
				theta0 = theta0%(2*pi)
				#alpha is a parameter representing the easiness that turtles have to follow the habitat gradient, the lower alpha is, the higher the diffusion will be

				#self.theta[k] = random.vonmisesvariate(theta0,alpha*gradh)
				self.theta[k] = random.vonmisesvariate(theta0,alpha*gradh)
                                #Save habitat
                                habT, habPP = self.pop.habitat(k,self.pop.Tenv_now[k],self.pop.PP_now[k],PPmax,Total=False)
                                hab = h_position
                                self.pop.hab[self.pop.tind[k],k] = hab
                                self.pop.habT[self.pop.tind[k],k] = habT
                                self.pop.habPP[self.pop.tind[k],k] = habPP
                                self.pop.xgrad[self.pop.tind[k],k] = xgradh
                                self.pop.ygrad[self.pop.tind[k],k] = ygradh

                                vect_speed_uswim = []
                                vect_speed_ucurrent = []
                                vect_speed_vswim = []
                                vect_speed_vcurrent = []

                        def zontr(j,i,k,side,ti,sign,time_drifting,TotalTime,paire,saving):

			        #Vectors start at 1 in Fortran, 0 in Python
				turtle = active_turtles[k-1]
				#print 'TURTLE ', turtle
				i0 = i-1     
				j0 = j-1
			        ti = ti-1
				side = side - 1 
				h = self.pop.hab[self.pop.tind[turtle],turtle]
                                vmax = self.pop.vmax[turtle]
                                
                                if self.umask[j0,i0] == 0 :
					#print "TI",ti,"SIDE",side,"I0",i0,"J0",j0
                                        h = self.pop.habitat(turtle,self.Temp[j0,ti],self.PP[j0,ti],PPmax)
                                        #print "UMASK NULL", self.umask[j0,i0]
                                        xgradh =sign*(h - self.pop.habitat(turtle,self.Temp[j0,side],self.PP[j0,side],PPmax,Total=True)*self.umask[j0,i0])/self.e1v[j0,i0]
                                        #print 'XGRAD', xgradh
                                        if (xgradh/gradmax) > 0 :
                                                xgrad=min((xgradh/gradmax),1)
                                        elif (xgradh/gradmax) < 0:
                                                xgrad=max((xgradh/gradmax),-1)
                                        else:
                                                xgrad=0.0
                                        xswim = vmax*xgrad

                                else :
                                	xswim = vmax*(1-h)*cos(self.theta[turtle])

                                	#xswim = vmax*cos(self.theta[turtle])
                                #print 'XSWIM',xswim
                                
                                #xswim = vmax*(1-h)*cos(self.theta[turtle])

                                if paire == 1:
                               		vect_speed_uswim.append(xswim)
                               		vect_speed_ucurrent.append(self.uu[j0,i0])
                               		#print 'paire1', self.uu[j0,i0],xswim*self.umask[j0,i0]
                                elif paire == 2 :
                               		vect_speed_uswim[-1]=(vect_speed_uswim[-1]+xswim)/2.
                               		#print 'CURRENT', ratio_posx,vect_speed_ucurrent[n-1],self.uu[j0,i0],self.uu[j0,i0]-vect_speed_ucurrent[n-1],vect_speed_ucurrent[n-1]+ratio_posx*(self.uu[j0,i0]-vect_speed_ucurrent[n-1])
                               		vect_speed_ucurrent[-1]=(vect_speed_ucurrent[-1] + self.uu[j0,i0])/2.
                               		#print 'paire2', self.uu[j0,i0],xswim*self.umask[j0,i0]
                            				   

                                if saving == 1:
                               		self.pop.u_current[self.pop.tind[turtle],turtle] = self.pop.u_current[self.pop.tind[turtle],turtle]+((time_drifting/86400.)*vect_speed_ucurrent[-1])
                               		self.pop.u_swim[self.pop.tind[turtle],turtle] = self.pop.u_swim[self.pop.tind[turtle],turtle]+((time_drifting/86400.)*vect_speed_uswim[-1])
                               		self.pop.u_tot[self.pop.tind[turtle],turtle] = self.pop.u_current[self.pop.tind[turtle],turtle] +  self.pop.u_swim[self.pop.tind[turtle],turtle]  
                               		#print 'time_drifiting ',self.pop.u_tot[self.pop.tind[turtle],turtle]*86.4#,  time_drifting/86400., (time_drifting/86400.)*vect_speed_ucurrent[-1], self.pop.u_current[self.pop.tind[turtle],turtle], (time_drifting/86400.)*vect_speed_uswim[-1],self.pop.u_swim[self.pop.tind[turtle],turtle]

				return (self.uu[j0,i0] + xswim)*self.e2u[j0,i0]
		
		
			def mertr(j,i,k,side,tj,sign,time_drifting,TotalTime,paire,saving):

				turtle = active_turtles[k-1]
				i0 = i-1    
				j0 = j-1
			        tj =tj-1
                                side = side -1
                                h = self.pop.hab[self.pop.tind[turtle],turtle]
				vmax = self.pop.vmax[turtle]
                                
                                if self.vmask[j0,i0]==0:
					#print "TJ",tj,"SIDE",side,"I0",i0,"J0",j0
                                        h = self.pop.habitat(turtle,self.Temp[tj,i0],self.PP[tj,i0],PPmax)
                                        #print "VMASK NULL",self.vmask[j0,i0]
                                        ygradh = sign*(h - self.pop.habitat(turtle,self.Temp[side,i0],self.PP[side,i0],PPmax,Total=True)*self.vmask[j0,i0])/self.e2u[j0,i0]
                                        #print 'YGRAD',ygradh

                                        if (ygradh/gradmax) > 0 :
                                                ygrad=min((ygradh/gradmax),1)
                                        elif (ygradh/gradmax) < 0:
                                                ygrad=max((ygradh/gradmax),-1)
                                        else:
                                                ygrad=0.0
                                        yswim = vmax*ygrad
                                
                                else :
                                        yswim = vmax*(1-h)*sin(self.theta[turtle])
                                        #yswim = vmax*sin(self.theta[turtle])
                                
                                #yswim = vmax*(1-h)*sin(self.theta[turtle])
                                #print 'YSWIM',yswim
                                if paire == 1:
                               		vect_speed_vswim.append(yswim)
                               		vect_speed_vcurrent.append(self.vv[j0,i0])
                                elif paire == 2 :
                               		#print 'CURRENT', ratio_posy,vect_speed_vcurrent[n-1],self.vv[j0,i0],self.vv[j0,i0]-vect_speed_vcurrent[n-1],vect_speed_vcurrent[n-1]+ratio_posy*(self.vv[j0,i0]-vect_speed_vcurrent[n-1])
                               		vect_speed_vswim[-1]=(vect_speed_vswim[-1] +yswim)/2.
                               		vect_speed_vcurrent[-1]=(vect_speed_vcurrent[-1] +self.vv[j0,i0])/2.

                                if saving == 1:

                               		self.pop.v_current[self.pop.tind[turtle],turtle] = self.pop.v_current[self.pop.tind[turtle],turtle]+((time_drifting/86400.)*vect_speed_vcurrent[-1])
                               		self.pop.v_swim[self.pop.tind[turtle],turtle] = self.pop.v_swim[self.pop.tind[turtle],turtle]+((time_drifting/86400.)*vect_speed_vswim[-1])
                               		self.pop.v_tot[self.pop.tind[turtle],turtle] = self.pop.v_current[self.pop.tind[turtle],turtle] +  self.pop.v_swim[self.pop.tind[turtle],turtle]  
                                	#print 'saving V',time_drifting/86400.,self.pop.v_current[self.pop.tind[turtle],turtle], self.pop.v_swim[self.pop.tind[turtle],turtle],  self.pop.v_tot[self.pop.tind[turtle],turtle]               
 



				#print 'mertr',(self.vv[j0,i0] + yswim)
				return (self.vv[j0,i0] + yswim)*self.e1v[j0,i0]
		

			self.zontr = zontr
			self.mertr = mertr





	def Advection2D(self,deltaT):
		"""
		Compute advection of the population by zonal and meridional currents of transport zontr & mertr during deltaT seconds
		
		deltaT : f, time step in seconds
		"""
		active_turtles = np.where(self.pop.state==1)[0]
		#Switch to Fortran array
		x = self.pop.x[active_turtles]+1
		y = self.pop.y[active_turtles]+1

		i0 = self.pop.i0[active_turtles]+1
		j0 = self.pop.j0[active_turtles]+1
		state = self.pop.state[active_turtles]
		
		adv.advect2d(x,y,i0,j0,state,deltaT,self.grid.tfac,self.grid.mask,self.grid.pivot,self.grid.key_periodic,self.grid.overlap,self.grid.key_jfold,self.zontr,self.mertr)

		self.pop.state[active_turtles] = state

		#For still active turtles
		active_turtles = np.where(self.pop.state==1)[0]
		#Switch to Python array
		self.pop.x[active_turtles] = x[state==1]-1
		self.pop.y[active_turtles] = y[state==1]-1

		self.pop.i0[active_turtles] = i0[state==1]-1
		self.pop.j0[active_turtles] = j0[state==1]-1

		#Increment time
		self.pop.tind[active_turtles]=self.pop.tind[active_turtles]+1
		self.pop.age_now[active_turtles] = self.pop.age_now[active_turtles]+deltaT/86400


		self.pop.age[self.pop.tind[active_turtles],active_turtles] = self.pop.age_now[active_turtles]

		#turtles in state 2 (out of domain) stay at the same position
		turtles_out = np.where(self.pop.state==2)[0]  
		for j in turtles_out:
			i = self.pop.tind[j]
			self.pop.trajlon[i+1:,j] = self.pop.trajlon[i,j]
			self.pop.trajlat[i+1:,j] = self.pop.trajlat[i,j]
			if self.pop.key_alltracers:
				self.pop.Tenv[i+1:,j] = self.pop.Tenv[i,j]
				self.pop.PP[i+1:,j] = self.pop.PP[i,j]

