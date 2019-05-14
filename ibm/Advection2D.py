#!/usr/bin/env python
#-*- coding:utf-8 -*-
"""
	Author : Romain Piton (intern at CLS Summer 2011 working with Phillipe Gaspar)
	Module description : Advection 2D aims at computing advection for a free-drifting particule experiencing ocean currents.
	Currents data are given in an irregular C-grid.
	Inputs : -positon x,y
		 -time t
		 -Time step deltaT
		 -Transports : zontr1,zontr2,mertr1,mertr2
	Outputs :-new position xnew,ynew
"""



import numpy as np
import math

def OrcaNorthPoleInteger(i,j,jToTest,pivot,key_jfold,xdim,ydim):
	"""
	Deal with periodicity at the North pole
	"""
		
	if ((key_jfold)&(jToTest>=ydim)):
		if pivot=='T':
			factor=2
		elif pivot=='F':
			factor=1
		return xdim+factor-i,2*(ydim)-factor-j
	else:
		return i,j


def OrcaNorthPoleFloat(i,j,jToTest,pivot,key_jfold,xdim,ydim):
	
	if ((key_jfold)&(jToTest>=ydim)):
		if pivot=='T':
			factor=2
		elif pivot=='F':
			factor=1
		return float(xdim+(factor-1)-i),float(2*ydim-(factor+1)-j)
	else:
		return i,j

def OrcaEastWestPeriodic(i,iToTestWest,iToTestEast,key_periodic,xdim):
	"""
	Deal with east/west periodicity of the grid
	"""	
	if (key_periodic):
		if (iToTestWest<1):
			return i + (xdim-2)
		elif (iToTestEast>=xdim):
			return i - (xdim-2)
		else:
			return i
	else:
		return i



def DontReachSide(position,grad,interp,TimeExit,index1,index2):
	"""
	Compute position of a particle that didn't reach the side of the cell
	"""
	
	if (abs(grad*TimeExit) >= 1e-07):
            new=position + interp*(np.exp(grad*TimeExit) - 1) / float(grad)
            
	else:
            new=position + interp*TimeExit*float(1+0.5*grad*TimeExit)
	
	if ((new<float(min(index1,index2))) or (new>float(max(index1,index2)))):
            print("bad value position corrected")
            new=np.round(new)
	
	return new
		


def Advect2D(x,y,i0,j0,deltaT,zontr,mertr,mesh,turtle):
	
	"""
	Indexes of position on C-grid and time
	"""

	xnew,ynew=np.float64(x),np.float64(y)	
	TotalTime,tx,ty =0., 0., 0.

	"""
	C-grid parameters
	"""
	
	pivot='T'
	key_periodic=True
	key_jfold=False
	xdim = mesh.glamu.shape[1]
	ydim = mesh.glamu.shape[0]
	
	"""
	Boundary cases (Repliment meridient, periodicite est_ouest), modèle Orca OPA-NEMO
	"""
	
	i0,j0=OrcaNorthPoleInteger(i0,j0,j0,pivot,key_jfold,xdim,ydim)
	i0=OrcaEastWestPeriodic(i0,i0,i0,key_periodic,xdim)

	while (TotalTime<deltaT):

		"""
		Linear interpolation of transports
		"""

                Finterp = np.float64(zontr(j0-1,i0-2,turtle) + (xnew-float(i0-1))*(zontr(j0-1,i0-1,turtle)-zontr(j0-1,i0-2,turtle)))
                Ginterp = np.float64(mertr(j0-2,i0-1,turtle) + (ynew-float(j0-1))*(mertr(j0-1,i0-1,turtle)-mertr(j0-2,i0-1,turtle)))

		"""
		Indexes i1,i2,j1,j2 determine the C-grid cell where the particule is located
		i2 determines the side where the particule is going
		"""
                
		i1=i0-1
		i2=i0

		if (Finterp<0):
			i2=i0-1
			i1=i0

		j1=j0-1
		j2=j0

		if (Ginterp<0):
			j2=j0-1
			j1=j0
                
				
		"""
		Transports Gradients
		"""
                
		GradF =  np.float64((zontr(j0-1,i2-1,turtle)-zontr(j0-1,i1-1,turtle))*float(i2-i1))
		GradG =  np.float64((mertr(j2-1,i0-1,turtle)-mertr(j1-1,i0-1,turtle))*float(j2-j1))

		tfac=mesh.tfac[j0-1,i0-1]

		"""
		Time to exit the cell in x and y directions
		"""
                #tx
                if (Finterp*zontr(j0-1,i2-1,turtle) <= 0):
                    tx = float(1E35)
		elif (abs(GradF/Finterp) <= 1E-11):
                    tx = np.float64(float(i2-xnew)/float(Finterp))
		else:
                    tx = np.float64(float(math.log(abs(zontr(j0-1,i2-1,turtle)))- math.log(abs(Finterp)))/float(GradF))


                #ty
		if (Ginterp*mertr(j2-1,i0-1,turtle) <= 0):
                    ty = float(1E35)
		elif (abs(GradG/Ginterp) <= 1E-11):
                    ty = np.float64(float(j2-ynew)/float(Ginterp))
		else:
                    ty = np.float64(float(math.log(abs(mertr(j2-1,i0-1,turtle)))- math.log(abs(Ginterp)))/float(GradG))

                                
		TimeExit=min(tx,ty)
		TotalTime_aux=TotalTime
		TotalTime=TotalTime + TimeExit*tfac

		#Sometimes, t may be negative. If time is negligible (<1s), it is only a numerical precision issue. Otherwise there is a problem and program stops.
		if (TimeExit<=0):
			if (abs(tfac*TimeExit)>1):
				print('t negative')
				print("F : ",tx,zontr(j0-1,i2-1,turtle),zontr(j0-1,i1-1,turtle), Finterp)
				print("G : ",ty,mertr(j2-1,i0-1,turtle),mertr(j1-1,i0-1,turtle))
				quit()
			else:
				TimeExit=0.
				TotalTime=deltaT
		
		#Due to non-divergence of 2D velocity fields, it may happen that particles are trapped inside a cell. In this case, the particle stay at the same position until the next time step		
		if TimeExit>1e34:
			TimeExit=(deltaT-TotalTime_aux)/tfac
			TotalTime=deltaT
				

		if ((TotalTime>deltaT)):
			TimeExit=(deltaT-TotalTime_aux)/tfac
			TotalTime=deltaT
		

		"""
		New position
		"""
                

		if (tx>TimeExit):
                    xnew=DontReachSide(xnew,GradF,Finterp,TimeExit,i1,i2)
                    
                if (ty>TimeExit):
                    ynew=DontReachSide(ynew,GradG,Ginterp,TimeExit,j1,j2)

                if (tx<=TimeExit):
			xnew=float(i2)
			if (i2>i1):
				i1=i2
				i2=i2+1
			else:
				i1=i2
				i2=i2-1
		if (ty<=TimeExit):
			ynew=float(j2)
			if (j2>j1):
				j1=j2
				j2=j2+1
			else:
				j1=j2
				j2=j2-1

		"""
		Treatment of Orca model singularities, for active particles
		"""
		if TotalTime < deltaT:
		    i0=max(i1,i2)
		    j0=max(j1,j2)

		    i1,j1=OrcaNorthPoleInteger(i1,j1,j0+1,pivot,key_jfold,xdim,ydim)
                    i2,j2=OrcaNorthPoleInteger(i2,j2,j0+1,pivot,key_jfold,xdim,ydim)
                    xnew,ynew=OrcaNorthPoleFloat(xnew,ynew,j0+1,pivot,key_jfold,xdim,ydim)
                    i0,j0=OrcaNorthPoleInteger(i0,j0,j0+1,pivot,key_jfold,xdim,ydim)
                    
		    min_tmp=min(i1,i2)
		    max_tmp=max(i1,i2)
		    xnew=OrcaEastWestPeriodic(xnew,min_tmp,max_tmp,key_periodic,xdim)
                    i0=OrcaEastWestPeriodic(i0,min_tmp,max_tmp,key_periodic,xdim)
                    
		#Temporarily for particles reaching north boundary
		if ynew>=ydim:
			print("Turtle reached north boundary, shifting to one cell south")
			ynew = ynew -1
			j0 = j0-1

		
		"""
		Coast crashes. This should never happen!!
		"""
		if (mesh.mask[j0-1,i0-1]==0) & (TotalTime<deltaT):
			print("coast crash turtle n. : ", turtle, "at x = ", xnew, ", y = ", ynew, "i1 = ", i1, "i2 = ", i2, "j1 = ", j1, "j2 = ", j2)
			print("TimeExit = ", TimeExit, Finterp, Ginterp)
			quit()

	
	return xnew,ynew, i0, j0
				


