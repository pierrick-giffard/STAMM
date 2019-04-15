!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! - Advection2D
!!
!!	Author : Amaury Dehecq (intern at CLS 2011-2012 working with Phillipe Gaspar)
!!	Module description : Advection2D aims at computing advection for a free-drifting particule experiencing ocean currents.
!!	Currents data are given in a C-grid.
!!      Advection scheme inspired from Ariane (http://stockage.univ-brest.fr/~grima/Ariane/ariane.html)
!!
!!      advect2d(x,y,i0_vec,j0_vec,state,deltaT,tfac_mat,mask,pivot,key_periodic,overlap,key_jfold,uu,vv) -> xnew, ynew
!!
!!	Inputs : - x, y : array, positions of all turtles in grid indices (x -> U grid, y -> V grid
!!               - i0_vec, j0_vec : array, coordinates of ambient cell for all turtles (T grid)
!!               - state : array, 1 if active, 2 if out of domain
!!		 - deltaT : float, time step
!!               - tfac_mat : matrix, area of each cell tfac = e1t*e2t
!!               - mask : matrix, land mask
!!               - pivot : char "T" or "F" depending on input grid
!!               - key_periodic : boolean, is the input grid east/west periodic?
!!               - overlap : int, number of overlapping points in east/west periodicity
!!               - key_jfold : boolean, is the input grid north/south periodic?
!!               - xdim, ydim : int, x & y dimension of uu/vv/mask (optionnal if used with f2py
!!               - npart : int, length of x & y vectors (optionnal if used with f2py)
!!		 - external functions uu & vv : zonal and meriodional transports
!!
!!	Outputs :- new positions xnew,ynew
!!
!!      Last modified : 13/07/2012
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE ADVECTION2D


CONTAINS

  SUBROUTINE Advect2D(x,y,i0_vec,j0_vec,state,deltaT,tfac_mat,mask,pivot,key_periodic,overlap,key_jfold,xdim,ydim,npart)
    
    IMPLICIT NONE

    !========================================================
    ! Functions defining the zonal and meridional transport
    ! for each grid point and each turtle
    !========================================================

    !f2py intent(callback) uu, vv
    external :: uu,vv
    double precision :: uu, vv

    !==========================
    !	C-grid parameters
    !==========================
    DOUBLE PRECISION, DIMENSION(ydim,xdim), intent(in) ::  mask, tfac_mat
    !f2py DOUBLE PRECISION, DIMENSION(ydim,xdim), intent(in) :: mask, tfac_mat
    INTEGER, INTENT(IN) :: xdim, ydim, npart, overlap
    
    !F2PY INTEGER, INTENT(IN) :: xdim, ydim, npart, overlap
    CHARACTER(len=1) :: pivot
    !f2py CHARACTER(len=1) :: pivot
    LOGICAL :: key_periodic, key_jfold
    !f2py LOGICAL :: key_periodic, key_jfold
    
    !================================================================
    ! Vectors containing the position and grid indice of each turtle
    !================================================================
    DOUBLE PRECISION, DIMENSION(npart), INTENT(INOUT) :: x,y
    !f2py DOUBLE PRECISION, DIMENSION(npart), INTENT(INOUT) :: x,y
    INTEGER, DIMENSION(npart), INTENT(INOUT) :: i0_vec, j0_vec, state
    !f2py INTEGER, DIMENSION(npart), INTENT(INOUT) :: i0_vec, j0_vec, state

    !============
    ! Time step
    !============
    DOUBLE PRECISION, INTENT(IN) :: deltaT 
    
    !f2py DOUBLE PRECISION, INTENT(IN) :: deltaT

    
   
    !=======================================
    !     Variables used in computation
    !=======================================
    DOUBLE PRECISION :: xnew, ynew
    INTEGER :: i0,j0,il, ir, jb, jt,tileft,tiright,tjbot,tjtop,left,right,bottom,top, k, iswap,paire1,paire2,nopaire,saving
    !INTEGER :: i0,j0,il, ir, jb, jt,k,iswap
    DOUBLE PRECISION :: u1, u2, v1,v2,TotalTime, tfac,time_1,time_2
    
    
    DO k=1,npart
      !WRITE(*,*) "    "
      !WRITE(*,*) "    "
      !WRITE(*,*) "INIT"       
      !WRITE(*,*) "    "
      !WRITE(*,*) "    "
      !WRITE(*,*) "    "

       xnew = x(k)
       ynew = y(k)
       i0 = i0_vec(k)
       j0 = j0_vec(k)
       iswap = state(k)        !set to 2 if particle out of domain
       ir = i0
       il = i0 - 1
       jt = j0
       jb = j0 - 1
       TotalTime=0.  

       !WRITE(*,*)"INIT INDEX", ir,il,jt,jb       
       !WRITE(*,*) "TURTLE ", k
       !======================================================================
       !Boundary cases (north pole, east/west periodicity),Orca OPA-NEMO model
       !======================================================================

       CALL OrcaNorthPoleInteger(ir,jt,jt,pivot,key_jfold,xdim,ydim,iswap)
       CALL OrcaEastWestPeriodicInteger(ir,ir,ir,key_periodic,overlap,xdim,iswap)

         
       DO WHILE (TotalTime<deltaT .and. iswap==1)
          !WRITE(*,*) "CURRENT INDEX",ir,il,jt,jb
          !WRITE(*,*) 'COUNTER',n
          time_1 = TotalTime
          ir = i0
          il = i0 - 1
          jt = j0
          jb = j0 - 1
          
          tileft = i0 - 1
          tiright = i0 + 1
          tjbot = j0 - 1
          tjtop = j0 + 1  
          left = 1
          right = -1
          bottom = 1
          top = -1
          paire1 = 1
          paire2 = 2 
          nopaire = 0
          saving = 0        
          !WRITE(*,*) "u1"
          !WRITE(*,*), 'TEST AVANT'
          u1 = uu(jt,il,k,tileft,i0,left,TotalTime,TotalTime,paire1,saving)  !uu(j0,i0-1)  in Ariane
          !WRITE(*,*) "u1",u1
          !WRITE(*,*),"u2"
          u2 = uu(jt,ir,k,tiright,i0,right,TotalTime,TotalTime,paire2,saving)  !uu(j0,i0)
          !WRITE(*,*) "u2",u2
          !WRITE(*,*),"v1"
          v1 = vv(jb,ir,k,tjbot,j0,bottom,TotalTime,TotalTime,paire1,saving)  !vv(j0-1,i0)
          !WRITE(*,*) "v1",v1
          !WRITE(*,*),"v2"
          v2 = vv(jt,ir,k,tjtop,j0,top,TotalTime,TotalTime,paire2,saving)  !vv(j0,i0)
          !WRITE(*,*), 'TEST APRES', u2
          !WRITE(*,*) "v2",v2
          !WRITE(*,*) "u1"
          !u1 = uu(jt,il,k)  !uu(j0,i0-1)  in Ariane
          !WRITE(*,*),"u2"
          !u2 = uu(jt,ir,k)  !uu(j0,i0)
          !WRITE(*,*),"v1"
          !v1 = vv(jb,ir,k)  !vv(j0-1,i0)
          !WRITE(*,*),"v2"
          !v2 = vv(jt,ir,k)  !vv(j0,i0)

          !WRITE(*,*) 'N',n
          tfac = tfac_mat(j0,i0)
          !WRITE(*,*) "xnew first ", xnew, " ynew first ", ynew

          call advect2d_element(xnew,ynew,i0,j0,deltaT,TotalTime,u1,u2,v1,v2,&
tfac,xdim,ydim,pivot,key_periodic,overlap,key_jfold,iswap)
          time_2 = TotalTime-time_1
          !WRITE(*,*) 'TotalTime', TotalTime,time_2
          saving = 1 
          u1 = uu(jt,il,k,tileft,i0,left,time_2,TotalTime,nopaire,saving)
          v1 = vv(jb,ir,k,tjbot,j0,bottom,time_2,TotalTime,nopaire,saving)


          !WRITE(*,*) "xnew adv", xnew, " ynew adv", ynew
          !=========================================
          !Coast crashes. This should never happen!!
          !=========================================
          
          IF ((mask(j0,i0)==0) .and. (TotalTime<deltaT) .and. iswap==1) THEN
             WRITE(*,*) "coast crash at i0 = ", i0, "j0 = ", j0, "active turtle n°", k
             STOP
          ENDIF
       END DO
       
       x(k)=xnew
       y(k)=ynew
       !WRITE(*,*),'X,Y FORTRAN', x(k),y(k)                          
       
       
       i0_vec(k)=i0
       j0_vec(k)=j0
       state(k) = iswap
    END DO
    
  END SUBROUTINE Advect2D
  




  SUBROUTINE Advect2D_element(x,y,i0,j0,deltaT,TotalTime,u1,u2,v1,v2,tfac,xdim,ydim,pivot,key_periodic,overlap,key_jfold,iswap)
    
    !Compute the elementary displacement until next exit of the cell or end of time step

    IMPLICIT NONE

    !================================================
    !	Indexes of position on C-grid and time
    !================================================
    DOUBLE PRECISION, INTENT(INOUT) :: x,y
    !f2py DOUBLE PRECISION, INTENT(INOUT) :: x,y
    INTEGER, INTENT(INOUT) :: i0, j0, iswap
    !f2py INTEGER, INTENT(INOUT) :: i0, j0, iswap

    !======================
    !	INPUT PARAMETERS
    !======================
    DOUBLE PRECISION, INTENT(IN) :: deltaT
    !f2py DOUBLE PRECISION, INTENT(IN) :: deltaT
    DOUBLE PRECISION, INTENT(INOUT) :: TotalTime
    !f2py DOUBLE PRECISION, INTENT(INOUT) :: TotalTime
    DOUBLE PRECISION, INTENT(IN) :: u1, u2, v1, v2, tfac
    !f2py DOUBLE PRECISION, INTENT(IN) :: u1, u2, v1, v2,tfac

    !u1 & u2 are the zonal transport at the left (resp. right) side of the ambient cell i0
    !v1 & v2 are the meridional transport at the bottom (resp. top) side of the ambient cell j0

    !==========================
    !	C-grid parameters
    !==========================

    INTEGER, INTENT(IN) :: xdim, ydim, overlap
    !F2PY INTEGER, INTENT(IN) :: xdim, ydim, overlap
    CHARACTER(len=1), INTENT(IN) :: pivot
    !f2py CHARACTER(len=1), INTENT(IN) :: pivot
    LOGICAL, INTENT(IN) :: key_periodic, key_jfold
    !F2PY LOGICAL, INTENT(IN) :: key_periodic, key_jfold
    

    !=======================================
    !     Variables used in computation
    !=======================================
    
    DOUBLE PRECISION :: xnew, ynew
    INTEGER :: i1, i2, j1, j2

    DOUBLE PRECISION :: TimeExit, TotalTime_aux,tx,ty, Finterp, Ginterp, GradF, GradG, uin, uout, vin, vout 
    INTEGER :: min_tmp, max_tmp

    
    !==================================
    !Linear interpolation of transports
    !==================================
    
    xnew = x
    ynew = y
    Finterp = u1 + (xnew - REAL(i0-1))*(u2-u1)
    !WRITE(*,*) "Finterp ", Finterp, " u1 ", u1, " u2 ", u2, "xnew", xnew, " i0-1 ", REAL(i0-1)
    Ginterp = v1 + (ynew - REAL(j0-1))*(v2-v1)
    !WRITE(*,*) "Ginterp ", Ginterp, " v1 ", v1, " v2 ", v2, "ynew", ynew, " j0-1 ", REAL(j0-1)
    !========================================================================
    !Indexes i1,i2,j1,j2 determine the C-grid cell where the particule is located
    !i2 determines the side where the particule is going
    !========================================================================
    
    i1=i0-1
    i2=i0
    uin = u1
    uout = u2

    IF (Finterp<0) THEN
       i2=i0-1
       i1=i0
       uin = u2
       uout = u1
    !WRITE(*,*)"uin ", u2, " uout ", u1
    ENDIF
    
    j1=j0-1
    j2=j0
    vin = v1
    vout = v2
    
    IF (Ginterp<0) THEN
       j2=j0-1
       j1=j0
       vin = v2
       vout = v1
    !WRITE(*,*) "vin ",v2, " vout ", v1
    ENDIF
    !WRITE(*,*) 'U SPEED',uin,uout,'V SPEED',vin,vout
    !=====================
    !Transports Gradients
    !=====================
    
    GradF = (uout - uin) * REAL(i2-i1)
    GradG = (vout - vin) * REAL(j2-j1)
    
    !WRITE(*,*) "GradF ", GradF, " GradG ", GradG
    !===========================================
    !Time to exit the cell in x and y directions
    !===========================================
    
    !Tx
    
    IF (Finterp*uout <= 0) THEN
       !WRITE(*,*) "tx = 1e35 "
       tx = 1.e35
    ELSEIF (ABS(GradF/Finterp) <= 1E-11) THEN
       tx = (REAL(i2) - xnew)/Finterp
       !WRITE(*,*) "ABS(GradF/Finterp)<=1E-11 ", tx
    ELSE
       tx = (LOG(ABS(uout)) - LOG(ABS(Finterp)))/GradF
       !WRITE(*,*) "Else", tx,(LOG(ABS(uout)) - LOG(ABS(Finterp))), uout, Finterp
    ENDIF
    
    !Ty
    IF (Ginterp*vout <= 0) THEN
       ty = 1.e35
       !WRITE(*,*) "ty = 1e35 "
    ELSEIF (ABS(GradG/Ginterp) <= 1E-11) THEN
       ty = (REAL(j2) - ynew)/Ginterp
       !WRITE(*,*) "ABS(GradG/Ginterp)<=1E-11 ", ty
    ELSE
       ty = (LOG(ABS(vout)) - LOG(ABS(Ginterp)))/GradG
       !WRITE(*,*) "Else", ty, (LOG(ABS(vout)) - LOG(ABS(Ginterp))), vout,Ginterp
    ENDIF
    
    TimeExit=min(tx,ty)
    !WRITE(*,*) "TX",tx,"TY",ty  
    TotalTime_aux=TotalTime
    TotalTime=TotalTime + TimeExit*tfac
    !WRITE(*,*) "TOTAL TIME", TotalTime,"TIME EXIT",TimeExit,"TIME EXIT X tfac",TimeExit*tfac
    !WRITE(*,*) "TimeExit ", TimeExit, "TotalTime", TotalTime
    !Sometimes, the particle loop undfinetly between 4 corner cells with a very short TimeExit. In that case we make the particle stay in the first cell
    
    !IF (TimeExit<=1e-9) THEN
        !TimeExit=(deltaT-TotalTime_aux)/tfac
        !TotalTime=deltaT

    IF (TimeExit<=1e-10) THEN
        !TimeExit=(deltaT-TotalTime_aux)/tfac
        !WRITE(*,*) "TIME EXIT 1e10"
        TimeExit=0.
        TotalTime=deltaT
    ENDIF

    !Sometimes, t may be negative. If time is negligible (<1s), it is only a numerical precision issue. Otherwise there is a problem and program stops.
    IF (TimeExit<=0) THEN
       IF (ABS(tfac*TimeExit)>1) THEN
          !WRITE(*,*) 't negative'
          !WRITE(*,*) "F : ",tx,uin,uout
          !WRITE(*,*) "G : ",ty,vin,vout
          STOP
       ELSE
          !WRITE(*,*) "t negligible ", TimeExit
          !WRITE(*,*) "TimeExit <= 0"
          TimeExit=0.
          TotalTime=deltaT
          !WRITE(*,*) "TotalTime = deltaT", deltaT
       ENDIF
    ENDIF
    
    !Due to non-divergence of 2D velocity fields, it may happen that particles are trapped inside a cell. In this case, the particle stay in the same cell until the next time step		
    IF (TimeExit>1e34) THEN
       !WRITE(*,*) "Trapped "
       TimeExit=(deltaT-TotalTime_aux)/tfac
       TotalTime=deltaT
    ENDIF
    
    !If TotalTime is greater than deltaT, the particle stops at time deltaT
    IF (TotalTime>deltaT) THEN
       !WRITE(*,*) "stopped "
       TimeExit=(deltaT-TotalTime_aux)/tfac
       TotalTime=deltaT
    ENDIF
    !WRITE(*,*) 'TOTAL TIME',TimeExit*tfac
    !============
    !New position
    !============
    
    IF (tx>TimeExit) THEN
       !WRITE(*,*) "Stay x "
       CALL sub_dont_reachside(xnew,GradF,Finterp,TimeExit,i1,i2)
    ENDIF
    
    IF (ty>TimeExit) THEN
       !WRITE(*,*) "Stay y "
       CALL sub_dont_reachside(ynew,GradG,Ginterp,TimeExit,j1,j2)
    ENDIF
    
    IF (tx<=TimeExit) THEN
       !WRITE(*,*) "Move x ", i1, i2
       xnew=real(i2)
       IF (i2>i1) THEN
          i1=i2
          i2=i2+1
       ELSE
          i1=i2
          i2=i2-1
       ENDIF
       !WRITE(*,*) "Moved x ",xnew,i1, i2
    ENDIF
    
    IF (ty<=TimeExit) THEN
       !WRITE(*,*) "Move y ", j1, j2
       ynew=float(j2)
       IF (j2>j1) THEN
          j1=j2
          j2=j2+1
       ELSE
          j1=j2
          j2=j2-1
       ENDIF
       !WRITE(*,*) "Moved y ",ynew, j1, j2
    ENDIF
          
    !===========================================================
    !Treatment of Orca model singularities, for active particles
    !===========================================================
    IF (TotalTime < deltaT) THEN
       !WRITE(*,*) "ACTIVE"
       i0=max(i1,i2)
       j0=max(j1,j2)
       
       CALL OrcaNorthPoleInteger(i1,j1,j0+1,pivot,key_jfold,xdim,ydim,iswap)
       CALL OrcaNorthPoleInteger(i2,j2,j0+1,pivot,key_jfold,xdim,ydim,iswap)
       CALL OrcaNorthPoleFloat(xnew,ynew,j0+1,pivot,key_jfold,xdim,ydim,iswap)
       CALL OrcaNorthPoleInteger(i0,j0,j0+1,pivot,key_jfold,xdim,ydim,iswap)
 
       min_tmp=min(i1,i2)
       max_tmp=max(i1,i2)
       !A.D : 16/08/2012  Bug at E/W periodicity
       !CALL OrcaEastWestPeriodicFloat(xnew,min_tmp+1,max_tmp,key_periodic,overlap,xdim,iswap)
       !CALL OrcaEastWestPeriodicInteger(i0,min_tmp+1,max_tmp,key_periodic,overlap,xdim,iswap)
       !replaced by:
       CALL OrcaEastWestPeriodicFloat(xnew,min_tmp,max_tmp,key_periodic,overlap,xdim,iswap)
       CALL OrcaEastWestPeriodicInteger(i0,min_tmp,max_tmp,key_periodic,overlap,xdim,iswap)
       
    ENDIF

    x = xnew
    y = ynew
    !WRITE (*,*) 'X Y FORTRAN',x,y          
  END SUBROUTINE Advect2D_element
  
  
  SUBROUTINE sub_dont_reachside(position, grad, interp, TimeExit, imin, imax)
    !===============================================================
    !Compute analytical position of a particle that didn't reach the 
    !side of a cell
    !===============================================================
    INTEGER, INTENT(in)  :: imin, imax
    DOUBLE PRECISION, INTENT(in)  :: grad, interp, TimeExit
    DOUBLE PRECISION, INTENT(inout) :: position
    
    IF (abs(grad*TimeExit) >= 1.e-07) THEN
       position = position + interp * ( EXP(grad*TimeExit) - 1) / grad
       
    ELSE
       position = position + interp * TimeExit * ( 1 + 0.5 * grad * TimeExit )
    ENDIF
    
    IF ((position < MIN(imin,imax)).OR.(position > MAX(imin,imax))) THEN
       print*, 'bad position value...(now corrected):',imin,position,imax
       position=ANINT(position)
       print*, 'New position value:',position,'-',position,interp,TimeExit,grad
    ENDIF
    
  END SUBROUTINE sub_dont_reachside
  
  
  SUBROUTINE OrcaNorthPoleInteger(i,j,jToTest,pivot,key_jfold,xdim,ydim,iswap)
    !=======================================
    !Deal with periodicity at the North pole
    !=======================================
    IMPLICIT NONE
    INTEGER, INTENT(INOUT) :: i,j, iswap
    INTEGER, INTENT(IN) :: jToTest, xdim, ydim
    LOGICAL, INTENT(IN) :: key_jfold
    CHARACTER(len=1) :: pivot
    INTEGER :: factor

    IF ((key_jfold) .and. (jToTest>=ydim)) THEN
       IF (pivot=="T") THEN
          factor=2
       ELSEIF (pivot=="F") THEN
          factor=1
       ENDIF

       i = xdim+factor-i
       j = 2*ydim-factor-j

    ELSEIF ((key_jfold .EQV. .FALSE.) .and. (j>=ydim))  THEN     !out of domain
       iswap = 2
    ENDIF
    
    IF (j<1) THEN            !out of domain
       iswap = 2
    ENDIF

  END SUBROUTINE OrcaNorthPoleInteger
  
  
  SUBROUTINE OrcaNorthPoleFloat(i,j,jToTest,pivot,key_jfold,xdim,ydim,iswap)
    !=======================================
    !Deal with periodicity at the North pole
    !=======================================
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(INOUT) :: i,j
    INTEGER, INTENT(INOUT) :: iswap
    INTEGER, INTENT(IN) :: jToTest, xdim, ydim
    LOGICAL, INTENT(IN) :: key_jfold
    CHARACTER(len=1) :: pivot
    INTEGER :: factor
    
    IF ((key_jfold) .and. (jToTest>=ydim)) THEN
       IF (pivot=="T") THEN
          factor=2
       ELSEIF (pivot=="F") THEN
          factor=1
       ENDIF
       i = REAL(xdim+(factor-1)-i, kind = KIND(1.0d0))
       j = REAL(2*ydim-(factor+1)-j, kind = KIND(1.0d0))

    ELSEIF ((key_jfold .EQV. .FALSE.) .and. (j>=ydim))  THEN     !out of domain
       iswap = 2
    ENDIF
    
    IF (j<1) THEN            !out of domain
       iswap = 2
    ENDIF

  END SUBROUTINE OrcaNorthPoleFloat
  
  
  SUBROUTINE OrcaEastWestPeriodicInteger(i,iToTestWest,iToTestEast,key_periodic,overlap,xdim,iswap)
    !===========================================
    !Deal with east/west periodicity of the grid
    !===========================================
    IMPLICIT NONE
    INTEGER, INTENT(INOUT) :: i, iswap
    INTEGER, INTENT(IN) :: iToTestWest, iToTestEast
    INTEGER, INTENT(IN) :: overlap, xdim
    LOGICAL, INTENT(IN) :: key_periodic

    IF (key_periodic) THEN
       IF (iToTestWest<1) THEN
          i = i + (xdim-overlap)
       ELSEIF (iToTestEast>=xdim) THEN
          i = i - (xdim-overlap)
       ENDIF
    ELSE
       IF (iToTestWest<1) THEN
          iswap = 2
       ELSEIF (iToTestEast>=xdim) THEN
          iswap = 2
       ENDIF 
    ENDIF
  END SUBROUTINE OrcaEastWestPeriodicInteger
  
  
  SUBROUTINE OrcaEastWestPeriodicFloat(i,iToTestWest,iToTestEast,key_periodic,overlap,xdim,iswap)
    !===========================================
    !Deal with east/west periodicity of the grid
    !===========================================
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(INOUT) :: i
    INTEGER, INTENT(INOUT) :: iswap
    INTEGER, INTENT(IN) :: iToTestWest, iToTestEast, xdim, overlap
    LOGICAL, INTENT(IN) :: key_periodic
    
    IF (key_periodic) THEN
       IF (iToTestWest<1) THEN
          i = i + REAL(xdim-overlap, kind=KIND(1.0d0))
       ELSEIF (iToTestEast>=xdim) THEN
          i = i - REAL(xdim-overlap, kind=KIND(1.0d0))
       ENDIF
    ELSE
       IF (iToTestWest<1) THEN
          iswap = 2
       ELSEIF (iToTestEast>=xdim) THEN
          iswap = 2
       ENDIF
    ENDIF
  END SUBROUTINE OrcaEastWestPeriodicFloat
  
END MODULE ADVECTION2D


