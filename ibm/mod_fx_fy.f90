!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! - mod_fx_fy
!!
!!	Author : Amaury Dehecq (intern at CLS 2011-2012 working with Phillipe Gaspar)
!!	Module description : Convert particles coordinates in grid indices to longitude/latitude coordinates
!!
!!	Inputs : - x,y : array, positions in grid indices for each particle (x->U grid, y -> V grid
!!		 - glamu, gphiv : 2D array, matrix of longitude, resp. latitude at each point of the grid
!!               - key_periodic : boolean, is the input grid east/west periodic?
!!               - overlap : int, number of overlapping points in east/west periodicity
!!               - xdim, ydim : int, x & y dimension of uu/vv (optionnal if used with f2py
!!               - npart : int, length of x & y vectors (optionnal if used with f2py)
!!
!!	Outputs :-new position xnew,ynew
!! 
!!      Last modified : 13/07/2012
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE FX_FY
  
CONTAINS
  
  FUNCTION FX(x, y , glamu, key_periodic,overlap,npart, xdim, ydim)
    IMPLICIT NONE
    
    !-arguments-!
    INTEGER, INTENT(IN) :: npart, xdim, ydim, overlap
    !f2py     INTEGER, INTENT(IN) :: npart, xdim, ydim, overlap
    DOUBLE PRECISION, DIMENSION(npart), INTENT(IN) :: x, y
    !f2py     DOUBLE PRECISION, DIMENSION(npart), INTENT(IN) :: x, y
    DOUBLE PRECISION, DIMENSION(ydim,xdim), INTENT(IN) :: glamu
    !f2py     DOUBLE PRECISION, DIMENSION(ydim,xdim), INTENT(IN) :: glamu
    LOGICAL, INTENT(IN) :: key_periodic
    !f2py     LOGICAL, INTENT(IN) :: key_periodic
    
    !_local variables-!
    INTEGER :: i, i1,i2, j1, j2
    DOUBLE PRECISION, DIMENSION(npart) :: fx
    !f2py depend(npart) fx
    DOUBLE PRECISION :: a, b, glamu_sw, glamu_nw, glamu_ne, glamu_se
    
    !============================================================
    !!Function to get the longitude corresponding to x,y position
    !============================================================

    DO i=1,npart
       !===================
       ! Longitude Indices 
       !===================
       i1   = INT(x(i))
       i2   = i1 + 1
       a    = x(i) - AINT(x(i))
       
       IF (key_periodic) THEN
          IF (i2 >= xdim) THEN
             i1   = i1-xdim+overlap
             i2   = i2-xdim+overlap
          ENDIF
       ELSE
          IF (i2 > xdim) i2 = i1
       ENDIF
       
       !==================
       ! Latitude Indices 
       !==================
       !y and j are in V grid, glamu is U grid => switch index in V grid to U grid
       j1 = INT(y(i)+0.5d0)
       j2 = j1 + 1
       b  = (y(i) + 0.5d0) - j1
       IF (j2 > ydim) j2 = j1
       
       !---------------!
       !- Periodicity -!
       !---------------!
       ! independant of where the longitude is cutted and 
       ! longitude conventions: 0/360 or -180/180 or 36.5/396.5
       
       !If glamu_SW > glamu_SE, we are at longitude cut, except sometimes where Orca grid is lightly "twisted" -> +0.5 
       IF (glamu(j1, i1) > glamu(j1, i2)+0.5) THEN   
          glamu_SW=glamu(j1, i1) - 360d0
       ELSE
          glamu_SW = glamu(j1,i1)
       ENDIF
       IF (glamu(j2, i1) > glamu(j2, i2)+0.5) THEN
          glamu_NW=glamu(j2, i1) - 360d0
       ELSE
          glamu_NW = glamu(j2, i1)
       ENDIF
       
       glamu_SE=glamu(j1, i2)
       glamu_NE=glamu(j2, i2)

       fx(i) = (1d0-a) * (1d0-b) * glamu_SW + &
               (1d0-a) *  b      * glamu_NW + &
               a       * (1d0-b) * glamu_SE + &
               a       *  b      * glamu_NE

       IF (fx(i) > 180d0) THEN
          fx(i) = fx(i) - 360d0
       END IF

!       IF (fx(i) < -180d0) THEN
!          fx(i) = fx(i) + 360d0
!       END IF
       

!       if (fx(i)<-180) then
!          print*, i1, j1,glamu_SW, glamu_NW, glamu_SE, glamu_NE, fx(i)
!       endif

    END DO
    
    RETURN
    
  END FUNCTION FX
  
  
  FUNCTION FY(x, y, gphiv, key_periodic, npart, xdim, ydim)
    IMPLICIT NONE
    
    !-arguments-!
    INTEGER, INTENT(IN) :: npart, xdim, ydim
    !f2py     INTEGER, INTENT(IN) :: npart, xdim, ydim
    DOUBLE PRECISION, DIMENSION(npart), INTENT(IN) :: x, y
    !f2py     DOUBLE PRECISION, DIMENSION(npart), INTENT(IN) :: x, y
    DOUBLE PRECISION, DIMENSION(ydim,xdim), INTENT(IN) :: gphiv
    !f2py     DOUBLE PRECISION, DIMENSION(ydim,xdim), INTENT(IN) :: gphiv
    LOGICAL, INTENT(IN) :: key_periodic
    !f2py     LOGICAL, INTENT(IN) :: key_periodic
    
    !- local variables -!
    DOUBLE PRECISION, DIMENSION(npart) :: fy    
    !f2py depend(npart) fy
    INTEGER :: i, i1, i2, j1, j2
    DOUBLE PRECISION :: a, b
    
    !============================================================
    !!Function to get the longitude corresponding to x,y position
    !============================================================

    DO i=1, npart
       !Switch index in U grid to V grid
       i1 = INT(x(i) + 0.5d0)
       i2 = i1 + 1d0
       a  = (x(i) + 0.5d0) - REAL(i1,kind=KIND(1.0d0))
       
       IF (key_periodic) THEN
          IF (i2.GE.xdim) THEN
             i1 = 1
             i2 = 2
          ENDIF
       ELSE
          IF (i2.GT.xdim) i2 = i1
       ENDIF
       
       j1 = INT(y(i))
       j2 = j1 + 1
       b  = y(i) - AINT(y(i))
       IF (j2.GT.ydim) j2=j1
       
       
       fy(i) =  (1d0-a) * (1d0-b) * gphiv(j1,i1) + &
                (1d0-a) *  b      * gphiv(j2,i1) + &
                a       * (1d0-b) * gphiv(j1,i2) + &
                a       * b       * gphiv(j2,i2)
       
    END DO
    
    RETURN 
    
  END FUNCTION FY
     
  
END MODULE FX_FY
