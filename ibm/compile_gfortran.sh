#!/bin/bash

# MACOSX 

#Run to compile Fortran-coded programs
rm -f *.so *.mod *.pyc || exit 1
f2py -c --fcompiler=gfortran -m mod_fx_fy mod_fx_fy.f90 || exit 1
f2py -c --fcompiler=gfortran -m Advection2D_f Advection2D.f90 || exit 1


