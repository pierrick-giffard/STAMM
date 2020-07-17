# -*- coding: utf-8 -*-
"""
Plot swimming velocity roses.
USE: python SwimmingRose.py output.nc
"""

# =============================================================================
# IMPORTS
# =============================================================================
import numpy as np
from netCDF4 import Dataset
import sys, os


#Personal librairies
sys.path.insert(1, os.path.join(sys.path[0], '../../LIB'))
import plot_lib as pl
import librumeau as brum

# =============================================================================
# AUTO PARAMETERS
# =============================================================================
tmax = -1 # last time
bins_classes = [0,0.05,0.1,0.2,0.3,0.5,1,2]
title = False


# =============================================================================
# OVERWRITE PARAMETERS
# =============================================================================
#tmax = 48 # last time in number of t_output
bins_classes = [0,0.05,0.1,0.2,0.3,0.5,1,2]
title = 'Swimming against wave'

# =============================================================================
# DATA
# =============================================================================
file = sys.argv[1]
outfile = file.replace(".nc","_SwimmingRose.png")

nc = Dataset(file)
us = nc.variables['u_swim'][:tmax, :]
vs = nc.variables['v_swim'][:tmax, :]

# =============================================================================
# CALCULATION
# =============================================================================

# Speed norm
speed = brum.compute_speed(us, vs)

# Direction
theta = np.arctan2(vs, us)

#Plot
pl.swimming_rose(theta, speed,  bins_classes,  save=outfile, title=title)
