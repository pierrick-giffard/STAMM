#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Define particle class with useful variables.
"""


from parcels import JITParticle, Variable
import numpy as np
from operator import attrgetter

def define_Turtle_Class(fieldset, param):
    class turtle(JITParticle):
        #Sampling
        if param['key_alltracers']:
            T = Variable('T', initial=fieldset.T)
            NPP = Variable('NPP', initial=fieldset.NPP)
        #Passive
        distance = Variable('distance', to_write=True, initial=0., dtype=np.float32)
        lat_dist = Variable('lat_dist', to_write=False, initial=0., dtype=np.float32)
        lon_dist= Variable('lon_dist', to_write=False, initial=0., dtype=np.float32)
        prev_lon = Variable('prev_lon', to_write=True, dtype=np.float32, initial=attrgetter('lon'))#tmp
        prev_lat = Variable('prev_lat', to_write=True, dtype=np.float32, initial=attrgetter('lat'))#tmp
        u_current = Variable('u_current', to_write=True, dtype=np.float32)
        v_current = Variable('v_current', to_write=True, dtype=np.float32)
        #Active
        SCL = Variable('SCL', to_write=False, dtype=np.float32)
        M = Variable('M', to_write=False, dtype=np.float32)
        vmax = Variable('vmax', to_write=False, dtype=np.float32)
        age = Variable('age', to_write=True, dtype=np.float32, initial=0.)
        PPmax = Variable('PPmax', to_write=False, dtype=np.float32)       
        habT = Variable('habT', to_write=True, dtype=np.float32)
        habPP = Variable('habPP', to_write=True, dtype=np.float32)
        hab = Variable('hab', to_write=True, dtype=np.float32)
        theta = Variable('theta', to_write=False, dtype=np.float32)
        u_swim = Variable('u_swim', to_write=True, dtype=np.float32)
        v_swim = Variable('v_swim', to_write=True, dtype=np.float32)
        xgradh = Variable('xgradh', to_write=True, dtype=np.float32)
        ygradh = Variable('ygradh', to_write=True, dtype=np.float32)

    return turtle