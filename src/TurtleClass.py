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
            T = Variable('T')
            NPP = Variable('NPP')
        #All particles
        active = Variable('active', to_write=True, initial=1, dtype=np.float32)
        prev_lon = Variable('prev_lon', to_write=False, dtype=np.float32, initial=attrgetter('lon'))
        prev_lat = Variable('prev_lat', to_write=False, dtype=np.float32, initial=attrgetter('lat'))
        u_current = Variable('u_current', to_write=True, dtype=np.float32)
        v_current = Variable('v_current', to_write=True, dtype=np.float32)
        age = Variable('age', to_write=True, dtype=np.float64, initial=0.)
        onland = Variable('onland', to_write=False, dtype=np.float32, initial=0.) #number of beachings in a row
        beached = Variable('beached', to_write=False, dtype=np.float32, initial=0.) #0=ocean, 1=onland
        if param['mode'] == 'passive':
            u_swim = Variable('u_swim', to_write=False, dtype=np.float32)
            v_swim = Variable('v_swim', to_write=False, dtype=np.float32)
        #Mortality
        if param['cold_death']:
            lethargy_time = Variable('lethargy_time', to_write=False, dtype=np.float32, initial=0.) #time spent under Tmin
            cold_death = Variable('cold_death', to_write=True, dtype=np.float32, initial=0)
        if param['cold_death'] or param['mode'] == 'active':
            Tmin = Variable('Tmin', to_write=False, dtype=np.float32, initial=fieldset.Tmin)
        #Active
        if param['mode'] == 'active':
            u_swim = Variable('u_swim', to_write=True, dtype=np.float32)
            v_swim = Variable('v_swim', to_write=True, dtype=np.float32)
            M = Variable('M', to_write=False, dtype=np.float32)
            vmax = Variable('vmax', to_write=False, dtype=np.float32)
            PPmax = Variable('PPmax', to_write=True, dtype=np.float32)       
            habT = Variable('habT', to_write=True, dtype=np.float32)
            habPP = Variable('habPP', to_write=True, dtype=np.float32)
            hab = Variable('hab', to_write=True, dtype=np.float32)
            theta = Variable('theta', to_write=False, dtype=np.float32)
            xgradh = Variable('xgradh', to_write=True, dtype=np.float32)
            ygradh = Variable('ygradh', to_write=True, dtype=np.float32)
            SCL = Variable('SCL', to_write=True, dtype=np.float32, initial=fieldset.SCL0)
            Topt = Variable('Topt', to_write=False, dtype=np.float32, initial=fieldset.Topt)
            if param['growth'] == 'Gompertz':
                K = Variable('K', to_write=False, dtype=np.float32, initial=fieldset.K0)
                
    return turtle