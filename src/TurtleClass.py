#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 18:05:20 2020

@author: pgiffard
"""
from parcels import FieldSet, ParticleSet, JITParticle, ParticleFile, plotTrajectoriesFile, Variable,ErrorCode,Field
import numpy as np


def define_Turtle_Class(fieldset):
    class turtle(JITParticle):
        T = Variable('T', initial=fieldset.T)
        u_swim = Variable('u_swim', to_write=True, dtype=np.float32)
        v_swim = Variable('v_swim', to_write=True, dtype=np.float32)
        SCL = Variable('SCL', to_write=True, dtype=np.float32)
        M = Variable('M', to_write=True, dtype=np.float32)
        vmax = Variable('vmax', to_write=True, dtype=np.float32)
        age = Variable('age', to_write=True, dtype=np.float32)
    return turtle