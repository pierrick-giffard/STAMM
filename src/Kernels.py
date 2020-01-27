#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
All kernels needed.
"""
from parcels import FieldSet, ParticleSet, JITParticle, ParticleFile, plotTrajectoriesFile, Variable,ErrorCode,Field
import numpy as np


def RK4():
    a=1

def AgeParticle(particle, fieldset, time, dt):
    particle.age += dt


def SampleP(particle, fieldset, time): 
    particle.T = fieldset.T[time, particle.depth, particle.lat, particle.lon]
    
def compute_SCL(particle, fieldset, time, mode, species):
    if mode == 'active':
        if species == 'leatherback':
            particle.SCL = 1.43*(1-np.exp(-0.226*(particle.age/365.+0.17)))
        elif species == 'green':
            particle.SCL = 0.983*(1-np.exp(-0.949*(particle.age/365.+0.074)))
        elif species == 'loggerhead':
            particle.SCL = 0
    
def compute_M(particle, fieldset, time, mode, species):
    if mode == 'active':
        if species == 'leatherback':
            particle.M = 112.31*(particle.SCL)**2.86
        else:
            particle.M = 0