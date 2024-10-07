# -*- coding: utf-8 -*-
"""
Created on Wed Jun  8 15:32:26 2022

@author: Tyrtl
"""

""" A module containing a number of commonly used functions by multiple levels of the model """

# Improvements required:
# - Temperature effects on ISA calculations

import numpy as np


# Calculate descent rate of payload under a parachute
def vSteady(mass, Cd, A, density):
    
    sinkRate = np.sqrt((9.81 * mass) / (0.5 * density * Cd * A))
    return sinkRate
    

# Calculate air properties from altitude:
# http://fisicaatmo.at.fcen.uba.ar/practicas/ISAweb.pdf
    
def isa(z, groundTemp, zGround):
    # z is ASL altitude, m
    # groundTemp is temperature at launch site, 0m AGL, C
    # zGround is ASL altitude of 0m AGL at launch site, m
    
    temp0 = 15 # C
    temp0K = temp0 + 273.15 # K
    
    # From ground temperature, determine temperature offset from ISA condition:
    tempMSL = groundTemp + (6.5 * (zGround / 1000)) # Temperature at sea level under these conditions
    dTemp = tempMSL - temp0
    
    temp0K = temp0 + 273.15 # K
    temp = (temp0 + dTemp) - (6.5 * (z/1000)) # C
    
    p0 = 101325 # Pa, ASL
    p = p0 * (1 - (0.0065 * (z / temp0K)))**5.2561 # Pa
    
    R = 287.05 # gas constant, J/kgK
    density = p / (R * (temp + 273.15))
    
    return temp, p, density

def in2m(inches):
    m = 0.0254 * inches
    return m

def sqft2sqm(sqft):
    sqm = 0.092903 * sqft
    return sqm

def oz2kg(oz):
    kg = 0.0283495 * oz
    return kg

def lbf2N(lbf):
    N = 4.44822 * lbf
    return N

def shootingMethod(guess, err): 
    gNew = guess[0] + ((err[0,0]/ (err[0,0] - err[0,1])) * (guess[1] - guess[0]))
    return gNew

def ke(m, v):
    KE = 0.5 * m * (v**2)
    return(KE)

def btAcc(u, s):
    # Assume constant acceleration:
    a = (u**2) / (2 * s)
    return a