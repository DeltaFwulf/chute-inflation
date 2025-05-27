import numpy as np
from math import sqrt, exp, sin

from parachute import *

"""Estimates the highest load imparted by a parachute during inflation using the extended Pflanz-Ludtke method."""

def openShockPL(parachute:Parachute, mSystem:float, vStretch:float, airDensity:float, angInit:float, extended:bool=True):

    """
    Uses the Extended Pflanz-Ludtke method to estimate open shock loads

    Please do not use extended method when system has a low ballistic parameter, A, as well as a low value of cx, 
    this leads to overshoots in the range of 0.1<=A<=10 when j == 6.
    """

    cds0 = parachute.cd * parachute.Aref # steady drag area
    vTerm = sqrt(2 * mSystem * 9.80665 / (airDensity * cds0))
    vRatio = vStretch / vTerm
   
    # Pflanz-Ludtke Variables #################################################################################################
    tFill = parachute.tFill(vStretch)
    nFill = tFill * vStretch / parachute.d0 # non-dimensional time
    A = vTerm**2 / (9.81 * parachute.d0 * nFill)
    B = A * (vStretch / vTerm)**2 # a combined parameter

    # Extended Pflanz-Ludtke Method ###########################################################################################
    j = parachute.infExp
    
    tm = tFill * ((j * (j + 1) * A) / (j+2))**(1/(j+1)) # time at which maximum load occurs # FIXME: this is > tMax, (seems to occur when vehicle continues to accelerate beyond tMax?)
    tMax = tFill * parachute.cx ** (1/j)
    Alx = (j+2) / (j * (j+1)) * parachute.cx**((j+1)/j)

    if A < Alx: # tm < tMax, vehicle decelerates significantly throughout inflation
        ck = ((j+2) / (2 * (j+1)))**2 * ((j*(j+1)*A)/(j+2))**(j/(j+1))
    else: # peak clk occurs at tm = tMax
        ck = (1 + (1 / (A*(j+1))) * parachute.cx**((j+1)/j))**-2 * parachute.cx

    if extended:
        C1 = sqrt(parachute.infExp) * (vTerm/vStretch)**2 * exp(-B)
        C2 = sqrt(parachute.infExp) * (vTerm/vStretch)**2 * (1 - exp(-B)) * sin(-angInit) * exp(-(A/6)*parachute.infExp**0.25)

        ck += ck + C1 + C2

    peakLoad = ck * 0.5 * airDensity * vStretch**2 * cds0

    return peakLoad, tFill, tm, tMax, A, vRatio


def openShock2D():
    """Uses a nondimensionalised 2D inflation simulation to estimate peak load"""
    pass