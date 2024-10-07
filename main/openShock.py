# -*- coding: utf-8 -*-
"""
Created on Thu Jun  9 12:10:27 2022

@author: Tyrtl
"""

import numpy as np
import fcns as fcns

# Current assumptions:
# The flight path angle is constant throughout inflation (this is because otherwise a system of coupled, nonlinear, nonautonomous 1st order ODES must be solved simultaneously yikeyz)
# 2 DOF simulation and extended Pflanz-Ludtke Method

# Necessary improvements:
# - j needs to change depending on geometry of parachute
# - tf needs to change from initial conditions and parachute scale, geometry etc.
# - coupled nonlinear ODE solver development would allow for arbitrary flight path angle etc.
# To solve the change in velocity and angle concurrently, solve them in a staggered method with a suitably small timestep (no unique solution is directly attainable and this is sufficiently accurate for the time being):

def openShock(mode, rocket, parachute, flight, v0, z0, nt):
    
    # Calculate air density:
    density = fcns.isa(z0, flight.lsTemp, flight.zGround)[2]
    
    # Inflation time, sec:
    tFill = parachute.D0 * parachute.nFill / v0
    
    # Steady descent rate, ve:
    ve = np.sqrt((2 * rocket.mass * 9.81) / (density * (parachute.Cd0 * parachute.S0)))
    
    # Calculate Steady Froude Number, Fre:
    Fre = (ve**2) / (9.81 * parachute.D0)
    
    # Ballistic parameter, A:
    A = Fre / parachute.nFill
    
    # Combined parameter, B:
    B = A * (v0 / ve)**2
    
    j = 6 # inflation exponent (constant, empirical)
    
    # Initial conditions:
    if(mode == "drogue" or mode == "reserve1"):
        gamma = 0
    else:
        gamma = -np.pi / 2
    
    vn0 = 1
    tn0 = 0
    tf = parachute.Cx**(1/j)
    
    vn = np.zeros((1,1))
    vn[0] = vn0
    
    # Nondimensional timestep:
    tn = np.linspace(tn0, tf, nt)
    dtn = tn[1] - tn[0]
    
    # Drag area function: (Cd * S) / (Cd * S)0
    dragAreaFcn = (tn / tf)**j
    
    # RK4 Method (1st order):
    k = np.zeros((1,4))
    y = np.zeros((1,3))
        
    for i in range(0, np.shape(tn)[0]-1):
        
        k[0, 0] = dVn(A, B, gamma, j, tn[i], vn[i])
        y[0, 0] = vn[i] + (k[0, 0] * (dtn/2))
        
        k[0, 1] = dVn(A, B, gamma, j, (tn[i] + (dtn / 2)), y[0, 0])
        y[0, 1] = vn[i] + (k[0, 1] * (dtn/2))
        
        k[0, 2] = dVn(A, B, gamma, j, (tn[i] + (dtn/ 2)), y[0, 1])
        y[0, 2] = vn[i] + (k[0, 2] * dtn)
        
        k[0, 3] = dVn(A, B, gamma, j, (tn[i] + dtn), y[0, 2])
        
        # New nondimensional velocity (at tn[i+1]):
        vnNew = vn[i] + ((dtn/6) * (k[0, 0] + (2 * (k[0, 1] + k[0, 2])) + k[0, 3]))
        vn = np.vstack([vn, vnNew])

    # Use vn and drag area function to calculate the non-dimensional loading:
    vn = np.transpose(vn)
    X = (vn**2) * dragAreaFcn
    
    # Determine maximum, Ck, redimensionalise for openshock load:
    fOpen = (density / 2) * (v0**2) * (parachute.Cd0 * parachute.S0) * X
    
    # Determine the peak load during inflation:
    fOpenShock = np.max(fOpen)
    return fOpenShock, fOpen, tFill

def dVn(A, B, gamma, j, t, y):
    
    
    # Added mass term, ma: (define parachute height, radius, density from altitude during deployment, porosity)
    #p = 0.2 # porosity
    #ka = 1.068 * (1 - (1.465 * p) - (0.25975 * p**2) + (1.2626* p**3))
    #ma = ka * density * (4/3) * np.pi * (parachute.D0/2)**2 * (parachute.D0 / 2)    
    
    # Assume constant angle (for now!)
    # Implement the added mass effects here somehow (rederive dynamics)
    dy = ((-1/B) * np.sin(gamma)) - ((1/A) * (t**j) * (y**2))
    
    return(dy)