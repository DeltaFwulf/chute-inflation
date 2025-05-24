# -*- coding: utf-8 -*-
"""
Created on Tue Jun  7 22:49:57 2022

@author: Tyrtl
"""

import numpy as np

def snatchLoad(rocket, bag, nosecone, lines, vRocket, density, mode, drogue):
    
    # if on drogue, compare Cd*A of nosecone and drogue, simulate larger one:
    if(mode == "drogue" and (bag.Cd * bag.A) < (nosecone.Cd * nosecone.A)):
        projectile = nosecone
        
    else:
        projectile = bag
        
    # include gravity?
    if(mode == "drogue"):
        g = 0
    else:
        g = 9.81
    
    # Timestep, s:
    dt = 0.001
    
    # Separation motion to beginning of linestretch:
    # Implement a shooting method for accurate length solution later
    (x, vRocket) = motion(projectile, dt, g, density, lines.length, vRocket, drogue, mode)
    
    vStretch = x[np.shape(x)[1], 1]
    
    #print("Rocket Velocity at line stretch: ", vRocket, " m/s")
    vCommon = ((rocket.mass * vRocket) + (projectile.mass * (vRocket - vStretch))) / (rocket.mass + projectile.mass)
    
    # Energy method of predicting shock loads:
    fSnatch = vStretch * np.sqrt(lines.k * (rocket.mass * projectile.mass) / (rocket.mass + projectile.mass)) # N
    #print("Max separation length: ", x[np.shape(x)[0] - 1, 0], " m")
    #print("Relative velocity at line stretch: ", x[np.shape(x)[0] - 1, 1], "m/s")
    
    return fSnatch, vCommon
    
def motion(projectile, dt, g, density, lineLength, vRocket, drogue, mode):
    
    # Simulate motion of projectile until line stretch:
    x = np.zeros((1,3), dtype=float)
    j = np.zeros((1,4), dtype=float)
    k = np.zeros((1,4), dtype=float)
    
    i = 0
    
    # Stop at line stretch:
    while(x[i, 0] < lineLength):
        
        # New State Vector:
        xNew = np.zeros((1,3))
        
        # Runge-Kutta 4 Method:
        if(mode == "main"):
            
            j[0,0] = dt * fcnSepMain(density, projectile, x[i], vRocket, g, drogue)
            j[0,1] = dt * fcnSepMain(density, projectile, [x[i,0], x[i,1] + (j[0,0]/2), x[i,2]], vRocket, g, drogue)
            j[0,2] = dt * fcnSepMain(density, projectile, [x[i,0], x[i,1] + (j[0,1]/2), x[i,2]], vRocket, g, drogue)
            j[0,3] = dt * fcnSepMain(density, projectile, [x[i,0], x[i,1] + j[0,2], x[i,2]], vRocket, g, drogue)
            
            xNew[0, 2] = fcnSepMain(density, projectile, x[i,:], vRocket, g, drogue)
            
        elif(mode == "drogue" or mode == "reserve1"):
            
            j[0,0] = dt * fcnSepDR(density, projectile, x[i], vRocket)
            j[0,1] = dt * fcnSepDR(density, projectile, [x[i,0], x[i,1] + (j[0,0]/2), x[i,2]], vRocket)
            j[0,2] = dt * fcnSepDR(density, projectile, [x[i,0], x[i,1] + (j[0,1]/2), x[i,2]], vRocket)
            j[0,3] = dt * fcnSepDR(density, projectile, [x[i,0], x[i,1] + j[0,2], x[i,2]], vRocket)
            
            xNew[0, 2] = fcnSepDR(density, projectile, x[i,:], vRocket)
        
        elif(mode == "reserve2"):
            
            j[0,0] = dt * fcnSepReserve2(density, projectile, x[i], vRocket, g)
            j[0,1] = dt * fcnSepReserve2(density, projectile, [x[i,0], x[i,1] + (j[0,0]/2), x[i,2]], vRocket, g)
            j[0,2] = dt * fcnSepReserve2(density, projectile, [x[i,0], x[i,1] + (j[0,1]/2), x[i,2]], vRocket, g)
            j[0,3] = dt * fcnSepReserve2(density, projectile, [x[i,0], x[i,1] + j[0,2], x[i,2]], vRocket, g)
            
            xNew[0, 2] = fcnSepReserve2(density, projectile, x[i,:], vRocket, g)
            
        else:
            print("Incorrect mode input somewhere bucko!")
            pass
            
        k[0,0] = dt * x[i,1]
        k[0,1] = dt * (x[i,1] + (j[0,0] / 2))
        k[0,2] = dt * (x[i,1] + (j[0,1] / 2))
        k[0,3] = dt * (x[i,1] + j[0,2])
        
        xNew[0, 0] = x[i,0] + ((1/6) * ((k[0,0] + k[0,3] + 2 * (k[0,1] + k[0,2]))))
        xNew[0, 1] = x[i,1] + ((1/6) * ((j[0,0] + j[0,3] + 2 * (j[0,1] + j[0,2]))))
        
        x = np.vstack([x, xNew])
        
        # Update rocket velocity:
        if(mode == "main"):
            # Allow rocket to accelerate under gravity over this period:
            vRocket += (dt * g)
        
        i += 1
        
    return x, vRocket
        
# Assumption: rocket does not accelerate during separation.        
def fcnSepDR(density, projectile, x, vRocket):
    
    vFreestream = vRocket - x[1]
    acc = (0.5 / projectile.mass) * ((density * (vFreestream**2) * projectile.Cd * projectile.A))
    
    return acc


def fcnSepMain(density, projectile, x, vRocket, g, drogue):
    
    vFreestream = vRocket - x[1]
    pDyn = 0.5 * density * (vFreestream**2)
    
    drag = pDyn * ((projectile.Cd * projectile.A) + (drogue.Cd0 + drogue.S0))
    
    # Added term g accounts for rocket in freefall, accelerating away from the decelerating projectile:
    acc = ((drag - (projectile.mass * g)) / projectile.mass) + g
    
    return acc

def fcnSepReserve2(density, projectile, x, vRocket, g):
    # The rocket is vertical, however is not accelerating as under drogue, also no parachute is pulling out the reserve:
    vFreestream = vRocket - x[1]
    pDyn = 0.5 * density * (vFreestream**2)
    
    drag = pDyn * (projectile.Cd * projectile.A)
    acc = (drag - (g * projectile.mass)) / projectile.mass
    return acc
    