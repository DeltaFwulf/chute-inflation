# -*- coding: utf-8 -*-
"""
Created on Sat Jun 11 10:53:27 2022

@author: Tyrtl
"""

import numpy as np
import fcns as fcns

# This module contains the instructions for creating different objects from these classes:

class Parachute:
    # Data in the CSV will be a mix of imperial and metric units, deal with them here:
    def __init__(self, chute_data):
        self.name = chute_data[0]
        self.supplier = chute_data[1]
        self.price = float(chute_data[2])
        self.geometry = chute_data[3]
        self.mass = fcns.oz2kg(float(chute_data[4]))
        self.dPack = fcns.in2m(float(chute_data[5]))
        self.lPack = fcns.in2m(float(chute_data[6]))
        self.CdP = float(chute_data[7])
        self.Cx = float(chute_data[8])
        self.Ap = fcns.sqft2sqm(float(chute_data[9]))
        self.Cd0 = float(chute_data[10])
        self.S0 = fcns.sqft2sqm(float(chute_data[11]))
        self.D0 = fcns.in2m(float(chute_data[12]))
        self.nFill = float(chute_data[13])
        self.maxLoad = fcns.lbf2N(float(chute_data[14]))
            
# Rocket data:
class Rocket:
    def __init__(self, mass, length, diameter):
        self.mass = mass # kg
        self.length = length # m
        self.diameter = diameter # m
        # For now, assume a cylindrical mortar:
        self.mortarDiameter = fcns.in2m(3)
        self.mortarLength = 0.125
        self.boatTailLength = 0.8 * self.diameter
            
class Nosecone:
    def __init__(self, mass, Cd, rocket):
        self.mass = mass # kg
        self.Cd = Cd
        self.A = np.pi * (rocket.diameter / 2) ** 2 # m^2
        
class Riser:
    def __init__(self, rocket, c):
        self.length = 2 * rocket.length
        self.k = 20814 / (2 * rocket.length)
        self.c = c
        
class Flight:
    # All altitudes converted to their ASL values from inputs in AGL
    def __init__(self):
        self.reactionTime = 2
        self.zGround = 300
        self.apogee = 5600 + self.zGround
        self.vApogee = 48
        self.zMainDeploy = 450 + self.zGround
        self.vMainDeploy = 23 # Initially this will be the dummy value; will be changed within the selection function later on
        # Assume freefall in vacuum over reaction time period:
        self.vReserve = np.sqrt((self.vApogee**2) + ((9.81 * self.reactionTime)**2))
        # Velocities:
        self.vMaxMain = 7
        self.vMaxDrogue = 46
        self.vMinDrogue = 23
        self.vMaxReserve = self.vMinDrogue * 0.7
        self.maxAccel = 15 * 9.81
        
        # Launch site temperature, C:
        self.lsTemp = 30#13.05
        
        # Air Properties throughout flight:
        self.densityApogee = fcns.isa(self.apogee, self.lsTemp, self.zGround)[2]
        self.densityMainDeploy = fcns.isa(self.zMainDeploy, self.lsTemp, self.zGround)[2]
        self.densityTouchDown = fcns.isa(self.zGround, self.lsTemp, self.zGround)[2]
        
        
# Define parameters of deployment bag:
class Bag:
    def __init__(self, parachute, mode):
        if(mode == "reserve"):
            self.Cd = 3
        else:
            self.Cd = 1.15
        
        self.A = np.pi * (parachute.dPack / 2)**2
        self.mass = parachute.mass
        self.diameter = parachute.dPack
        self.length = parachute.lPack
        
        
class Config:
    def __init__(self, bestDrogue, bestMain, bestReserve, rocket, flight):
        
        # Parachutes:
        self.bestDrogue = bestDrogue
        self.bestMain = bestMain
        self.bestReserve = bestReserve
        
        # Peak load during flight:
        self.peakLoad = np.max([bestDrogue.openShock, bestDrogue.snatchLoad, bestMain.snatchLoad, bestMain.openShock, bestReserve.snatchLoad1, bestReserve.snatchLoad2, bestReserve.openShock1, bestReserve.openShock2])
        
        # Peak aerodynamic acceleration:
        self.maxAccel = self.peakLoad / rocket.mass
        
        # Touchdown Velocities:
        self.vTouchDown0 = bestMain.vTouchDown
        self.vTouchDown1 = fcns.vSteady(rocket.mass, bestReserve.Cd0, bestReserve.S0, flight.densityTouchDown)
        self.vTouchDown2 = np.sqrt(2 * rocket.mass * 9.81 / (flight.densityTouchDown * ((bestDrogue.Cd0 * bestDrogue.S0) + (bestReserve.Cd0 * bestReserve.S0))))
        
        # Max kinetic energy at touchdown:
        self.KE0 = fcns.ke(rocket.mass, self.vTouchDown0)
        self.KE1 = fcns.ke(rocket.mass, self.vTouchDown1)
        self.KE2 = fcns.ke(rocket.mass, self.vTouchDown2)
        
        # Ideal (lowest) touchdown acceleration under boattail, if all energy dissipated through boattail:
        self.btAccel0 = fcns.btAcc(self.vTouchDown0, rocket.boatTailLength)
        self.btAccel1 = fcns.btAcc(self.vTouchDown1, rocket.boatTailLength)
        self.btAccel2 = fcns.btAcc(self.vTouchDown2, rocket.boatTailLength)
        
        # Cost of configuration:
        self.cost = bestDrogue.price + bestMain.price + bestReserve.price
        
        # Critical Velocities throughout flight:
        self.vDrogueDeploy = bestDrogue.vDeploy
        self.vMainDeploy = bestMain.vDeploy
        self.vReserveDeploy1 = bestReserve.vDeploy1
        self.vReserveDeploy2 = bestReserve.vDeploy2