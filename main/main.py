# -*- coding: utf-8 -*-
"""
Created on Tue Jun  7 19:21:36 2022

@author: Tyrtl
"""

import numpy as np
import fcns as fcns
import snatch_load as sl
import openShock as os
import matplotlib.pyplot as plt
import classes as cl

# This function is the main executed file, calling all sub-functions.
def select_chutes():
    
    # Flight data: (at 83 degrees with windspeed avg 5 m/s)
    flight = cl.Flight()
    rocket = cl.Rocket(26.67, 3.7, 0.15)
    nosecone = cl.Nosecone(0.01, 1, rocket)
    riser = cl.Riser(rocket, 0)
    
    # Parachute data:SSS
    with open(r"C:\Users\Tyrtl\University of Bath\Bath Rocket Team - Recovery - 2022\Computational Model\Python\chute_data\chute_info.csv") as file_name:
        chute_array = np.loadtxt(file_name, delimiter=",", dtype=str)
    
    # Number of timesteps for calculating inflation behaviour
    ntMain = 100 
    ntDrogue = 100
    
    # Declare numpy arrays:
    vTouchDown = np.zeros((1, np.shape(chute_array)[0] - 1))
    peakLoadMain = np.zeros((1,np.shape(chute_array)[0] - 1))
    
    # Air density values:
    densityApogee = fcns.isaDensity(flight.apogee + flight.zGround)
    densityMainDeploy = fcns.isaDensity(flight.zMainDeploy + flight.zGround)
    densityGround = fcns.isaDensity(flight.zGround)     
    
    # From all parachutes, select best main:
    for i in range(1, np.shape(chute_array)[0]): # Change the index by one later, is untidy.

        main = cl.Parachute(chute_array[i])
        vTouchDown[0, i-1] = fcns.vSteady(rocket.mass, main.CdP, main.Ap, densityGround)

        # Filter for possible main chutes:
        if(vTouchDown[0, i-1] < flight.vMaxMain):
                    
            mainBag = cl.Bag(main)
            
            # Placeholder drogue parachute for comparison:
            class drogue:
                pass
            drogue.Cd = 1.6 # ul
            drogue.A = 1 # m^2
            drogue.mass = 0.1 # kg
            
            # Snatch loading:
            (mainSnatch, mainVS) = sl.snatchLoad(rocket, mainBag, nosecone, riser, flight.vMainDeploy, densityMainDeploy, "main", drogue)
            
            # Open Shock:
            openShockMain = os.openShock("main", rocket, main, mainVS, flight.zMainDeploy, ntMain)[0]
            
            # Is this parachute the best main so far?
            peakLoadMain[0, i-1] = np.max([mainSnatch, openShockMain])
            if(peakLoadMain[0, i-1] == np.min(peakLoadMain[np.nonzero(peakLoadMain)])):
                bestMain = main
                bestMain.touchDown = vTouchDown[0, i-1]
    
    # Select best drogue:
    # Rationale: The slower the rocket is going when deploying the main, the lower the main shocks will be, however the larger that the drogue forces will be;
    # Optimum hits a balance point if possible
    # Under current best main, simulate a full descent under all eligible drogues, calculate peak loads:
    
    # Relevant Velocities (for alternative comparisons):
    vDrogue = np.zeros((1, np.shape(chute_array)[0] - 1))
    
    # Relevant Loads:
    drogueSnatch = np.empty((1, np.shape(chute_array)[0] - 1))
    mainSnatch = np.empty((1, np.shape(chute_array)[0] - 1))
    drogueInfLoad = np.empty((np.shape(chute_array)[0] - 1, ntDrogue))
    mainInfLoad = np.empty((np.shape(chute_array)[0] - 1, ntMain))
    drogueOpenShock = np.empty((1, np.shape(chute_array)[0] - 1))
    mainOpenShock = np.empty((1, np.shape(chute_array)[0] - 1))
    peakLoad = np.zeros((1, np.shape(chute_array)[0] - 1))
    
    # Nondimensional times (convert to true inflation times later):
    drogueTn = np.zeros((np.shape(chute_array)[0] - 1, ntDrogue))
    mainTn = np.zeros((np.shape(chute_array)[0] - 1, ntMain))
    
    # Deployment bag information for main chute:
    mainBag = cl.Bag(bestMain)
    
    for i in range(0, np.shape(chute_array)[0] - 1):
        
        drogue = cl.Parachute(chute_array[i+1])
        vDrogue = fcns.vSteady(rocket.mass, drogue.Cd0, drogue.S0, densityMainDeploy)
        print(vDrogue)
        
        # Filter eligible drogues:
        if(vDrogue >= flight.vMinDrogue and vDrogue <= flight.vMaxDrogue):
            
            drogueBag = cl.Bag(drogue)
            
            # Drogue Snatch load:
            (drogueSnatch[0, i], drogueVS) = sl.snatchLoad(rocket, drogueBag, nosecone, riser, flight.vApogee, densityApogee, "drogue", drogue)
            # Drogue Open Shock Load:
            (drogueOpenShock[0, i], drogueInfLoad[i, :], drogueTn[i, :]) = os.openShock("drogue", rocket, drogue, drogueVS, flight.apogee, ntDrogue)
            # Main Snatch Load:
            (mainSnatch[0, i], mainVS) = sl.snatchLoad(rocket, mainBag, nosecone, riser, vDrogue, densityMainDeploy, "main", bestMain)
            # Main Open Shock Load:
            (mainOpenShock[0, i], mainInfLoad[i, :], mainTn[i, :]) = os.openShock("main", rocket, main, mainVS, flight.zMainDeploy, ntMain)
            
    
            # Is this the best drogue so far? (current metric is loading, though will become a switch in the future):
            peakLoad[0, i] = np.max([drogueSnatch[0, i], drogueOpenShock[0, i], mainSnatch[0, i], mainOpenShock[0, i]])
            if(peakLoad[0, i] == np.min(peakLoad[np.nonzero(peakLoad)])):
                bestDrogue = drogue
                bestDrogue.vMainDeploy = vDrogue
    
    
    
    
    
    # Select best reserve:
    # Best reserve must pack into 
                
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    # Useful information:
    print("Best main: ", bestMain.name) 
    #print("Best Drogue: ", bestDrogue.name)
    #print("Open Shock: ", bestMain.openShock, " N")
    #print("Main Snatch Load: ", bestMain.snatchLoad, " N") 
    print("Touchdown Velocity: ", bestMain.touchDown, " m/s")
    #print("Maximum Acceleration: ", bestMain.peakLoad / (9.81 * rocket.mass), " g")
    """
    # Plot the drag loading during inflation:
    plt.style.use('classic')       
    
    fig,ax = plt.subplots()
    ax.plot(tn[0], fOpenMain[0], linewidth=2.0)
    ax.set(xlim=(0, tn.max(1)[0]), xticks = np.linspace(0,tn.max(1)[0], 11), ylim=(0, fOpenMain.max(1)[0]), yticks = np.linspace(0, fOpenMain.max(1)[0], 11))

    plt.show       
    """
select_chutes()