# -*- coding: utf-8 -*-
"""
Created on Tue Jun  7 19:21:36 2022

@author: Tyrtl
"""

# Improvements:
# - Include temperature effects on isa calculations x
# - Replace current inflation time calculations with that of Ludtke's method
# - If possible, find a better number for drag area exponent, j
# - Include multiple variables that can be optimised for i.e. vTouchdown, cost, mass, etc.
# - Allow for simulation of descent in entirety
# - Better organise results: generate an object of the 'best flight', with all relevant parameters. x
# - Include mass of reserve hatch when calculating snatch load.
# - Update vRocket within RK4 algorithm for snatchLoad. x
# - Account for rocket's changing mass as parachutes are deployed.

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
    rocket = cl.Rocket(17.761, 2.5, 0.15)
    nosecone = cl.Nosecone(0.01, 1, rocket)
    riser = cl.Riser(rocket, 0)
    
    # Parachute data:
    #with open(r"C:\Users\Tyrtl\University of Bath\Bath Rocket Team - Recovery - 2022\Computational Model\Python\chute_data\chute_info.csv") as fileName:
    with open(r"C:\Users\Tyrtle13king\University of Bath\Bath Rocket Team - Recovery - 2022\Computational Model\Python\chute_data\chute_info.csv") as fileName:    
        chute_array = np.loadtxt(fileName, delimiter=",", dtype=str)
    
    # Number of timesteps for calculating inflation behaviour
    ntMain = 100 
    ntDrogue = 100
    ntReserve = 100
    
    """Step 1: select main:"""
    
    # Current best peak load value for main under arbitrary drogue:
    minPeakLoadMain = np.inf
    
    for i in range(1, np.shape(chute_array)[0]): # Change the index by one later, is untidy.

        main = cl.Parachute(chute_array[i])
        vTouchDown = fcns.vSteady(rocket.mass, main.CdP, main.Ap, flight.densityTouchDown)

        # Filter for possible main chutes:
        if(vTouchDown < flight.vMaxMain):
                    
            mainBag = cl.Bag(main, "main")
            
            # Placeholder drogue parachute for comparison:
            class drogue:
                pass
            drogue.Cd0 = 1.6 # ul
            drogue.S0 = 1 # m^2
            drogue.mass = 0.1 # kg
            
            # Snatch loading:
            (mainSnatch, mainVS) = sl.snatchLoad(rocket, mainBag, nosecone, riser, flight.vMainDeploy, flight.densityMainDeploy, "main", drogue)
            
            # Open Shock:
            openShockMain = os.openShock("main", rocket, main, flight, mainVS, flight.zMainDeploy, ntMain)[0]
            
            # Is this parachute the best main so far?
            #peakLoadMain[0, i-1] = np.max([mainSnatch, openShockMain])
            peakLoadMain = np.max([mainSnatch, openShockMain])
            
            if(peakLoadMain == np.min([peakLoadMain, minPeakLoadMain])):
                bestMain = main
                bestMain.vTouchDown = vTouchDown
                minPeakLoadMain = peakLoadMain
                         
    # Determine max deployment velocity for this main (implement shooting method... now?): YES! >:)
    tol = 0.01 # m/s
    err = np.ones((1,2))
    vDeployMain = np.array([23, 30])
    peakLoad = np.zeros((1,2))
    
    # Deployment bag information for main chute:
    mainBag = cl.Bag(bestMain, "main")
    
    while(np.abs(err[0,1]) > tol):
        
        for i in range(0, 2):
            # Calculate peak load on bestMain:
            (snatchLoadMain, mainVS) = sl.snatchLoad(rocket, mainBag, nosecone, riser, vDeployMain[i], flight.densityMainDeploy, "main", drogue)
            openShockMain = os.openShock("main", rocket, bestMain, flight, mainVS, flight.zMainDeploy, ntMain)[0]
            peakLoad[0, i] = np.max([snatchLoadMain, openShockMain])
            
        # Calculate error term:
        err = peakLoad - np.min([bestMain.maxLoad, rocket.mass * flight.maxAccel])
        # Produce new guess:
        vNew = fcns.shootingMethod(vDeployMain, err)   
        vDeployMain = np.array([vDeployMain[1], vNew])        
    
    # Set the upper limit for drogue descent rate at main deployment:
    flight.vMaxDrogue = vNew
    
    """Step 2: select drogue:"""
    
    # Relevant Loads:
    drogueSnatch = np.empty((1, np.shape(chute_array)[0] - 1))
    mainSnatch = np.empty((1, np.shape(chute_array)[0] - 1))
    drogueInfLoad = np.empty((np.shape(chute_array)[0] - 1, ntDrogue))
    mainInfLoad = np.empty((np.shape(chute_array)[0] - 1, ntMain))
    drogueOpenShock = np.empty((1, np.shape(chute_array)[0] - 1))
    mainOpenShock = np.empty((1, np.shape(chute_array)[0] - 1))
    peakLoad = np.zeros((1, np.shape(chute_array)[0] - 1))
    
    minPeakLoad = np.inf
    
    for i in range(0, np.shape(chute_array)[0] - 1):
        
        drogue = cl.Parachute(chute_array[i+1])
        vDrogue = fcns.vSteady(rocket.mass, drogue.Cd0, drogue.S0, flight.densityMainDeploy)
        
        # Filter eligible drogues:
        if(vDrogue >= flight.vMinDrogue and vDrogue <= flight.vMaxDrogue):
            
            drogueBag = cl.Bag(drogue, "drogue")
            
            # Drogue Snatch load:
            (drogueSnatch[0, i], drogueVS) = sl.snatchLoad(rocket, drogueBag, nosecone, riser, flight.vApogee, flight.densityApogee, "drogue", drogue)
            # Drogue Open Shock Load:
            (drogueOpenShock[0, i], drogueInfLoad[i, :], drogueTFill) = os.openShock("drogue", rocket, drogue, flight, drogueVS, flight.apogee, ntDrogue)
            # Main Snatch Load:
            (mainSnatch[0, i], mainVS) = sl.snatchLoad(rocket, mainBag, nosecone, riser, vDrogue, flight.densityMainDeploy, "main", bestMain)
            # Main Open Shock Load:
            (mainOpenShock[0, i], mainInfLoad[i, :], mainTFill) = os.openShock("main", rocket, bestMain, flight, mainVS, flight.zMainDeploy, ntMain)
            
    
            # Is this the best drogue so far? (current metric is loading, though will become a switch in the future):
            peakLoad = np.max([drogueSnatch[0, i], drogueOpenShock[0, i], mainSnatch[0, i], mainOpenShock[0, i]])
            
            # Will the drogue survive this descent?
            if(np.max([drogueSnatch[0, i], drogueOpenShock[0, i]]) > drogue.maxLoad):
                peakLoad = np.inf
            
            # Is this the best drogue so far?
            if(peakLoad < minPeakLoad):
                
                # Update Selection:
                bestDrogue = drogue
                
                # Loads:
                bestDrogue.snatchLoad = drogueSnatch[0, i]
                bestDrogue.openShock = drogueOpenShock[0, i]
                bestMain.snatchLoad = mainSnatch[0, i]
                bestMain.openShock = mainOpenShock[0, i]
                bestDrogue.infLoad = drogueInfLoad[i, :]
                bestMain.infLoad = mainInfLoad[i, :]
                
                # Inflation times:
                bestDrogue.tFill = drogueTFill
                bestMain.tFill = mainTFill
                
                # Deployment velocities:
                bestDrogue.vDeploy = drogueVS
                bestMain.vDeploy = mainVS
                
                minPeakLoad = peakLoad
                
 
    """Step 3: select reserve"""
    
    # Best reserve must pack into a mortar and not exceed the maximum acceleration laid out in the rocket object
    drogueBag = cl.Bag(bestDrogue, "drogue")
    
    reserveSnatch1 = np.zeros((1, np.shape(chute_array)[0] - 1))
    reserveOpenShock1 = np.zeros((1, np.shape(chute_array)[0] - 1))
    reserveSnatch2 = np.zeros((1, np.shape(chute_array)[0] - 1))
    reserveOpenShock2 = np.zeros((1, np.shape(chute_array)[0] - 1))
    reserveInfLoad1 = np.zeros((np.shape(chute_array)[0] - 1, ntReserve))
    reserveInfLoad2 = np.zeros((np.shape(chute_array)[0] - 1, ntMain))
    
    minPeakLoadReserve = np.inf
    
    for i in range(0, np.shape(chute_array)[0] - 1):
        
        reserve = cl.Parachute(chute_array[i+1])
        reserveBag = cl.Bag(reserve, "reserve")
        
        # Calculate landing velocity under scenario 1:
        reserveVTouchDown = fcns.vSteady(rocket.mass, reserve.Cd0, reserve.S0, flight.densityTouchDown)
        
        # Does this parachute fit into the mortar?, does it also provide sufficiently low landing velocity:
        if(reserveBag.diameter <= rocket.mortarDiameter and reserveBag.length <= rocket.mortarLength and reserveVTouchDown < flight.vMaxReserve):
            
            """Scenario 1: The drogue fails to deploy; the reserve must instead deploy, acting as the sole decelerator from apogee to touchdown."""
            # Snatch Load, Vs:
            (reserveSnatch1[0, i], reserveVS1) = sl.snatchLoad(rocket, reserveBag, nosecone, riser, flight.vReserve, flight.densityApogee, "reserve1", bestDrogue)
            # Reserve Open Shock:
            (reserveOpenShock1[0, i], reserveInfLoad1[i, :], reserve1TFill) = os.openShock("reserve1", rocket, reserve, flight, reserveVS1, flight.apogee, ntReserve)
            
            """Scenario 2: The drogue deploys successfully, however the main fails to release. The reserve is released to assist the drogue in its deceleration, and touchdown occurs under both parachutes."""
            # Snatch Load, Vs:
            (reserveSnatch2[0, i], reserveVS2) = sl.snatchLoad(rocket, reserveBag, nosecone, riser, flight.vMainDeploy, flight.densityMainDeploy, "reserve2", bestDrogue)
            # Reserve Open Shock:
            (reserveOpenShock2[0, i], reserveInfLoad2[i, :], reserve2TFill) = os.openShock("reserve2", rocket, reserve, flight, reserveVS2, flight.zMainDeploy, ntReserve)
            
            # Peak loads:
            reservePeakLoad = np.max([reserveSnatch1[0, i], reserveOpenShock1[0, i], reserveSnatch2[0, i], reserveOpenShock2[0, i], bestDrogue.snatchLoad, bestDrogue.openShock])
            
            # Does this reserve survive?
            if(reservePeakLoad > reserve.maxLoad):
                reservePeakLoad = np.inf
                
            # Is this reserve the best so far?
            if(reservePeakLoad < minPeakLoadReserve):
                bestReserve = reserve
                bestReserve.snatchLoad1 = reserveSnatch1[0, i]
                bestReserve.snatchLoad2 = reserveSnatch2[0, i]
                bestReserve.openShock1 = reserveOpenShock1[0, i]
                bestReserve.openShock2 = reserveOpenShock2[0, i]
                
                bestReserve.infLoad1 = reserveInfLoad1[i,:]
                bestReserve.infLoad2 = reserveInfLoad2[i,:]
                
                bestReserve.tFill1 = reserve1TFill
                bestReserve.tFill2 = reserve2TFill
                
                bestReserve.vDeploy1 = reserveVS1
                bestReserve.vDeploy2 = reserveVS2
                
                minPeakLoadReserve = reservePeakLoad
        
        
    # Best configuration: (all relevant information will be stored here)
    bestConfig = cl.Config(bestDrogue, bestMain, bestReserve, rocket, flight)

    """ Step 4: findings"""
    
    # Plot the dynamic loads on both main and drogue chutes during inflation:
    tInfDrogue = np.linspace(0, bestDrogue.tFill, ntDrogue)
    tInfMain = np.linspace(0, bestMain.tFill, ntMain)
    tInfReserve1 = np.linspace(0, bestReserve.tFill1, ntReserve)
    tInfReserve2 = np.linspace(0, bestReserve.tFill2, ntReserve)
    
    # Plot the drag loading during inflation:
    plt.style.use('classic')       
    fig,ax = plt.subplots()
    drogueForce, = ax.plot(tInfDrogue, bestDrogue.infLoad, linewidth=2.0, label = "Drogue")
    mainForce, = ax.plot(tInfMain, bestMain.infLoad, linewidth=2.0, label = "Main")
    reserveForce1, = ax.plot(tInfReserve1, bestReserve.infLoad1, linewidth=2.0, label = "Reserve 1")
    reserveForce2, = ax.plot(tInfReserve2, bestReserve.infLoad2, linewidth=2.0, label = "Reserve 2")
    ax.legend(handles = [drogueForce, mainForce, reserveForce1, reserveForce2], loc = 2)
    plt.xlabel("Time, s")
    plt.ylabel("Drag, N")
    plt.show       
    
    # Show selection:
    print("Best Main: ", bestConfig.bestMain.name)
    print("Best Drogue: ", bestConfig.bestDrogue.name)
    print("Best Reserve: ", bestConfig.bestReserve.name)
    
    # Peak Loads and breakdown:
    print("Drogue Snatch: ", bestDrogue.snatchLoad, " N")
    print("Drogue OS: ", bestDrogue.openShock, " N")
    print("Reserve Snatch 1: ", bestReserve.snatchLoad1, " N")
    print("Reserve OS1: ", bestReserve.openShock1, " N")
    print("Main Snatch: ", bestMain.snatchLoad, " N")
    print("Main OS: ", bestMain.openShock, " N")
        
    print("Peak Load: ", "{:.2f}".format(bestConfig.peakLoad), " N")
    print("Max Aerodynamic Acceleration: ", "{:.2f}".format(bestConfig.maxAccel), " m/s^2")
    
    # Touchdown Velocities:
    print("Scenario 0 Touchdown Velocity: ", "{:.2f}".format(bestConfig.vTouchDown0), " m/s")
    print("Scenario 1 Touchdown Velocity: ", "{:.2f}".format(bestConfig.vTouchDown1), " m/s")
    print("Scenario 2 Touchdown Velocity: ", "{:.2f}".format(bestConfig.vTouchDown2), " m/s")
    
    # Kinetic Energies and minimum possible landing accelerations:
    print("Scenario 0 KE: ", "{:.2f}".format(bestConfig.KE0), " J ", "| Ideal Deceleration: ", "{:.2f}".format(bestConfig.btAccel0), " m/s^2")
    print("Scenario 1 KE: ", "{:.2f}".format(bestConfig.KE1), " J ", "| Ideal Deceleration: ", "{:.2f}".format(bestConfig.btAccel1), " m/s^2")
    print("Scenario 2 KE: ", "{:.2f}".format(bestConfig.KE2), " J ", "| Ideal Deceleration: ", "{:.2f}".format(bestConfig.btAccel2), " m/s^2")
    
    # Velocities:
    print("Main deploy velocity: ", "{:.2f}".format(bestConfig.vMainDeploy), " m/s")
    print("Drogue deployment velocity: ", "{:.2f}".format(bestConfig.vDrogueDeploy), " m/s")
    print("Reserve deployment velocity, scenario 1: ", "{:.2f}".format(bestConfig.vReserveDeploy1), " m/s")
    print("Reserve deployment velocity, scenario 2: ", "{:.2f}".format(bestConfig.vReserveDeploy2), " m/s")
    
    # Price of system:
    print("Cost of parachutes: $", "{:.2f}".format(bestConfig.cost))
    
    return bestConfig
    
select_chutes()