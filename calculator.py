from math import sqrt

from openShock import openShockPL
from parachute import *
from atmosphere import standardAtmosphere


def chuteCalc():

    zGround = 0

    parachute = BlackCatIris(d0=0.6096)

    mParachute = 0.072
    mPayload = 0.6
    
    vStretch = 20
    zStretch = 2000

    rhoInflation = standardAtmosphere(zStretch + zGround)[2]
    rhoTouchdown = standardAtmosphere(zGround)[2]

    angInit = -90 * pi / 180
    extended = True

    mSystem = mParachute + mPayload
    vTerm = sqrt(2 * mSystem * 9.80665 / (rhoTouchdown * (parachute.cd * parachute.Aref)))
    impactEnergy = 0.5 * mSystem * vTerm**2
    
    peakLoad, tFill, tm, tMax, A, vRatio = openShockPL(parachute, mSystem, vStretch, rhoInflation, angInit, extended)

    # TODO: right format this to line up at the colon
    print('')
    print("SETUP INFORMATION ==================================")
    print(f"parachute type:             {parachute.type}")
    print(f"parachute nominal diameter: {parachute.d0} m")
    print(f"parachute drag coefficient: {parachute.cd}")
    print(f"System mass:                {'%.3f' % mSystem} kg")
    print(f"Ground altitude (ASL):      {zGround} m")
    print(f"Deployment altitude (AGL):  {zStretch} m")
    print(f"Air density (inflation):    {'%.3f' % rhoInflation} kg/m^3")
    print(f"Air density (ground):       {'%.3f' % rhoTouchdown} kg/m^3")
    print('')

    print("OPEN SHOCK RESULTS =================================")
    print(f"Velocity at line stretch:   {'%.3f' % vStretch} m/s")
    print(f"Open Shock Load:            {'%.3f' % peakLoad} N")
    print(f"Parachute fill time:        {'%.3f' % tFill} s")
    print(f"Peak load time:             {'%.3f' % tm} s")
    print(f"Peak drag-area time:        {'%.3f' % tMax} s")
    print(f"Ballistic parameter:        {'%.3f' % A}")
    print(f"Velocity ratio:             {'%.3f' % vRatio}")
    print('')

    print("POST-INFLATION =====================================")
    print(f"Touchdown velocity:         {'%.3f' % vTerm} m/s")
    print(f"Impact energy:              {'%.3f' % impactEnergy} J\n")



chuteCalc()