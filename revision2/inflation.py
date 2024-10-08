from math import sqrt, pi, sin, cos
import numpy as np
import matplotlib.pyplot as plt


"""This script is used to run a 2-DOF simulation of nondimensionalised inflation of a parachute using an extended Pflanz-Ludtke method.


An annoying thing about this system of equations is that we cannot use the RK methods without a modification:
 - we can get the mid values of each variable by taking the other's initial value up until the half timestep, then iterating the solutions of the next timestep values until they converge. 
   however, it is unclear if modified RK4 is actually any better than simply running Euler with a fine timestep. If so, I would prefer to use it as Euler is stinky poopoo and I don't like it

- turns out, I have invented spicy rk4, but to validate it we must:
    - run this simulation with Euler and a TINY timestep (convergence study)
    - see which of the midpoint values better approximates the converged euler case
"""



def dVelRatio(nv, nt, ang, params):
    """Calculates the derivative of the velocity ratio at a given non-dimensional time t/tf"""

    A = params['A']
    B = params['B']
    j = params['j']

    return -(sin(ang) / B) - (nt**j * nv**2 / A)



def dAngle(ang, nv, params):
    """Calculates the angular derivative at a given non-dimensional time t/tf"""
    
    return  -(1/params['B']) * cos(ang) / nv



def rk4mutated(nvFcn, angFcn, nv0:float, ang0:float, t0:float, h:float, params:dict):
    """Computes the new angle and velocity ratio using a modified RK4 solver"""

    k1 = nvFcn(nv0, t0, ang0, params)
    m1 = angFcn(ang0, nv0, params)

    k2 = nvFcn(nv0 + k1*h/2, t0 + h/2, ang0 + m1*h/2, params)
    m2 = angFcn(ang0 + m1*h/2, nv0 + k1*h/2, params)

    # refine midpoint value of nv, ang: I am currently using the mean of the two slopes, but it may be better to only use the <slope>2 instead
    k3 = nvFcn(nv0 + k2*h/2, t0+h/2, ang0 + m2*h/2, params) 
    m3 = angFcn(ang0 + m2*h/2, nv0 + k2*h/2, params)

    k4 = nvFcn(nv0 + k3*h, t0+h, ang0 + m3*h, params)
    m4 = angFcn(ang0 + k3*h, nv0 + k3*h, params)

    nv = nv0 + h*(k1 + 2*k2 + 2*k3 + k4)/6
    ang = ang0 + h*(m1 + 2*m2 + 2*m3 + m4)/6

    return nv, ang



def inflation():
    """non-dimensional inflation of a parachute"""

    # flight parameters:
    vStretch = 30 # estimate using line stretch simulation
    density = 0.06373  # XXX: use altitude to get an approximation here
    tFill = 0.5 # XXX: obtain from drop testing
    g = 9.81
    ang0 = 0 # initial flight path angle 

    mSystem = 6

    # Parachute-specific data
    D0 = 3.5
    S0 = pi * D0**2 / 4
    Cd0 = 1.9
    CdS0 = Cd0 * S0
    j = 6
    
    # Parameters for ODEs
    nFill = tFill * vStretch / D0 # non-dimensional time
    vTerm = sqrt(2 * mSystem * g / (density * CdS0))

    print(f"terminal velocity = {vTerm}")

    A = vTerm**2 / (g * D0 * nFill)
    B = A * (vStretch / vTerm)**2 # a combined parameter

    params = {'A':A, 'B':B, 'j':j}

    dnt = 1e-3
    ntFinal = 1 # FIXME: this should really be Cx^1/j
    nt = np.arange(0, ntFinal, dnt)
    
    nv = np.zeros((nt.size), float) # this is v/vs
    ang = np.zeros((nt.size), float)
    x = np.zeros((nt.size), float) # the maximum of this array is Ck

    nv[0] = 1
    ang[0] = ang0

    # calculate the non-dimensional inflation behaviour
    for i in range(1, nt.size):

        nv[i], ang[i] = rk4mutated(dVelRatio, dAngle, nv[i-1], ang[i-1], nt[i-1], dnt, params)
        x[i] = nv[i]**2 * nt[i]**params['j']

        print(f"t/tf: {nt[i]}, nv: {nv[i]}, ang{ang[i]}")

    plt.figure(0)
    plt.plot(nt, nv, '-k')

    plt.xlabel("t/tf")
    plt.ylabel("v/vs")
       
    # convert to real forces and plot parameters over inflation period (v *= vs, t *= tFill)
    v = nv * vStretch
    t = nt * tFill

    plt.figure(20)
    plt.plot(t, v, '-b')
    plt.xlabel("time, s")
    plt.ylabel("velocity, m/s")

    plt.show()

    # estimate peak loading using base Pflanz-Ludtke method (see if we underpredict)


    # estimate peak loading using extended Pflanz-Ludtke method (see if this improves things)



inflation()