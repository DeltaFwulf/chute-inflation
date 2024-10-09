from math import sqrt, pi, sin, cos, exp
import numpy as np
import matplotlib.pyplot as plt
from parachute import Parachute


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



def inflatePflanzLudtke():
    """non-dimensional inflation of a parachute"""

    parachute = Parachute('flat-circular',
                          D0=3.5,
                          Cd=1.9
                          )

    # flight parameters:
    vStretch = 10 # estimate using line stretch simulation
    density = 0.42
    ang0 = 0 # 0 is horizontal, -90 is straight down
    mSystem = 6

    # Estimate inflation time with Knacke's method (constant inflation volume):
    tFill = parachute.tFill(vStretch, density)
    CdS0 = parachute.Cd * parachute.S0
        
    # Parameters for ODEs
    nFill = tFill * vStretch / parachute.D0 # non-dimensional time
    vTerm = sqrt(2 * mSystem * 9.81 / (density * CdS0))

    print(f"terminal velocity = {vTerm}")

    A = vTerm**2 / (9.81 * parachute.D0 * nFill)
    B = A * (vStretch / vTerm)**2 # a combined parameter

    params = {'A':A, 'B':B, 'j':parachute.j}

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

        #print(f"t/tf: {nt[i]}, nv: {nv[i]}, ang{ang[i]}")

    # convert to real forces and plot parameters over inflation period (v *= vs, t *= tFill)
    v = nv * vStretch
    t = nt * tFill

    plt.figure(10)
    plt.plot(t, v, '-b')
    plt.xlabel("Time, s")
    plt.ylabel("Velocity, m/s")

    plt.figure(20)
    plt.plot(t, ang*180/pi, '-k')
    plt.xlabel("Time, s")
    plt.ylabel("Angle, degree")

    # estimate peak loading using base Pflanz-Ludtke method (see if we underpredict)
    Ckl = np.max(x)
    fMaxBasic = Ckl * 0.5 * density * vStretch**2 * CdS0

    # estimate peak loading using extended Pflanz-Ludtke method (see if this improves things)
    C1 = sqrt(parachute.j) * (vTerm/vStretch)**2 * exp(-B)
    C2 = sqrt(parachute.j) * (vTerm/vStretch)**2 * (1 - exp(-B)) * sin(-ang0) * exp(-(A/6)*parachute.j**0.25)

    Ck_ext = Ckl + C1 + C2
    fMaxExt = Ck_ext * 0.5 * density * vStretch**2 * CdS0

    print(f"peak acceleration (extended) = {fMaxExt/(mSystem * 9.81)} g")
    print(f"peak load (basic model): {fMaxBasic} N")
    print(f"peak load (extended model): {fMaxExt} N")

    plt.show()

inflatePflanzLudtke()