from math import sqrt, pi, sin, cos, exp
import numpy as np
import matplotlib.pyplot as plt
from parachute import *


"""This is a research script to help myself understand the explanation of Extended Pflanz-Ludtke Method, by comparing numerically solved trajectories to EPL relations (uncorrected and corrected)"""

def dimensionalRK4(v0:float, ang0:float, t0:float, dt:float, etc:dict):

    pass


def nonDimensionalRK4(nv0:float, ang0:float, n0:float, dn:float, A, B, j):
    """Computes the new angle and velocity ratio using a modified RK4 solver
    
    An annoying thing about this system of equations is that we cannot use the RK methods without a modification:
    - we can get the mid values of each variable by taking the other's initial value up until the half timestep, then iterating the solutions of the next timestep values until they converge. 
    however, it is unclear if modified RK4 is actually any better than simply running Euler with a fine timestep. If so, I would prefer to use it as Euler is stinky poopoo and I don't like it
 
    The RK4 Algorithm has been altered here to 'stair-step' the intermediate values in order to couple them as well as I can over a given timestep. 
    # TODO: find notes and explain what I did here lol
    """

    def dnv_dn(nv, nt, ang, A, B, j):
        """Calculates the derivative of the velocity ratio at a given non-dimensional time t/tf"""
        return -(sin(ang) / B) - (nt**j * nv**2 / A)


    def dAng_dn(ang, nv, B):
        """Calculates the angular derivative at a given non-dimensional time t/tf"""
        return  -(1/B) * cos(ang) / nv


    k1 = dnv_dn(nv0, n0, ang0, A, B, j)
    m1 = dAng_dn(ang0, nv0, B)

    k2 = dnv_dn(nv0 + k1*dn/2, n0 + dn/2, ang0 + m1*dn/2, A, B, j)
    m2 = dAng_dn(ang0 + m1*dn/2, nv0 + k1*dn/2, B)

    # refine midpoint value of nv, ang: I am currently using the mean of the two slopes, but it may be better to only use the <slope>2 instead
    k3 = dnv_dn(nv0 + k2*dn/2, n0+dn/2, ang0 + m2*dn/2, A, B, j) 
    m3 = dAng_dn(ang0 + m2*dn/2, nv0 + k2*dn/2, B)

    k4 = dnv_dn(nv0 + k3*dn, n0+dn, ang0 + m3*dn, A, B, j)
    m4 = dAng_dn(ang0 + k3*dn, nv0 + k3*dn, B)

    nv = nv0 + dn*(k1 + 2*k2 + 2*k3 + k4)/6
    ang = ang0 + dn*(m1 + 2*m2 + 2*m3 + m4)/6

    return nv, ang



def maxLoad(parachute:Parachute, masses:dict, Ain:float, vRatio:float, airDensity:float, angInit:float):

    """For a given range of ballistic parameters and cx, outputs expected maximum values of load between a 2dof simulation, pl and extended pl
    
    
    outputs:
    - A, ck by method
    """
    
    mSystem = masses['parachute'] + masses['payload']
    cds0 = parachute.cd * parachute.Aref # steady drag area
    vTerm = sqrt(2 * mSystem * 9.81 / (airDensity * cds0))
    vStretch = vRatio * vTerm # enforced fixed value
    A = Ain
    B = A * (vStretch / vTerm)**2 # a combined parameter
   
    # Pflanz-Ludtke Variables #################################################################################################

    tFill = parachute.tFill(vStretch)
    nFill = tFill * vStretch / parachute.d0 # non-dimensional time
   
    #A = vTerm**2 / (9.81 * parachute.d0 * nFill)
    
    # 2 DOF numerical solution ################################################################################################

    dTau = 1e-3
    tau = np.arange(0, parachute.cx**(1/parachute.infExp), dTau) # tau = t / tfill
    
    nv = np.zeros((tau.size), float) # non-dimensional velocity ratio, (v/vStretch)
    ang = np.zeros((tau.size), float)
    x = np.zeros((tau.size), float)

    nv[0] = 1
    ang[0] = angInit

    for i in range(1, tau.size):
        nv[i], ang[i] = nonDimensionalRK4(nv[i-1], ang[i-1], tau[i-1], dTau, A, B, parachute.infExp)
        x[i] = nv[i]**2 * tau[i]**parachute.infExp

    # v = nv * vStretch
    # t = tau * tFill

    ck_sim = np.max(x)
    
    # Extended Pflanz-Ludtke Method ###########################################################################################
    j = parachute.infExp
    
    tm = ((j * (j + 1) * A) / (j+2))**(1/(j+1)) * tFill # time at which maximum load occurs
    Alx = (j+2) / (j * (j+1)) * parachute.cx**((j+1)/j)

    if A < Alx: # tm < tMax, vehicle decelerates significantly throughout inflation
        ck_pl = ((j+2) / (2 * (j+1)))**2 * ((j*(j+1)*A)/(j+2))**(j/(j+1))
    else: # peak clk occurs at tm = tMax
        ck_pl = (1 + (1 / (A*(j+1))) * parachute.cx**((j+1)/j))**-2 * parachute.cx

    # estimate peak loading using extended Pflanz-Ludtke method (see if this improves things)
    C1 = sqrt(parachute.infExp) * (vTerm/vStretch)**2 * exp(-B)
    C2 = sqrt(parachute.infExp) * (vTerm/vStretch)**2 * (1 - exp(-B)) * sin(-angInit) * exp(-(A/6)*parachute.infExp**0.25)

    ck_ext = ck_pl + C1 + C2

    F_sim = ck_sim * 0.5 * airDensity * vStretch**2 * cds0
    F_pl = ck_pl * 0.5 * airDensity * vStretch**2 * cds0
    F_ext = ck_ext * 0.5 * airDensity * vStretch**2 * cds0

    # compare ck obtained from 2-dof simulation against that obtained by extended pflanz-ludtke
    # print(f"Ck obtained from 2-DOF simulation: {ck_sim}, load: {F_sim} N")
    # print(f"Ck obtained from Pflanz-Ludtke: {ck_pl}, load: {F_pl} N")
    # print(f"Ck obtained from extended Pflanz-Ludtke: {ck_ext}, load: {F_ext} N")

    return ck_sim, ck_pl, ck_ext



def compareMethods():

    # Decelerator Configuration ###############################################################################################
    parachute = FlatCircular(d0=1)
    #parachute = Ribbon(d0=1)
    mParachute = 0.1
    mPayload = 10
    masses = {'parachute':mParachute, 'payload':mPayload}

    # Initial flight state ####################################################################################################
    vRatio = 2.5 # enforce this to validate plots seen in paper
    airDensity = 1.225
    angInit = -90 * pi / 180 # 0 is horizontal, -90 is straight down

    nPoints = 50
    ind = np.linspace(-2, 2, nPoints, float)
    A = np.zeros(nPoints, float)
    ck_sim = np.zeros((A.size), float)
    ck_pl = np.zeros((A.size), float)
    ck_ext = np.zeros((A.size), float)

    for i in range(0, A.size):

        A[i] = 10**ind[i]
        ck_sim[i], ck_pl[i], ck_ext[i] = maxLoad(parachute, masses, A[i], vRatio, airDensity, angInit)

    fig, ax = plt.subplots()
    ax.plot(A, ck_sim, '-k')
    ax.plot(A, ck_pl, 'og')
    ax.plot(A, ck_ext, '*r')
    ax.set_ylabel('Ck')
    ax.set_xlabel('A')
    ax.set_xscale('log')
    ax.legend(['simulation', 'pflanz-ludtke', 'extended'])

    plt.show()



compareMethods()