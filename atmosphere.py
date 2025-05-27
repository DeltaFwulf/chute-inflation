from math import exp

def standardAtmosphere(alt:float, dT:float=0):

    if alt < 11000:

        T = 288.19 - 0.00649 * alt + dT
        P = 101290 * (T / 288.08)**5.256

    elif alt < 25000:

        T = 216.69 + dT
        P = 22650 * exp(1.73 - 0.000157*alt)

    else:

        T = 141.94 + 0.00299*alt + dT
        P = 2488 * (T/216.6)**-11.388

    rho = P / (286.9 * T)

    return T, P, rho