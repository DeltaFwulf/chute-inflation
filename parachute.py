from math import pi

# TODO: have import of specific cx, cd data from parachute name and a reference file,
#       with geometric values remaining constant depending on type

# TODO: investigate altitude effects of inflation time (constant volume method)


class Parachute():

    def __init__(self, d0:float, porosity:float=0):

        self.type = 'generic'
        self.infExp = 0
        self.d0 = d0
        self.porosity = porosity
        self.cd = 0
        self.Aref = pi * self.d0**2 / 4
        self.cx = 0 # wind tunnel value


    # NOTE: there are two methods specified for tFill in DTIC ADA247666. The equations chosen were those obtained via field testing.
    def tFill(self, v0:float) -> float:

        return 0

           

class FlatCircular(Parachute):

    def __init__(self, d0:float, porosity:float=0):

        self.type = 'flat circular'
        self.infExp = 6
        self.d0 = d0
        self.porosity = porosity
        self.cd = 0.8 # 0.75 to 0.8
        self.Aref = pi * self.d0**2 / 4
        self.cx = 1.7

    def tFill(self, v0:float) -> float:
        """Uses the constant volume relationship for parachutes"""

        if self.porosity < 0.25:
            t0 = 2.5 * self.d0/(v0**0.85)

        else:
            t0 = 4 * self.d0 / (v0**0.85)

        return t0


        
class Ribbon(Parachute):

    def __init__(self, d0:float, porosity:float=0):

        self.type = 'ribbon'
        self.infExp = 1
        self.d0 = d0
        self.porosity = porosity
        self.cd = 0.55
        self.Aref = pi * self.d0**2 / 4
        self.cx = 1.3 # 1.05 to 1.3


    def tFill(self, v0:float) -> float:

        if self.porosity < 0.25:
            tFill = 8 * self.d0 / v0
        else:
            tFill = 0.65 * self.porosity * self.d0 / v0

        return tFill


        
class RingSlot(Parachute):

    def __init__(self, d0:float, porosity:float=0):

        self.type = 'ringslot'
        self.infExp = 1
        self.d0 = d0
        self.porosity = porosity
        self.cd = 0.65 # 0.56 to 0.65
        self.Aref = pi * self.d0**2 / 4


    def tFill(self, v0:float) -> float:
        return 14 * self.d0 / v0
    


class ExtendedSkirt(Parachute):

    def __init__(self, d0:float, porosity:float=0):

        self.type = 'extended skirt'
        self.infExp = 2
        self.d0 = d0
        self.porosity = porosity
        self.cd = 0.9
        self.Aref = pi * self.d0**2 / 4
        self.cx = 1.05

    
    def tFill(self, v0:float) -> float:
        return 12 * self.d0 / v0
    


class Cross(Parachute):

    def __init__(self, d0:float, porosity:float=0):

        self.type = 'cross'
        self.infExp = 3
        self.d0 = d0
        self.porosity = porosity
        self.cd = 0.85 # 0.6 to 0.85
        self.Aref = pi * self.d0**2 / 4
        self.cx = 1.2 # 1.1 to 1.2

    
    def tFill(self, v0:float) -> float:
        return 8.7 * self.cd * self.Aref / (v0**0.9)
    


class BlackCatIris(Parachute):

    # this has been approximated as an annular parachute for the time being, until drop tests can confirm parachute performance

    def __init__(self, d0:float, porosity:float=0):

        self.type = 'blackcat rocketry iris'
        self.infExp = 6
        self.d0 = d0
        self.porosity = porosity
        self.cd = 1.34
        self.Aref = pi * self.d0**2 / 4
        self.cx = 1.4

    
    def tFill(self, v0:float) -> float:
        
        if self.porosity < 0.25:
            t0 = 2.5 * self.d0/(v0**0.85)

        else:
            t0 = 4 * self.d0 / (v0**0.85)

        return t0