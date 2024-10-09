from math import pi



class Parachute():

    inflationExponents = {'flat-circular':6, 
                          'flat-circular-porous':6,
                          'ribbon':1,
                          'ring-slot':1,
                          'ext-skirt':2,
                          'cross':3
                          }
    

    # NOTE: there are two methods specified for tFill in DTIC ADA247666. The equations chosen were those obtained via field testing.
    def tFill(self, V0:float, density:float) -> float:

        match self.geometry:

            case 'flat-circular':
                t0 = 2.5 * self.D0/(V0**0.85)
            
            case 'flat-circular-porous':
                t0 = 4 * self.D0 / (V0**0.85)
            
            case 'ribbon':

                if self.porosity is None:
                    t0 = 8 * self.D0 / V0
                else:
                    t0 = 0.65 * self.porosity * self.D0 / V0
            
            case 'ringslot':
                t0 = 14 * self.D0 / V0
            
            case 'ext-skirt':
                t0 = 12 * self.D0 / V0
            
            case 'cross':
                t0 = 8.7 * self.Cd * self.S0 / (V0**0.9)

        return t0 * (1.225 / density) # this accounts for altitude effects


    def __init__(self, geometry:str, D0:float, Cd:float, porosity:float=None):

        self.geometry = geometry

        self.D0 = D0
        self.S0 = pi * D0**2 / 4
        self.Cd = Cd
        self.porosity = porosity

        self.j = Parachute.inflationExponents[geometry]
        
