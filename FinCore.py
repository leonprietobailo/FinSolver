import math
import numpy as np

class Fin:
    def __init__(self, dim, mat, bound, mL = None):
        self.dim: Dimension = dim
        self.mat: Material = mat
        self.bound: Boundary = bound
        self.mL: float = mL
        self.cache: Cache = Cache()

    @property
    def computeConvectiveTemperature(self):
        s = self
        if(s.cache.getConvectiveTemperature is not None):
            return s.cache.getConvectiveTemperature
        else:
            _, _, m, _, hmk, xa = s.__computeGenerals()
            T = s.bound.Tinf + (s.bound.Tb - s.bound.Tinf)*((np.cosh(m*(s.dim.L-xa)))+ hmk*np.sinh((m*(s.dim.L -xa))))/(np.cosh(m*s.dim.L) + hmk*np.sinh(m*s.dim.L))
            s.cache.setConvectiveTemperature(T)
            return T

    @property
    def computeConvectiveHeat(self):
        s = self
        if(s.cache.getConvectiveHeat is not None):
            return s.cache.getConvectiveHeat
        else:
            _, _, m, M, hmk, _= s.__computeGenerals()
            Q =  M * math.sinh(m*s.dim.L) + (hmk*math.cosh(m*s.dim.L) + hmk*math.sinh(m*s.dim.L))
            s.cache.setConvectiveHeat(Q)
            return Q

    @property
    def computeAdiabaticTemperature(self):
        s = self
        if(s.cache.getAdiabaticTemperature is not None):
            return s.cache.getAdiabaticTemperature
        else:
            _, _, m, _, _, xa = s.__computeGenerals()
            mL = s.mL if s.mL is not None else m * s.dim.L
            T = s.bound.Tinf + (s.bound.Tb - s.bound.Tinf)*(np.cosh(mL - m * xa)/np.cosh(mL))
            s.cache.setAdiabaticTemperature(T)
            return T

    @property
    def computeAdiabaticHeat(self):
        s = self
        if(s.cache.getAdiabaticHeat is not None):
            return s.cache.getAdiabaticHeat
        else:
            _, _, m, M, _, _= s.__computeGenerals()
            mL = s.mL if s.mL is not None else m * s.dim.L
            Q =  M*(math.tanh(mL))
            s.cache.setAdiabaticHeat(Q)
            return Q

    @property
    def computeEfficiency(self):
        s = self
        _, _, m, _, _, _= s.__computeGenerals()
        return 1/m/s.dim.L

    @property
    def computeAdiabaticEffectiveness(self):
        s = self
        _, A, _, _, _, _= s.__computeGenerals()
        return s.computeAdiabaticHeat / s.mat.h / A / (s.bound.Tb - s.bound.Tinf)
    
    @property
    def computeConvectiveEffectiveness(self):
        s = self
        _, A, _, _, _, _= s.__computeGenerals()
        return s.computeConvectiveHeat / s.mat.h / A / (s.bound.Tb - s.bound.Tinf)

    @property
    def computeInfiniteFinTemperature(self):
        s = self
        
        if(s.cache.getInfiniteTemperature is not None):
            return s.cache.getInfiniteTemperature
        else:
            _, _, m, _, _, xa = s.__computeGenerals()
            T = s.bound.Tinf + (s.bound.Tb - s.bound.Tinf)*(np.exp(-m*xa))
            s.cache.setInfiniteTemperature(T)
            return T
        
    @property
    def computeInfiniteFinHeat(self):
        s = self
        if(s.cache.getInfiniteHeat is not None):
            return s.cache.getInfiniteHeat
        else:
            _, _, _, M, _, _ = s.__computeGenerals()
            Q = M * (s.bound.Tb - s.bound.Tinf)
            print(s.computeML)
            s.cache.setInfiniteHeat(Q)
            return Q
        
    @property
    def computeML(self):
        s = self
        _, _, m, _, _, _ = s.__computeGenerals()
        return s.dim.L*m

# Private methods

    def __computeGenerals(self):
        s = self
        P =  2 * (s.dim.H + s.dim.W)
        A = s.dim.H * s.dim.W
        m = math.sqrt(s.mat.h*P/s.mat.k/A)
        M = math.sqrt(s.mat.k*A*s.mat.h*P)
        hmk = s.mat.h/m/s.mat.k
        xa = s.dim.discretize
        return P, A, m, M, hmk, xa


# Support classes

class Dimension:
    def __init__(self, Length, Width, Height, Divisions):
        self._Length = Length
        self._Width = Width
        self._Height = Height
        self._Divisions = Divisions

    @property
    def L(self):
        return self._Length

    @property
    def W(self):
        return self._Width

    @property
    def H(self):
        return self._Height

    @property
    def N(self):
        return self._Divisions

    @property
    def discretize(self):
        return np.linspace(0,self._Length,self._Divisions)

class Material:
    def __init__(self, k, h):
        self._k = k
        self._h = h

    @property
    def k(self):
        return self._k

    @property
    def h(self):
        return self._h

class Boundary:
    def __init__(self, Tb, Tinf):
        self._Tb = Tb
        self._Tinf = Tinf

    @property
    def Tb(self):
        return self._Tb

    @property
    def Tinf(self):
        return self._Tinf
    
class Cache:
    def __init__(self):
        self._convectiveTemperature = None
        self._adiabaticTemperature = None
        self._convectiveHeat = None
        self._adiabaticHeat = None
        self._infiniteTemperature = None
        self._infiniteHeat = None

    @property
    def getConvectiveTemperature(self):
        return self._convectiveTemperature
    
    @property
    def getAdiabaticTemperature(self):
        return self._adiabaticTemperature
    
    @property
    def getConvectiveHeat(self):
        return self._convectiveHeat
    
    @property
    def getAdiabaticHeat(self):
        return self._adiabaticHeat
    
    @property
    def getInfiniteTemperature(self):
        return self._infiniteTemperature
    
    @property
    def getInfiniteHeat(self):
        return self._infiniteHeat
    
    def setConvectiveTemperature(self, convectiveTemperature):
        self._convectiveTemperature = convectiveTemperature
        
    def setConvectiveHeat(self, convectiveHeat):
        self._convectiveHeat = convectiveHeat
        
    def setAdiabaticTemperature(self, adiabaticTemperature):
        self._adiabaticTemperature = adiabaticTemperature
        
    def setAdiabaticHeat(self, adiabaticHeat):
        self._adiabaticHeat = adiabaticHeat

    def setInfiniteTemperature(self, infiniteTemperature):
        self._infiniteTemperature = infiniteTemperature

    def setInfiniteHeat(self, infiniteHeat):
        self._infiniteHeat = infiniteHeat