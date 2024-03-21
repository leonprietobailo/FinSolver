import math
import numpy as np

class Fin:
    def __init__(self, dim, mat, bound):
        self.dim: Dimension = dim;
        self.mat: Material = mat;
        self.bound: Boundary = bound;

    # @property
    # def getDimensions(self):
    #     return self.Dimensions;

    # @property
    # def getMaterial(self):
    #     return self.Material;

    # @property
    # def getBoundary(self):
    #     return self.Boundary;

    @property
    def computeConvectiveTemperature(self):
        s = self;
        _, _, m, _, xa = s.__computeGenerals();
        return s.bound.Tinf + (s.bound.Tb - s.bound.Tinf)*((np.cosh(m*(s.dim.L-xa))) + s.mat.b*m*s.mat.k*np.sinh((m*(s.dim.L - xa))))/(np.cosh(m*s.dim.L) + s.mat.b*m*s.mat.k*np.sinh(m*s.dim.L));

    @property
    def computeConvectiveQ(self):
        s = self;
        _, _, m, M, _= s.__computeGenerals();
        return M * math.sinh(m*s.dim.L) + (s.mat.b*m*s.mat.k*math.cosh(m*s.dim.L) + s.mat.b*m*s.mat.k*math.sinh(m*s.dim.L));

    @property
    def computeAdiabaticTemperature(self):
        s = self;
        _, _, m, _, xa = s.__computeGenerals();
        return s.bound.Tinf + (s.bound.Tb - s.bound.Tinf)*(np.cosh(m*(s.dim.L - xa))/np.cosh(m*s.dim.L));

    @property
    def computeAdiabaticQ(self):
        s = self;
        _, _, m, M, _= s.__computeGenerals();
        return M*(math.tanh(m*s.dim.L))

    def __computeGenerals(self):
        s = self;
        P =  2 * (s.dim.H + s.dim.L);
        A = s.dim.H * s.dim.L;
        m = math.sqrt(s.mat.b*P/s.mat.k/A);
        M = math.sqrt(s.mat.k*A*s.mat.b*P);
        xa = s.dim.discretize;
        return P, A, m, M, xa;

class Dimension:
    def __init__(self, Length, Width, Height, Divisions):
        self._Length = Length;
        self._Width = Width;
        self._Height = Height;
        self._Divisions = Divisions;

    @property
    def L(self):
        return self._Length;

    @property
    def W(self):
        return self._Width;

    @property
    def H(self):
        return self._Height;

    @property
    def N(self):
        return self._Divisions;

    @property
    def discretize(self):
        return np.linspace(0,self._Length,self._Divisions);

class Material:
    def __init__(self, k, b):
        self._k = k;
        self._b = b;

    @property
    def k(self):
        return self._k;

    @property
    def b(self):
        return self._b;

class Boundary:
    def __init__(self, Tb, Tinf):
        self._Tb = Tb;
        self._Tinf = Tinf;

    @property
    def Tb(self):
        return self._Tb;

    @property
    def Tinf(self):
        return self._Tinf;