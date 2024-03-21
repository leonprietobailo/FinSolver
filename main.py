import numpy as np
import ClassUtils as cu
import matplotlib.pyplot as plt
# Boundary Conditions
Tb = 100;       # ºC
Tinf = 20;      # ºC

# Geometry
L = 100e-3;     # m
H = 1e-3;       # m
W = 5e-3;       # m
N = 100;

# Material properties
k = 385;        # W/m·K
b = 25;         # W/m²·K

dim = cu.Dimension(L,W,H,N);
mat = cu.Material(k, b);
bound = cu.Boundary(Tb, Tinf);

fin = cu.Fin(dim, mat, bound);

adiabaticTemperature = fin.computeAdiabaticTemperature;
convectiveTemperature = fin.computeConvectiveTemperature;

plt.plot(adiabaticTemperature, label="Adiabatic Temperature")
plt.plot(convectiveTemperature, label="Convective Temperature")
plt.xlabel("Distance [m]")
plt.ylabel("Value [ºC]")
plt.title("Adiabatic Vs Convective Temperature")
plt.show();



# def init():


# def main():



# def computePerimeter(w,t):
#     return 2*(w+t);

# def computeArea(w,t):
#     return w*t;

# def computeM(h, k, A, P):
#     return math.sqrt(h*P/k/A);

# def discretizeFin(L, N):
#     return np.linspace(0,L,N);

# def computeConvectiveTemperature(Tinf, Tb,h, k, A, P, L, N):
#     xa = discretizeFin(L, N);
#     m = computeM(h,k,A,P);
#     Ta1 = Tinf + (Tb -Tinf)*((np.cosh(m*(L-xa))) + h*m*k*np.sinh((m*(L -xa))))/(np.cosh(m*L) + h*m*k*np.sinh(m*L))

# def computeConvectiveQ(h, k, A, P, L):
#     m = computeM(h,k,A,P);
#     M = math.sqrt(k*A*h*P);
#     qf1 = M * math.sinh(m*L) + (h*m*k*math.cosh(m*L) + h*m*k*math.sinh(m*L));

# def computeAdiabaticTemperature(Tinf, Tb,h, k, A, P, L, N):
#     xa = discretizeFin(L, N);
#     m = computeM(h,k,A,P);
#     Ta2 = Tinf + (Tb - Tinf)*(np.cosh(m*(L - xa))/np.cosh(m*L))

# def computeAdiabaticQ(h, k, A, P, L):
#     m = computeM(h,k,A,P);
#     M = math.sqrt(k*A*h*P);
#     qf2 = M*(math.tanh(m*L))