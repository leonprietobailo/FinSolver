import numpy as np
import FinCore as core
import matplotlib.pyplot as plt

# Boundary Conditions
Tb = 100       # ºC
Tinf = 20      # ºC

# Geometry
L = 100e-3     # m
t = 1e-3       # m
w = 5e-3       # m
N = 100

# Material properties
k = 385        # W/m·K
h = 25         # W/m²·K

# Class instantiation
dim = core.Dimension(L,w,t,N)
mat = core.Material(k, h)
bound = core.Boundary(Tb, Tinf)

fin = core.Fin(dim, mat, bound)

# Exercise 1.
print("*************** EXERCISE: ONE ***************")
ex1_adiabaticTemperature = fin.computeAdiabaticTemperature
ex1_convectiveTemperature = fin.computeConvectiveTemperature
ex1_adiabaticQ = fin.computeAdiabaticHeat
ex1_convectiveQ = fin.computeConvectiveHeat
print("Adiabatic Heat Loss: ", ex1_adiabaticQ, "Convective Heat Loss: ", ex1_convectiveQ)

fig = plt.figure()
ax1 = fig.add_subplot(1, 1, 1)
ax1.set_title("Adiabatic and Convective borders.")
ax1.set_xlabel("Distance [m]")
ax1.set_ylabel("Temperature [ºC]")
ax1.plot(np.linspace(0, L, N),ex1_adiabaticTemperature, label="Adiabatic Hypothesis")
ax1.plot(np.linspace(0, L, N),ex1_convectiveTemperature, label="Convective Hypothesis")
ax1.legend()
fig.savefig('output/ex1.svg')

# Exercise 2
print("*************** EXERCISE: TWO ***************")
ex2_efficiency = fin.computeEfficiency
ex2_adiabaticEffectiveness = fin.computeAdiabaticEffectiveness
ex2_convectiveEffectiveness = fin.computeConvectiveEffectiveness
print("Efficiency:", ex2_efficiency,"Adiabatic Effectiveness:", ex2_adiabaticEffectiveness, "Convective Effectivenes:", ex2_convectiveEffectiveness)

# Exercise 3
print("*************** EXERCISE: THREE ***************")
ex3_temperature = fin.computeInfiniteFinTemperature
ex3_heat = fin.computeInfiniteFinHeat

fig = plt.figure()
ax1 = fig.add_subplot(1, 1, 1)
ax1.plot(np.linspace(0, L, N),ex1_adiabaticTemperature, label="Adiabatic Hypothesis")
ax1.plot(np.linspace(0, L, N),ex1_convectiveTemperature, label="Convective Hypothesis")
ax1.plot(np.linspace(0, L, N), ex3_temperature, label="Infinite Hypothesis")
ax1.set_xlabel("Distance [m]")
ax1.set_ylabel("Temperature [ºC]")
ax1.legend()
ax1.set_title("Adiabatic, Convective and Infinite borders.")
fig.savefig('output/ex3.svg')

# Exercise 4
print("*************** EXERCISE: FOUR ***************")
ex4_ml = fin.computeML
print("ML Factor:", ex4_ml)

# Exercise 5
print("*************** EXERCISE: FIVE ***************")

ex5_L1 = 1
ex5_L2 = 0.01

ex5_dim_L1 = core.Dimension(ex5_L1, w, t, N)
ex5_dim_L2 = core.Dimension(ex5_L2, w, t, N)

ex5_fin_L1 = core.Fin(ex5_dim_L1, mat, bound)
ex5_fin_L2 = core.Fin(ex5_dim_L2, mat, bound)

ex5_adiabaticTemperature_L1 = ex5_fin_L1.computeAdiabaticTemperature
ex5_adiabaticTemperature_L2 = ex5_fin_L2.computeAdiabaticTemperature

ex5_convectiveTemperature_L1 = ex5_fin_L1.computeConvectiveTemperature
ex5_convectiveTemperature_L2 = ex5_fin_L2.computeConvectiveTemperature

ex5_temperature_L1 = ex5_fin_L1.computeInfiniteFinTemperature
ex5_temperature_L2 = ex5_fin_L2.computeInfiniteFinTemperature

ex5_adiabaticQ_L1 = ex5_fin_L1.computeAdiabaticHeat
ex5_adiabaticQ_L2 = ex5_fin_L2.computeAdiabaticHeat

ex5_convectiveQ_L1 = ex5_fin_L1.computeConvectiveHeat
ex5_convectiveQ_L2 = ex5_fin_L2.computeConvectiveHeat

ex5_heat_L1 = ex5_fin_L1.computeInfiniteFinHeat
ex5_heat_L2 = ex5_fin_L2.computeInfiniteFinHeat

ex5_efficiency_L1 = ex5_fin_L1.computeEfficiency
ex5_efficiency_L2 = ex5_fin_L2.computeEfficiency

ex5_adiabaticEffectiveness_L1 = ex5_fin_L1.computeAdiabaticEffectiveness
ex5_adiabaticEffectiveness_L2 = ex5_fin_L2.computeAdiabaticEffectiveness

ex5_convectiveEffectiveness_L1 = ex5_fin_L1.computeConvectiveEffectiveness
ex5_convectiveEffectiveness_L2 = ex5_fin_L2.computeConvectiveEffectiveness

ex5_text = "For L = {} [m] | Adiabatic Heat: {}, Convective Heat: {}, Infinite Heat: {}."

print(ex5_text.format(ex5_L1, ex5_adiabaticQ_L1, ex5_convectiveQ_L1, ex5_heat_L1))
print(ex5_text.format(ex5_L2, ex5_adiabaticQ_L2, ex5_convectiveQ_L2, ex5_heat_L2))

fig = plt.figure()
fig.suptitle("Adiabatic, Convective and Infinite borders.")

ax1 = fig.add_subplot(1, 2, 1)
ax1.plot(np.linspace(0, ex5_L1, N),ex5_adiabaticTemperature_L1, label="Adiabatic Hypothesis")
ax1.plot(np.linspace(0, ex5_L1, N),ex5_convectiveTemperature_L1, label="Convective Hypothesis")
ax1.plot(np.linspace(0, ex5_L1, N),ex5_temperature_L1, label="Infinite Hypothesis")
ax1.legend()
ax1.set_xlabel("Distance [m]")
ax1.set_ylabel("Temperature [ºC]")
ax1.set_title("L = 1 [m]")

ax2 = fig.add_subplot(1, 2, 2)
ax2.plot(np.linspace(0, ex5_L2, N),ex5_adiabaticTemperature_L2, label="Adiabatic Hypothesis")
ax2.plot(np.linspace(0, ex5_L2, N),ex5_convectiveTemperature_L2, label="Convective Hypothesis")
ax2.plot(np.linspace(0, ex5_L2, N),ex5_temperature_L2, label="Infinite Hypothesis")
ax2.legend()
ax2.set_xlabel("Distance [m]")
ax2.set_ylabel("Temperature [ºC]")
ax2.set_title("L = 0.01 [m]")
fig.savefig('output/ex5.svg')

print("*************** EXERCISE: SIX ***************")
ex6_t1 = 1e-2
ex6_t2 = 0.1

ex6_dim_t1 = core.Dimension(L, w, ex6_t1, N)
ex6_dim_t2 = core.Dimension(L, w, ex6_t2, N)

ex6_fin_t1 = core.Fin(ex6_dim_t1, mat, bound)
ex6_fin_t2 = core.Fin(ex6_dim_t2, mat, bound)

ex6_adiabaticTemperature_t1 = ex6_fin_t1.computeAdiabaticTemperature
ex6_adiabaticTemperature_t2 = ex6_fin_t2.computeAdiabaticTemperature

ex6_convectiveTemperature_t1 = ex6_fin_t1.computeConvectiveTemperature
ex6_convectiveTemperature_t2 = ex6_fin_t2.computeConvectiveTemperature

ex6_temperature_t1 = ex6_fin_t1.computeInfiniteFinTemperature
ex6_temperature_t2 = ex6_fin_t2.computeInfiniteFinTemperature

ex6_adiabaticQ_t1 = ex6_fin_t1.computeAdiabaticHeat
ex6_adiabaticQ_t2 = ex6_fin_t2.computeAdiabaticHeat

ex6_convectiveQ_t1 = ex6_fin_t1.computeConvectiveHeat
ex6_convectiveQ_t2 = ex6_fin_t2.computeConvectiveHeat

ex6_heat_t1 = ex6_fin_t1.computeInfiniteFinHeat
ex6_heat_t2 = ex6_fin_t2.computeInfiniteFinHeat

ex6_efficiency_t1 = ex6_fin_t1.computeEfficiency
ex6_efficiency_t2 = ex6_fin_t2.computeEfficiency

ex6_adiabaticEffectiveness_t1 = ex6_fin_t1.computeAdiabaticEffectiveness
ex6_adiabaticEffectiveness_t2 = ex6_fin_t2.computeAdiabaticEffectiveness

ex6_convectiveEffectiveness_t1 = ex6_fin_t1.computeConvectiveEffectiveness
ex6_convectiveEffectiveness_t2 = ex6_fin_t2.computeConvectiveEffectiveness

ex6_text = "For t = {} [m] | Adiabatic Heat: {}, Convective Heat: {}, Infinite Heat: {}."

print(ex6_text.format(ex6_t1, ex6_adiabaticQ_t1, ex6_convectiveQ_t1, ex6_heat_t1))
print(ex6_text.format(ex6_t2, ex6_adiabaticQ_t2, ex6_convectiveQ_t2, ex6_heat_t2))

fig = plt.figure()
fig.suptitle("Adiabatic, Convective and Infinite borders.")

ax1 = fig.add_subplot(1, 2, 1)
ax1.plot(np.linspace(0, L, N),ex6_adiabaticTemperature_t1, label="Adiabatic Hypothesis")
ax1.plot(np.linspace(0, L, N),ex6_convectiveTemperature_t1, label="Convective Hypothesis")
ax1.plot(np.linspace(0, L, N),ex6_temperature_t1, label="Infinite Hypothesis")
ax1.legend()
ax1.set_xlabel("Distance [m]")
ax1.set_ylabel("Temperature [ºC]")
ax1.set_title("t = 1e-2 [m]")

ax2 = fig.add_subplot(1, 2, 2)
ax2.plot(np.linspace(0, L, N),ex6_adiabaticTemperature_t2, label="Adiabatic Hypothesis")
ax2.plot(np.linspace(0, L, N),ex6_convectiveTemperature_t2, label="Convective Hypothesis")
ax2.plot(np.linspace(0, L, N),ex6_temperature_t2, label="Infinite Hypothesis")
ax2.legend()
ax2.set_xlabel("Distance [m]")
ax2.set_ylabel("Temperature [ºC]")
ax2.set_title("t = 0.1 [m]")
fig.savefig('output/ex6.svg')

print("*************** EXERCISE: SEVEN ***************")
# Material properties
ex7_k = 17        # W/m·K

ex7_mat_Ti = core.Material(ex7_k, h)

ex7_fin_Ti = fin
ex7_fin_t2 = core.Fin(dim, ex7_mat_Ti, bound)

ex7_adiabaticTemperature_Ti = ex7_fin_Ti.computeAdiabaticTemperature
ex7_adiabaticTemperature_t2 = ex7_fin_t2.computeAdiabaticTemperature

ex7_convectiveTemperature_Ti = ex7_fin_Ti.computeConvectiveTemperature
ex7_convectiveTemperature_t2 = ex7_fin_t2.computeConvectiveTemperature

ex7_temperature_Ti = ex7_fin_Ti.computeInfiniteFinTemperature
ex7_temperature_t2 = ex7_fin_t2.computeInfiniteFinTemperature

ex7_adiabaticQ_Ti = ex7_fin_Ti.computeAdiabaticHeat
ex7_adiabaticQ_t2 = ex7_fin_t2.computeAdiabaticHeat

ex7_convectiveQ_Ti = ex7_fin_Ti.computeConvectiveHeat
ex7_convectiveQ_t2 = ex7_fin_t2.computeConvectiveHeat

ex7_heat_Ti = ex7_fin_Ti.computeInfiniteFinHeat
ex7_heat_t2 = ex7_fin_t2.computeInfiniteFinHeat

ex7_efficiency_Ti = ex7_fin_Ti.computeEfficiency
ex7_efficiency_t2 = ex7_fin_t2.computeEfficiency

ex7_adiabaticEffectiveness_Ti = ex7_fin_Ti.computeAdiabaticEffectiveness
ex7_adiabaticEffectiveness_t2 = ex7_fin_t2.computeAdiabaticEffectiveness

ex7_convectiveEffectiveness_Ti = ex7_fin_Ti.computeConvectiveEffectiveness
ex7_convectiveEffectiveness_t2 = ex7_fin_t2.computeConvectiveEffectiveness

ex7_text = "Using = {} [m] | Adiabatic Heat: {}, Convective Heat: {}, Infinite Heat: {}."

print(ex7_text.format("Copper [Cu]", ex7_adiabaticQ_Ti, ex7_convectiveQ_Ti, ex7_heat_Ti))
print(ex7_text.format("Titanium [Ti]", ex7_adiabaticQ_t2, ex7_convectiveQ_t2, ex7_heat_t2))

fig = plt.figure()
fig.suptitle("Adiabatic, Convective and Infinite borders.")

ax1 = fig.add_subplot(1, 2, 1)
ax1.plot(np.linspace(0, L, N),ex7_adiabaticTemperature_Ti, label="Adiabatic Hypothesis")
ax1.plot(np.linspace(0, L, N),ex7_convectiveTemperature_Ti, label="Convective Hypothesis")
ax1.plot(np.linspace(0, L, N),ex7_temperature_Ti, label="Infinite Hypothesis")
ax1.legend()
ax1.set_xlabel("Distance [m]")
ax1.set_ylabel("Temperature [ºC]")
ax1.set_title("Copper [Cu]")

ax2 = fig.add_subplot(1, 2, 2)
ax2.plot(np.linspace(0, L, N),ex7_adiabaticTemperature_t2, label="Adiabatic Hypothesis")
ax2.plot(np.linspace(0, L, N),ex7_convectiveTemperature_t2, label="Convective Hypothesis")
ax2.plot(np.linspace(0, L, N),ex7_temperature_t2, label="Infinite Hypothesis")
ax2.legend()
ax2.set_xlabel("Distance [m]")
ax2.set_ylabel("Temperature [ºC]")
ax2.set_title("Titanium [Ti]")
fig.savefig('output/ex7.svg')

print("*************** EXERCISE: EIGHT ***************")
# Material properties
ex8_b1 = 30        # W/m²·K
ex8_b2 = 35        # W/m²·K

ex8_mat_V1 = core.Material(k, ex8_b1)
ex8_mat_V2 = core.Material(k, ex8_b2)

ex8_fin_Ti = core.Fin(dim, ex8_mat_V1, bound)
ex8_fin_t2 = core.Fin(dim, ex8_mat_V2, bound)

ex8_adiabaticTemperature_Ti = ex8_fin_Ti.computeAdiabaticTemperature
ex8_adiabaticTemperature_t2 = ex8_fin_t2.computeAdiabaticTemperature

ex8_convectiveTemperature_Ti = ex8_fin_Ti.computeConvectiveTemperature
ex8_convectiveTemperature_t2 = ex8_fin_t2.computeConvectiveTemperature

ex8_temperature_Ti = ex8_fin_Ti.computeInfiniteFinTemperature
ex8_temperature_t2 = ex8_fin_t2.computeInfiniteFinTemperature

ex8_adiabaticQ_Ti = ex8_fin_Ti.computeAdiabaticHeat
ex8_adiabaticQ_t2 = ex8_fin_t2.computeAdiabaticHeat

ex8_convectiveQ_Ti = ex8_fin_Ti.computeConvectiveHeat
ex8_convectiveQ_t2 = ex8_fin_t2.computeConvectiveHeat

ex8_heat_Ti = ex8_fin_Ti.computeInfiniteFinHeat
ex8_heat_t2 = ex8_fin_t2.computeInfiniteFinHeat

ex8_efficiency_Ti = ex8_fin_Ti.computeEfficiency
ex8_efficiency_t2 = ex8_fin_t2.computeEfficiency

ex8_adiabaticEffectiveness_Ti = ex8_fin_Ti.computeAdiabaticEffectiveness
ex8_adiabaticEffectiveness_t2 = ex8_fin_t2.computeAdiabaticEffectiveness

ex8_convectiveEffectiveness_Ti = ex8_fin_Ti.computeConvectiveEffectiveness
ex8_convectiveEffectiveness_t2 = ex8_fin_t2.computeConvectiveEffectiveness

ex8_text = "For k = {} [W/(m^2·K)] | Adiabatic Heat: {}, Convective Heat: {}, Infinite Heat: {}."

print(ex8_text.format(ex8_b1, ex8_adiabaticQ_Ti, ex8_convectiveQ_Ti, ex8_heat_Ti))
print(ex8_text.format(ex8_b2, ex8_adiabaticQ_t2, ex8_convectiveQ_t2, ex8_heat_t2))

fig = plt.figure()
fig.suptitle("Adiabatic, Convective and Infinite borders.")

ax1 = fig.add_subplot(1, 2, 1)
ax1.plot(np.linspace(0, L, N),ex8_adiabaticTemperature_Ti, label="Adiabatic Hypothesis")
ax1.plot(np.linspace(0, L, N),ex8_convectiveTemperature_Ti, label="Convective Hypothesis")
ax1.plot(np.linspace(0, L, N),ex8_temperature_Ti, label="Infinite Hypothesis")
ax1.legend()
ax1.set_xlabel("Distance [m]")
ax1.set_ylabel("Temperature [ºC]")
ax1.set_title("V = 7.5 [m/s]")

ax2 = fig.add_subplot(1, 2, 2)
ax2.plot(np.linspace(0, L, N),ex8_adiabaticTemperature_t2, label="Adiabatic Hypothesis")
ax2.plot(np.linspace(0, L, N),ex8_convectiveTemperature_t2, label="Convective Hypothesis")
ax2.plot(np.linspace(0, L, N),ex8_temperature_t2, label="Infinite Hypothesis")
ax2.legend()
ax2.set_xlabel("Distance [m]")
ax2.set_ylabel("Temperature [ºC]")
ax2.set_title("V = 20 [m/s]")
fig.savefig('output/ex8.svg')

print("*************** EXERCISE: NINE ***************")

ex9_mL1 = 1
ex9_mL2 = 3
ex9_mL3 = 8

ex9_fin_mL1 = core.Fin(dim, mat, bound, ex9_mL1)
ex9_fin_mL2 = core.Fin(dim, mat, bound, ex9_mL2)
ex9_fin_mL3 = core.Fin(dim, mat, bound, ex9_mL3)

ex9_adiabaticTemperature_mL1 = ex9_fin_mL1.computeAdiabaticTemperature
ex9_adiabaticTemperature_mL2 = ex9_fin_mL2.computeAdiabaticTemperature
ex9_adiabaticTemperature_mL3 = ex9_fin_mL3.computeAdiabaticTemperature
ex9_infiniteTemperatrue = fin.computeInfiniteFinTemperature

ex9_text = "Using = {} | Adiabatic Heat: {}, Convective Heat: {}, Infinite Heat: {}."

fig = plt.figure()
ax1 = fig.add_subplot(1, 1, 1)
ax1.set_title("Adiabatic, Convective and Infinite borders.")
ax1.plot(np.linspace(0, L, N),ex9_adiabaticTemperature_mL1, label="Adiabatic Hypothesis | mL = 1")
ax1.plot(np.linspace(0, L, N),ex9_adiabaticTemperature_mL2, label="Adiabatic Hypothesis | mL = 2")
ax1.plot(np.linspace(0, L, N),ex9_adiabaticTemperature_mL3, label="Adiabatic Hypothesis | mL = 8")
ax1.plot(np.linspace(0, L, N),ex9_infiniteTemperatrue, label="Adiabatic Hypothesis | mL = ∞ ")
ax1.legend()
ax1.set_xlabel("Distance [m]")
ax1.set_ylabel("Temperature [ºC]")
fig.savefig('output/ex9.svg')

plt.show()