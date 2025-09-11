"""
Numerical Solution for Structure Growth with Different Equations of State

This code numerically integrates the differential equation for the growth of density
perturbations in a universe with matter and dark energy. It compares the growth
factors for standard Cold Dark Matter (w=0) and the proposed T-field gradient
component (w=-1/3) from the paper "Quantum Time as a Scalar Field".

Author: Tiziano Demaria
Date: [Current Date]
License: MIT
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# Define the Hubble parameter: H(a)^2 / H0^2 = Om/a^3 + Ol
def H2_over_H02(a, Om=0.3):
    Ol = 1 - Om
    return Om / a**3 + Ol

# Define the growth equation for a component with equation of state w
def growth_equation(y, a, w, Om):
    delta, ddelta_da = y
    # Convert derivatives wrt a to derivatives wrt time
    # Second order ODE: d²δ/da² + [dlnH/da + 3/a] dδ/da - [3/2 Om(a) (1-w)(1+3w) / (a^5 H(a)^2/H0^2)] δ = 0
    H2 = H2_over_H02(a, Om)
    dlnH_da = ( -3*Om/(a**4) ) / (2 * H2) # Derivative of ln(H) wrt a
    coefficient = (3/2) * Om * (1 - w) * (1 + 3*w) / (a**5 * H2)
    d2delta_da2 = - (dlnH_da + 3/a) * ddelta_da + coefficient * delta
    return [ddelta_da, d2delta_da2]

# Integration parameters
a_initial = 0.01
a_final = 1.0
a_values = np.linspace(a_initial, a_final, 1000)

# Initial conditions: delta=1, ddelta/da=0 at a_initial
y_initial = [1, 0]

# Integrate for w = 0 (CDM)
solution_cdm = odeint(growth_equation, y_initial, a_values, args=(0.0, 0.3))
delta_cdm = solution_cdm[:, 0]
D_cdm = delta_cdm / delta_cdm[0] # Growth factor D(a) = δ(a)/δ(initial)

# Integrate for w = -1/3 (T-field gradient energy)
solution_w = odeint(growth_equation, y_initial, a_values, args=(-1/3, 0.3))
delta_w = solution_w[:, 0]
D_w = delta_w / delta_w[0] # Growth factor D(a)

# Plot
plt.figure(figsize=(8, 5))
plt.plot(a_values, D_cdm, 'r--', label='CDM ($w=0$)')
plt.plot(a_values, D_w, 'b-', label='$T$-field Gradients ($w=-1/3$)')
plt.xlabel('Scale Factor ($a$)')
plt.ylabel('Growth Factor ($D(a) = \\delta(a)/\\delta(a_{initial})$')
plt.title('Growth of Structure: Comparison of Equations of State')
plt.legend()
plt.grid(True)
plt.show()
