"""
K-essence Dynamics Simulation
Converted from MATLAB/Octave -> Python (numpy + matplotlib)

Author: Tiziano Demaria
Date: 19 september 2025
License: MIT
"""
import numpy as np
import matplotlib.pyplot as plt

# Parameters
alpha = 1       # Arbitrary scaling parameter (doesn't affect w)
M = 1           # Mass scale parameter
Y_min = 0.1     # Minimum Y value
Y_max = 5       # Maximum Y value
n_points = 1000 # Number of points

# Create Y array
Y = np.linspace(Y_min, Y_max, n_points)

# Calculate energy density ρ and pressure p
rho = alpha * (Y - M)**2 * (3 - (Y - M)/Y)
p = alpha * (Y - M)**2

# Calculate equation of state parameter w
w = p / rho

# Alternative direct calculation: w = (Y - M)/(Y + M)
w_direct = (Y - M) / (Y + M)

# Create figure (single large plot)
plt.figure(figsize=(10, 7))
plt.plot(Y, w, 'b-', linewidth=2, label='w = p/ρ')
plt.plot(Y, w_direct, 'r--', linewidth=1.5, label='w = (Y-M)/(Y+M)')

# Labels and title
plt.xlabel('Y = √(2X)', fontsize=12)
plt.ylabel('Equation of State Parameter w', fontsize=12)
plt.title('K-essence Equation of State: w = p/ρ', fontsize=14, fontweight='bold')
plt.legend(loc='lower right', fontsize=10)
plt.grid(True)

# Reference lines
plt.axhline(0, color='k', linestyle='--', label='w = 0 (Dust/Matter)')
plt.axhline(1, color='k', linestyle='--', label='w = 1 (Stiff Fluid)')
plt.axvline(M, color='k', linestyle='--', label='Y = M')

plt.show()

# Display key results
print("K-essence Dynamics Simulation Results:")
print("Lagrangian: P(X) = α(√(2X) - M)²")
print("Equation of State: w = (Y - M)/(Y + M)\n")
print(f"At Y = M: w = {(M - M)/(M + M):.4f} (Dark Matter regime)")
print(f"At Y = 5M: w = {(5*M - M)/(5*M + M):.4f} (approaching stiff fluid)")
print(f"At Y = 0.1M: w = {(0.1*M - M)/(0.1*M + M):.4f}")
