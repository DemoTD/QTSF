"""
k-essence Lagrangian Simulation
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

# Equation of state parameter w = p/ρ
w = p / rho

# Alternative direct calculation
w_direct = (Y - M) / (Y + M)

# Create figure
fig, axs = plt.subplots(2, 2, figsize=(12, 8))

# ---- Plot 1: Equation of state ----
axs[0, 0].plot(Y, w, 'b-', linewidth=2, label='w = p/ρ')
axs[0, 0].plot(Y, w_direct, 'r--', linewidth=1.5, label='w = (Y-M)/(Y+M)')
axs[0, 0].axhline(0, color='k', linestyle='--', label='w = 0 (Dust/Matter)')
axs[0, 0].axhline(1, color='k', linestyle='--', label='w = 1 (Stiff Fluid)')
axs[0, 0].axvline(M, color='k', linestyle='--', label='Y = M')
axs[0, 0].set_xlabel('Y = √(2X)')
axs[0, 0].set_ylabel('Equation of State Parameter w')
axs[0, 0].set_title('K-essence Equation of State: w = p/ρ')
axs[0, 0].legend(loc='southeast'.replace('southeast','lower right'))
axs[0, 0].grid(True)

# ---- Plot 2: Energy density ----
axs[0, 1].semilogy(Y, rho, 'r-', linewidth=2)
axs[0, 1].axvline(M, color='k', linestyle='--', label='Y = M')
axs[0, 1].set_xlabel('Y = √(2X)')
axs[0, 1].set_ylabel('Energy Density ρ')
axs[0, 1].set_title('Energy Density vs Y')
axs[0, 1].grid(True)

# ---- Plot 3: Pressure ----
axs[1, 0].plot(Y, p, 'g-', linewidth=2)
axs[1, 0].axvline(M, color='k', linestyle='--', label='Y = M')
axs[1, 0].set_xlabel('Y = √(2X)')
axs[1, 0].set_ylabel('Pressure p')
axs[1, 0].set_title('Pressure vs Y')
axs[1, 0].grid(True)

# ---- Plot 4: Zoom near Y = M ----
Y_zoom = np.linspace(M*0.8, M*1.2, 200)
w_zoom = (Y_zoom - M) / (Y_zoom + M)
axs[1, 1].plot(Y_zoom, w_zoom, 'm-', linewidth=2)
axs[1, 1].axvline(M, color='k', linestyle='--', label='Y = M')
axs[1, 1].axhline(0, color='k', linestyle='--', label='w = 0')
axs[1, 1].set_xlabel('Y = √(2X)')
axs[1, 1].set_ylabel('Equation of State Parameter w')
axs[1, 1].set_title('Zoomed View Near Y = M')
axs[1, 1].grid(True)

# Add overall title
fig.suptitle('K-essence Dynamics: P(X) = α(√(2X) - M)²', fontsize=14, fontweight='bold')
plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.show()

# ---- Display key results ----
print("K-essence Dynamics Simulation Results:")
print("Lagrangian: P(X) = α(√(2X) - M)²")
print("Equation of State: w = (Y - M)/(Y + M)\n")
print(f"At Y = M: w = {(M - M)/(M + M):.4f} (Dark Matter regime)")
print(f"At Y = 5M: w = {(5*M - M)/(5*M + M):.4f} (approaching stiff fluid)")
print(f"At Y = 0.1M: w = {(0.1*M - M)/(0.1*M + M):.4f}")
