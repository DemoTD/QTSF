"""
Retrocausal Suppression - Time Evolution Figures
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
import numpy as np
import matplotlib.pyplot as plt

# Parameters
gamma_range = np.linspace(0.01, 0.2, 6)   # Decoherence rates
N_env = 50                                # Environmental degrees of freedom
t_max = 50                                # Maximum time
dt = 0.1                                  # Time step
t = np.arange(0, t_max + dt, dt)

# Retrocausal coupling strength (small off-diagonal term)
beta_retro = 0.1

# Initialize results
beta_eff = np.zeros((len(gamma_range), len(t)))
purity = np.zeros((len(gamma_range), len(t)))

# Colors for plotting (Matlab's "lines" equivalent)
colors = plt.cm.tab10(np.linspace(0, 1, len(gamma_range)))

# Main simulation loop
for g_idx, gamma in enumerate(gamma_range):
    # Initial state: superposition state |+>
    psi0 = np.array([1, 1]) / np.sqrt(2)
    rho0 = np.outer(psi0, psi0.conj())  # Density matrix
    rho = rho0.copy()

    for t_idx, current_time in enumerate(t):
        # Decoherence factor
        decoherence_factor = np.exp(-gamma * current_time * N_env)

        # Effective retrocausal coupling
        beta_eff[g_idx, t_idx] = beta_retro * decoherence_factor

        # Update density matrix (simplified model)
        rho[0, 1] = rho0[0, 1] * decoherence_factor
        rho[1, 0] = rho0[1, 0] * decoherence_factor
        rho[0, 0] = 0.5 * (1 + (1 - decoherence_factor))
        rho[1, 1] = 0.5 * (1 - (1 - decoherence_factor))

        # Purity Tr(ρ²)
        purity[g_idx, t_idx] = np.trace(rho @ rho).real

# --- Plotting ---
fig, axs = plt.subplots(2, 1, figsize=(10, 8))

# Plot 1: Suppression of retrocausal term
for g_idx, gamma in enumerate(gamma_range):
    axs[0].plot(t, beta_eff[g_idx, :], linewidth=2.5, color=colors[g_idx],
                label=f'γ = {gamma:.2f}')
axs[0].set_xlabel('Time')
axs[0].set_ylabel('Effective Retrocausal Coupling β_eff')
axs[0].set_title('(a) Time Evolution of Retrocausal Suppression')
axs[0].set_yscale('log')
axs[0].legend(loc='upper right')
axs[0].grid(True)
axs[0].tick_params(labelsize=12)

# Plot 2: Purity decay
for g_idx, gamma in enumerate(gamma_range):
    axs[1].plot(t, purity[g_idx, :], linewidth=2.5, color=colors[g_idx],
                label=f'γ = {gamma:.2f}')
axs[1].set_xlabel('Time')
axs[1].set_ylabel('Purity Tr(ρ²)')
axs[1].set_title('(b) Decoherence Process: Loss of Quantum Coherence')
axs[1].legend(loc='lower left')
axs[1].grid(True)
axs[1].tick_params(labelsize=12)

plt.tight_layout()
plt.savefig("retrocausal_time_evolution.png", dpi=300)
plt.show()
