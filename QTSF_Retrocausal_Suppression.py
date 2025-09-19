"""
Retrocausal Suppression via Decoherence Simulation
Converted from MATLAB/Octave -> Python (numpy + matplotlib)

Author: Tiziano Demaria
"""

import numpy as np
import matplotlib.pyplot as plt

# ---------------------------
# Helper: safe exponent (avoids underflow warnings)
# ---------------------------
def safe_exp(x, lower_clip=-700.0, upper_clip=700.0):
    """Compute exp(x) while clipping x to avoid overflow/underflow warnings.
    Default clips at ±700 which keeps values finite and tiny but non-zero."""
    return np.exp(np.clip(x, lower_clip, upper_clip))


# ---------------------------
# Parameters
# ---------------------------
gamma_range = np.linspace(0.01, 0.2, 6)  # Decoherence rates
N_env = 50                               # Number of environmental degrees of freedom
t_max = 50.0                             # Maximum time
dt = 0.1                                 # Time step
t = np.arange(0.0, t_max + dt, dt)       # time vector

beta_retro = 0.1                         # Retrocausal coupling strength (small off-diagonal term)

# Preallocate results
beta_eff = np.zeros((len(gamma_range), len(t)))
purity = np.zeros((len(gamma_range), len(t)))

# Color palette
colors = plt.get_cmap('tab10')(np.linspace(0, 1, len(gamma_range)))

# ---------------------------
# Main simulation loop
# ---------------------------
for g_idx, gamma in enumerate(gamma_range):
    # Initial |+> state
    psi0 = np.array([1.0, 1.0]) / np.sqrt(2.0)
    rho0 = np.outer(psi0, psi0.conj())  # 2x2 density matrix

    # Start with initial rho (we'll replace entries each timestep as in original model)
    rho = rho0.copy()

    for t_idx, current_time in enumerate(t):
        # Decoherence effect: exponential decay of off-diagonals (with safe exp)
        decoherence_exponent = -gamma * current_time * N_env
        decoherence_factor = safe_exp(decoherence_exponent)

        # Effective retrocausal term (suppressed by decoherence)
        beta_eff[g_idx, t_idx] = beta_retro * decoherence_factor

        # Update density matrix (simple toy model)
        rho[0, 1] = rho0[0, 1] * decoherence_factor
        rho[1, 0] = rho0[1, 0] * decoherence_factor
        rho[0, 0] = 0.5 * (1.0 + (1.0 - decoherence_factor))
        rho[1, 1] = 0.5 * (1.0 - (1.0 - decoherence_factor))

        # Purity Tr(rho^2)
        purity[g_idx, t_idx] = np.real_if_close(np.trace(rho @ rho))

# ---------------------------
# PLOTTING: 2x2 subplot layout
# ---------------------------
fig, axes = plt.subplots(2, 2, figsize=(14, 10))
ax1, ax2, ax3, ax4 = axes.flat
fig.suptitle('Retrocausal Suppression via Environmental Decoherence', fontsize=16, fontweight='bold')

# Plot 1: Suppression of retrocausal term over time (log y-scale)
for g_idx, gamma in enumerate(gamma_range):
    ax1.plot(t, beta_eff[g_idx, :], linewidth=2, color=colors[g_idx],
             label=rf'$\gamma$ = {gamma:.2f}')
ax1.set_xlabel('Time')
ax1.set_ylabel(r'Effective Retrocausal Coupling $\beta_{\mathrm{eff}}$')
ax1.set_title('Suppression of Retrocausal Effects by Decoherence')
ax1.legend(loc='upper right')
ax1.grid(True)
ax1.set_yscale('log')

# Plot 2: Purity decay over time
for g_idx, gamma in enumerate(gamma_range):
    ax2.plot(t, purity[g_idx, :], linewidth=2, color=colors[g_idx],
             label=rf'$\gamma$ = {gamma:.2f}')
ax2.set_xlabel('Time')
ax2.set_ylabel(r'Purity $\mathrm{Tr}(\rho^2)$')
ax2.set_title('Decoherence Process: Loss of Quantum Coherence')
ax2.legend(loc='lower left')
ax2.grid(True)

# Plot 3: Exponential suppression vs environment size (semilogy)
N_range = np.arange(1, 101)  # 1 ... 100
gamma_fixed = 0.1
time_fixed = 10.0

# compute suppression safely (clip exponent)
suppression = safe_exp(-gamma_fixed * time_fixed * N_range)

ax3.semilogy(N_range, suppression, 'b-', linewidth=2, label='Decoherence Factor')
ax3.semilogy(N_range, beta_retro * suppression, 'r--', linewidth=2, label=r'$\beta_{\mathrm{eff}} = \beta \times$ Decoherence')

ax3.set_xlabel('Environmental Degrees of Freedom (N)')
ax3.set_ylabel('Suppression Factor')
ax3.set_title('Exponential Suppression with System Size')
ax3.grid(True)
ax3.legend(loc='lower left')

# Add markers for specific system sizes (including macroscopic)
macroscopic_N = 1e23
marker_N = np.array([1.0, 10.0, 100.0, macroscopic_N])
marker_suppression = safe_exp(-gamma_fixed * time_fixed * marker_N)

# To avoid plotting true zeros on log axes, the safe_exp clipped values are non-zero and tiny.
ax3.loglog(marker_N, marker_suppression, 'ko', markerfacecolor='k', markersize=6)
# annotate near the 100 marker
ax3.text(marker_N[2] * 1.1, marker_suppression[2] * 5.0, 'Macroscopic', fontsize=9)

# Plot 4: Suppression vs Entropy/Complexity
S_range = np.linspace(0.0, 10.0, 200)
beta_suppressed = beta_retro * safe_exp(-S_range)

ax4.plot(S_range, beta_suppressed, 'm-', linewidth=2)
ax4.set_xlabel('Entropy S (or log(N))')
ax4.set_ylabel(r'Effective Retrocausal Coupling $\beta_{\mathrm{eff}}$')
ax4.set_title('Suppression vs Entropy/Complexity')
ax4.grid(True)
ax4.set_yscale('log')

# Add regime labels & shading
quantum_cutoff = 1.0
classical_start_idx = np.searchsorted(S_range, quantum_cutoff)
if classical_start_idx < len(S_range):
    ax4.fill_between(S_range[:classical_start_idx], beta_suppressed[:classical_start_idx],
                     color='c', alpha=0.2, edgecolor='none')
    ax4.fill_between(S_range[classical_start_idx:], beta_suppressed[classical_start_idx:],
                     color='y', alpha=0.2, edgecolor='none')
    # place labels
    ax4.text(np.mean(S_range[:classical_start_idx]), beta_retro / 2,
             'Quantum Regime', ha='center')
    ax4.text(np.mean(S_range[classical_start_idx:]), beta_retro / 100,
             'Classical Regime', ha='center')

plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.show()

# ---------------------------
# Display key results (prints)
# ---------------------------
print('Retrocausal Suppression via Decoherence - Simulation Results:')
print(f'Initial retrocausal coupling: beta = {beta_retro:.3f}\n')

print('Suppression at different scales (gamma_fixed = {:.2f}, time_fixed = {:.1f}):'
      .format(gamma_fixed, time_fixed))
def compute_beta_eff_for_N(N):
    return beta_retro * safe_exp(-gamma_fixed * time_fixed * N)

print('Microscopic (N=1): beta_eff ≈ {:.3f}'.format(compute_beta_eff_for_N(1.0)))
print('Mesoscopic (N=1e3): beta_eff ≈ {:.3e}'.format(compute_beta_eff_for_N(1e3)))
print('Macroscopic (N=1e23): beta_eff ≈ {:.3e}'.format(compute_beta_eff_for_N(1e23)))
print('Suppression factor (macroscopic): {:.1e}'.format(safe_exp(-gamma_fixed * time_fixed * 1e23)))
