"""
Retrocausal Suppression - Scaling Figures
Converted from MATLAB/Octave -> Python (numpy + matplotlib)

Author: Tiziano Demaria
Date: 19 september 2025
License: MIT
"""
# Retrocausal Suppression - Scaling Figures

import numpy as np
import matplotlib.pyplot as plt

# Parameters
beta_retro = 0.1      # Initial retrocausal coupling
gamma_fixed = 0.1     # Fixed decoherence rate
time_fixed = 10       # Fixed time

# Create figure for scaling behavior
fig, axes = plt.subplots(2, 1, figsize=(10, 8))

# -------------------------------------------------
# Plot 1: Exponential suppression vs environment size
# -------------------------------------------------
N_range = np.logspace(0, 4, 100)  # From 1 to 10^4
suppression = np.exp(-gamma_fixed * time_fixed * N_range)

axes[0].loglog(N_range, suppression, 'b-', linewidth=2.5, label='Decoherence Factor')
axes[0].loglog(N_range, beta_retro * suppression, 'r--', linewidth=2.5, label=r'$\beta_{eff} = \beta \times$ Decoherence')

axes[0].set_xlabel('Environmental Degrees of Freedom (N)')
axes[0].set_ylabel('Magnitude')
axes[0].set_title('(c) Exponential Suppression with System Size')
axes[0].legend(loc='lower left')
axes[0].grid(True, which="both", ls="--")

# Add markers and labels
micro_N = 1
meso_N = 1e3
axes[0].plot(micro_N, np.exp(-gamma_fixed * time_fixed * micro_N), 'ko', markersize=10, markerfacecolor='k')
axes[0].plot(meso_N, np.exp(-gamma_fixed * time_fixed * meso_N), 'ko', markersize=10, markerfacecolor='k')

axes[0].text(micro_N*1.5, 1e-1, 'Microscopic', fontsize=10, fontweight='bold')
axes[0].text(meso_N*0.8, 1e-10, 'Mesoscopic', fontsize=10, fontweight='bold')

# -------------------------------------------------
# Plot 2: Quantum vs Classical regime
# -------------------------------------------------
S_range = np.linspace(0, 15, 100)
beta_suppressed = beta_retro * np.exp(-S_range)

axes[1].semilogy(S_range, beta_suppressed, 'm-', linewidth=2.5)
axes[1].set_xlabel('Entropy S (or log(N))')
axes[1].set_ylabel(r'Effective Retrocausal Coupling $\beta_{eff}$')
axes[1].set_title('(d) Suppression vs Entropy/Complexity')
axes[1].grid(True, which="both", ls="--")

# Shading for quantum/classical regimes
quantum_cutoff = 3
classical_start = np.argmax(S_range > quantum_cutoff)

if classical_start > 0:
    # Quantum regime shading
    axes[1].fill_between([S_range[0], S_range[classical_start]],
                         1e-20, 1, color=[0.8, 0.9, 1], alpha=0.3)
    # Classical regime shading
    axes[1].fill_between([S_range[classical_start], S_range[-1]],
                         1e-20, 1, color=[1, 0.9, 0.8], alpha=0.3)

    # Labels
    axes[1].text(np.mean(S_range[:classical_start]), beta_retro/5, 'Quantum Regime',
                 fontsize=11, fontweight='bold', ha='center')
    axes[1].text(np.mean(S_range[classical_start:]), beta_retro/500, 'Classical Regime',
                 fontsize=11, fontweight='bold', ha='center')

    # Vertical line at transition
    axes[1].plot([quantum_cutoff, quantum_cutoff], [1e-20, 1], 'k--', linewidth=1.5)
    axes[1].annotate('Quantum-Classical Transition',
                 xy=(quantum_cutoff, beta_suppressed[classical_start]),
                 xytext=(quantum_cutoff+1, beta_retro/2),
                 arrowprops=dict(arrowstyle='->', color='black'),
                 fontsize=9)

# -------------------------------------------------
# Save as high-resolution PNG
# -------------------------------------------------
plt.tight_layout()
plt.savefig('retrocausal_scaling.png', dpi=300)
plt.show()

# -------------------------------------------------
# Display key results
# -------------------------------------------------
print("Retrocausal Suppression Scaling Results:")
print(f"Initial retrocausal coupling: β = {beta_retro:.3f}\n")
print("Suppression at different scales:")
print(f"Microscopic (N=1): β_eff ≈ {beta_retro * np.exp(-gamma_fixed * time_fixed * 1):.3f}")
print(f"Mesoscopic (N=1e3): β_eff ≈ {beta_retro * np.exp(-gamma_fixed * time_fixed * 1e3):.3e}")
print(f"Macroscopic (N=1e23): β_eff ≈ {beta_retro * np.exp(-gamma_fixed * time_fixed * 1e23):.3e}")
print(f"Suppression factor: {np.exp(-gamma_fixed * time_fixed * 1e23):.1e}")
