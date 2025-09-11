"""
Generate Figure A1: Required Underdensity to Resolve Hubble Tension

This code produces the visualization for the relationship between the fraction of
dark matter constituted by T-field gradients (f) and the required underdensity (δ)
to explain the observed Hubble tension (ΔH₀/H₀ ≈ 0.08).

The relationship is derived from the equation in Appendix A1 of the paper:
δ_∇T ≈ 0.19 / f

Author: Tiziano Demaria
Date: [Current Date]
"""

import matplotlib.pyplot as plt
import numpy as np

# Parameters from the paper's calculation
H0_DISCREPANCY = 0.08  # Observed ΔH₀/H₀
CONSTANT_TERM = 0.19   # Derived constant: 0.19 = 0.08 / (0.42 * (Ω_DM/Ω_m))

# Create array of f values (fraction of DM from T-field gradients)
f = np.linspace(0.1, 1.0, 100)

# Calculate required underdensity: δ ≈ 0.19 / f
required_underdensity = CONSTANT_TERM / f

# Create figure
plt.figure(figsize=(10, 6))

# Plot the main relationship
plt.plot(f, required_underdensity, linewidth=3, color='royalblue', 
         label='Model: |δ∇T| = 0.19/f')

# Highlight specific values discussed in the paper
plt.scatter([1.0, 0.5], [0.19, 0.38], color='crimson', s=100, zorder=5, 
            label='Specific solutions')
plt.text(1.02, 0.17, "f=1.0, δ=-19%", fontsize=11, ha='left')
plt.text(0.52, 0.36, "f=0.5, δ=-38%", fontsize=11, ha='left')

# Add a horizontal line for the observed Hubble tension
plt.axhline(y=H0_DISCREPANCY, color='green', linestyle='--', alpha=0.7, 
            label='Observed ΔH₀/H₀ = 0.08')

# Customize the plot
plt.xlabel('Fraction of Dark Matter from T-field Gradients (f)', fontsize=12)
plt.ylabel('Required Underdensity (|δ∇T|)', fontsize=12)
plt.title('Required Underdensity to Resolve Hubble Tension', fontsize=14)
plt.grid(True, alpha=0.3)
plt.legend(loc='upper right')

# Set axis limits and ticks
plt.xlim(0.1, 1.0)
plt.ylim(0.0, 0.5)
plt.xticks(np.arange(0.1, 1.1, 0.1))

# Add explanatory text box
textstr = '\n'.join([
    'Model prediction: |δ∇T| = 0.19/f',
    '',
    'For the observed Hubble tension',
    'ΔH₀/H₀ ≈ 0.08, the required underdensity',
    'depends on the fraction f of dark matter',
    'constituted by T-field gradients.'
])
props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
plt.text(0.15, 0.40, textstr, fontsize=11, bbox=props, va='top')

plt.tight_layout()

# Save the figure (uncomment to save)
# plt.savefig('required_underdensity_plot.png', dpi=300, bbox_inches='tight')

plt.show()

# Print the derivation for clarity
print("Derivation of the relationship:")
print("From Appendix A1: ΔH₀/H₀ ≈ 0.42 · f · δ∇T")
print(f"With observed ΔH₀/H₀ = {H0_DISCREPANCY}:")
print(f"{H0_DISCREPANCY} ≈ 0.42 · f · δ∇T")
print("Therefore: δ∇T ≈ 0.08 / (0.42 · f) ≈ 0.19 / f")