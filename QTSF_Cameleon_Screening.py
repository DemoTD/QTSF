# ----------------------------------
# MathLab to adapted Python' script
# Cameleon Screening
# ----------------------------------
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
from matplotlib.colors import LogNorm

# -----------------------------
# Parameters (natural units)
# -----------------------------
Lambda = 1.0    # Energy scale
lambda_c = 1.0  # Coupling constant
n = 2           # Exponent in the potential

# Density range (log scale: cosmic → Earth)
rho_cosmo = 1e-27
rho_earth = 1e0
n_rho = 20
rho_m = np.logspace(np.log10(rho_cosmo), np.log10(rho_earth), n_rho)

# Field range
T_min_range = 0.1
T_max_range = 10
n_T = 1000
T = np.linspace(T_min_range, T_max_range, n_T)

# Arrays for minima and mass
T_minima = np.zeros_like(rho_m)
m_eff = np.zeros_like(rho_m)
V_eff_min = np.zeros_like(rho_m)

# Colormap
cmap = get_cmap('jet', n_rho)
colors = [cmap(i) for i in range(n_rho)]

# -----------------------------
# PLOT 1: Effective Potential
# -----------------------------
fig1, ax1 = plt.subplots(figsize=(10, 7))

for i, rho in enumerate(rho_m):
    V_pot = Lambda**4 * np.exp((Lambda**n) / (T**n))
    V_coupling = lambda_c * T * rho
    V_eff = V_pot + V_coupling

    # Minimum
    idx = np.argmin(V_eff)
    T_min = T[idx]
    V_min = V_eff[idx]
    T_minima[i] = T_min
    V_eff_min[i] = V_min

    # Plot selected densities
    if i % 4 == 0 or i == n_rho - 1:
        ax1.plot(T, V_eff, color=colors[i], lw=1.5, label=f'ρ_m = {rho:.1e}')
        ax1.plot(T_min, V_min, 'o', color=colors[i], ms=6, mfc=colors[i])

ax1.set_xlabel("Field Value T")
ax1.set_ylabel(r"Effective Potential $V_{eff}(T)$")
ax1.set_title("Chameleon Effective Potential for Different Matter Densities")
ax1.set_yscale("log")
ax1.grid(True, which="both")
ax1.legend(loc="upper left")

# Colorbar: use LogNorm so ticks are in density (not log-density)
sm = plt.cm.ScalarMappable(cmap='jet', norm=LogNorm(vmin=rho_cosmo, vmax=rho_earth))
sm.set_array([])  # required for colorbar to work with a ScalarMappable not attached to data
cbar_ticks = np.logspace(np.log10(rho_cosmo), np.log10(rho_earth), 5)
cbar = fig1.colorbar(sm, ax=ax1, ticks=cbar_ticks)
cbar.set_ticklabels([f"{t:.1e}" for t in cbar_ticks])
cbar.set_label("Matter Density ρ_m")

fig1.suptitle("Chameleon Screening Mechanism:\n"
              r"$V_{eff}(T) = \Lambda^4 \exp(\Lambda^n/T^n) + \lambda T \rho_m$")
fig1.tight_layout(rect=[0, 0, 1, 0.95])

# -----------------------------
# Compute Effective Mass
# -----------------------------
for i, rho in enumerate(rho_m):
    T_min = T_minima[i]
    # avoid negative/zero when T_min is very small
    delta = min(1e-4, 0.1 * max(T_min, 1e-8))
    if T_min - delta <= 0:
        delta = 0.5 * T_min  # ensure positive
        if delta <= 0:
            delta = 1e-8

    V_plus = Lambda**4 * np.exp((Lambda**n) / ((T_min + delta)**n)) + lambda_c * (T_min + delta) * rho
    V_min_val = Lambda**4 * np.exp((Lambda**n) / (T_min**n)) + lambda_c * T_min * rho
    V_minus = Lambda**4 * np.exp((Lambda**n) / ((T_min - delta)**n)) + lambda_c * (T_min - delta) * rho

    # Finite difference second derivative
    m_eff_sq = (V_plus - 2*V_min_val + V_minus) / (delta**2)
    m_eff[i] = np.sqrt(abs(m_eff_sq))

# -----------------------------
# PLOT 2 & 3: Mass and Range
# -----------------------------
fig2, axes = plt.subplots(1, 2, figsize=(10, 6))
fig2.suptitle("Chameleon Field Properties vs. Density")

# Plot 2: Mass vs density
ax = axes[0]
ax.loglog(rho_m, m_eff, 'b-o', lw=2, ms=4)
ax.set_xlabel("Matter Density ρ_m")
ax.set_ylabel("Effective Mass m_eff")
ax.set_title("Chameleon Mass vs Environmental Density")
ax.grid(True, which="both")

# Reference densities
rho_galaxy = 1e-24
rho_solar = 1e-18
rho_lab = 1e-12
ymin, ymax = ax.get_ylim()
ax.axvline(rho_galaxy, color="r", ls="--")
ax.axvline(rho_solar, color="g", ls="--")
ax.axvline(rho_lab, color="m", ls="--")
ax.text(rho_galaxy, ymax/5, "Galactic", color="r", fontsize=8)
ax.text(rho_solar, ymax/5, "Solar System", color="g", fontsize=8)
ax.text(rho_lab, ymax/5, "Laboratory", color="m", fontsize=8)

# Plot 3: Compton wavelength (force range)
ax2 = axes[1]
compton_range = 1.0 / m_eff
ax2.loglog(rho_m, compton_range, 'k-s', lw=2, ms=4, mfc='yellow')
ax2.set_xlabel("Matter Density ρ_m")
ax2.set_ylabel("Compton Wavelength (Range)")
ax2.set_title("Chameleon Range vs Environmental Density")
ax2.grid(True, which="both")

fig2.tight_layout(rect=[0, 0, 1, 0.95])
plt.show()
