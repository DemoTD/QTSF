% Chameleon Screening Mechanism Simulation - OCTAVE COMPATIBLE VERSION
% Demonstrates how the effective mass changes with matter density

clear; clc; close all;

% Define parameters (using natural units where appropriate)
Lambda = 1.0;       % Energy scale of the potential
lambda = 1.0;       % Coupling constant to matter
n = 2;              % Exponent in the potential

% Define density range (logarithmic scale from cosmic to terrestrial densities)
rho_cosmo = 1e-27;  % Cosmological density scale
rho_earth = 1e0;    % Earth density scale
n_rho = 20;         % Number of density values
rho_m = logspace(log10(rho_cosmo), log10(rho_earth), n_rho);

% Field range to explore
T_min_range = 0.1;
T_max_range = 10;
n_T = 1000;
T = linspace(T_min_range, T_max_range, n_T);

% Preallocate arrays
T_minima = zeros(size(rho_m));
m_eff = zeros(size(rho_m));
V_eff_min = zeros(size(rho_m));

% Colors for different densities (using jet instead of parula for Octave)
colors = jet(n_rho);

% --- PLOT 1: EFFECTIVE POTENTIAL (Full figure) ---
figure('Position', [100, 100, 1000, 700]); % <--- New, dedicated figure
hold on;

for i = 1:n_rho
    % Calculate effective potential
    V_pot = Lambda^4 * exp((Lambda^n)./(T.^n));
    V_coupling = lambda * T * rho_m(i);
    V_eff = V_pot + V_coupling;

    % Find minimum numerically
    [V_min, idx] = min(V_eff);
    T_min = T(idx);
    T_minima(i) = T_min;
    V_eff_min(i) = V_min;

    % Plot only selected densities for clarity
    if mod(i, 4) == 1 || i == n_rho
        plot(T, V_eff, 'Color', colors(i, :), 'LineWidth', 1.5, ...
            'DisplayName', sprintf('ρ_m = %.1e', rho_m(i)));
        plot(T_min, V_min, 'o', 'Color', colors(i, :), 'MarkerSize', 6, ...
            'MarkerFaceColor', colors(i, :));
    end
end

xlabel('Field Value T');
ylabel('Effective Potential V_{eff}(T)');
title('Chameleon Effective Potential for Different Matter Densities');
legend('Location', 'northwest');
grid on;
set(gca, 'YScale', 'log');
c = colorbar('Ticks', linspace(0, 1, 5), ...
             'TickLabels', arrayfun(@(x) sprintf('%.1e', x), ...
             linspace(rho_cosmo, rho_earth, 5), 'UniformOutput', false));
colormap(jet);
ylabel(c, 'Matter Density ρ_m');

% Add a super title to the first figure
suptitle({'Chameleon Screening Mechanism:', 'V_{eff}(T) = Λ⁴exp(Λⁿ/Tⁿ) + λTρ_m'});

% --- PLOT 2 & 3: MASS & RANGE (New figure with subplots) ---
% Calculate effective masses
for i = 1:n_rho
    % Numerical second derivative around minimum
    T_min = T_minima(i);
    delta = 1e-4;

    % Calculate potential at T_min ± delta
    V_plus = Lambda^4 * exp((Lambda^n)/((T_min + delta)^n)) + lambda * (T_min + delta) * rho_m(i);
    V_min_val = Lambda^4 * exp((Lambda^n)/(T_min^n)) + lambda * T_min * rho_m(i);
    V_minus = Lambda^4 * exp((Lambda^n)/((T_min - delta)^n)) + lambda * (T_min - delta) * rho_m(i);

    % Second derivative (finite difference)
    m_eff_sq = (V_plus - 2*V_min_val + V_minus) / (delta^2);
    m_eff(i) = sqrt(abs(m_eff_sq));
end

figure('Position', [100, 100, 1000, 600]); % <--- New figure for the remaining two plots
sgtitle('Chameleon Field Properties vs. Density'); % <--- Overall title for this figure

% Plot 2: Effective mass vs matter density
subplot(1, 2, 1);
loglog(rho_m, m_eff, 'b-o', 'LineWidth', 2, 'MarkerSize', 4, 'MarkerFaceColor', 'b');
xlabel('Matter Density ρ_m');
ylabel('Effective Mass m_{eff}');
title('Chameleon Mass vs Environmental Density');
grid on;

% Add reference lines for different environments
hold on;
rho_galaxy = 1e-24;      % Galactic density
rho_solar = 1e-18;       % Solar system density
rho_lab = 1e-12;         % Laboratory density

yL = ylim;
plot([rho_galaxy rho_galaxy], yL, 'r--', 'LineWidth', 1);
plot([rho_solar rho_solar], yL, 'g--', 'LineWidth', 1);
plot([rho_lab rho_lab], yL, 'm--', 'LineWidth', 1);

text(rho_galaxy, m_eff(1)*10, ' Galactic', 'Color', 'r', 'FontSize', 8);
text(rho_solar, m_eff(1)*10, ' Solar System', 'Color', 'g', 'FontSize', 8);
text(rho_lab, m_eff(1)*10, ' Laboratory', 'Color', 'm', 'FontSize', 8);

% Plot 3: Compton wavelength (range of force
