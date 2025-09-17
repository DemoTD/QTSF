% Retrocausal Suppression - Time Evolution Figures
% Figures 1 & 2: Time evolution of suppression and purity

clear; clc; close all;

% Parameters
gamma_range = linspace(0.01, 0.2, 6);  % Decoherence rates
N_env = 50;                            % Number of environmental degrees of freedom
t_max = 50;                            % Maximum time
dt = 0.1;                              % Time step
t = 0:dt:t_max;

% Retrocausal coupling strength (small off-diagonal term)
beta_retro = 0.1;

% Initialize results
beta_eff = zeros(length(gamma_range), length(t));
purity = zeros(length(gamma_range), length(t));

% Colors for plotting
colors = lines(length(gamma_range));

% Create figure for time evolution
figure('Position', [100, 100, 1000, 800]);

% Main simulation loop
for g_idx = 1:length(gamma_range)
    gamma = gamma_range(g_idx);

    % Initial state: superposition state
    psi0 = [1; 1]/sqrt(2);  % |+> state
    rho0 = psi0 * psi0';

    rho = rho0;

    for t_idx = 1:length(t)
        % Simple decoherence model: off-diagonal damping
        current_time = t(t_idx);

        % Decoherence effect: exponential decay of off-diagonals
        decoherence_factor = exp(-gamma * current_time * N_env);

        % Effective retrocausal term (suppressed by decoherence)
        beta_eff(g_idx, t_idx) = beta_retro * decoherence_factor;

        % Update density matrix (simplified model)
        rho(1,2) = rho0(1,2) * decoherence_factor;
        rho(2,1) = rho0(2,1) * decoherence_factor;
        rho(1,1) = 0.5 * (1 + (1 - decoherence_factor));
        rho(2,2) = 0.5 * (1 - (1 - decoherence_factor));

        % Calculate purity
        purity(g_idx, t_idx) = trace(rho^2);
    end
end

% Plot 1: Suppression of retrocausal term over time
subplot(2, 1, 1);
hold on;
for g_idx = 1:length(gamma_range)
    plot(t, beta_eff(g_idx, :), 'LineWidth', 2.5, ...
        'Color', colors(g_idx, :), ...
        'DisplayName', sprintf('γ = %.2f', gamma_range(g_idx)));
end
xlabel('Time');
ylabel('Effective Retrocausal Coupling β_{eff}');
title('(a) Time Evolution of Retrocausal Suppression');
legend('Location', 'northeast');
grid on;
set(gca, 'YScale', 'log', 'FontSize', 12);

% Plot 2: Purity decay over time
subplot(2, 1, 2);
hold on;
for g_idx = 1:length(gamma_range)
    plot(t, purity(g_idx, :), 'LineWidth', 2.5, ...
        'Color', colors(g_idx, :), ...
        'DisplayName', sprintf('γ = %.2f', gamma_range(g_idx)));
end
xlabel('Time');
ylabel('Purity Tr(ρ²)');
title('(b) Decoherence Process: Loss of Quantum Coherence');
legend('Location', 'southwest');
grid on;
set(gca, 'FontSize', 12);

% Save as high-resolution PNG
exportgraphics(gcf, 'retrocausal_time_evolution.png', 'Resolution', 300);
