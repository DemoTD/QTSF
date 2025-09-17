% Retrocausal Suppression via Decoherence Simulation
% Toy model demonstrating how environmental decoherence suppresses retrocausal effects

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

% Create figure
figure('Position', [100, 100, 1200, 800]);

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
subplot(2, 2, 1);
hold on;
for g_idx = 1:length(gamma_range)
    plot(t, beta_eff(g_idx, :), 'LineWidth', 2, ...
        'Color', colors(g_idx, :), ...
        'DisplayName', sprintf('γ = %.2f', gamma_range(g_idx)));
end
xlabel('Time');
ylabel('Effective Retrocausal Coupling β_{eff}');
title('Suppression of Retrocausal Effects by Decoherence');
legend('Location', 'northeast');
grid on;
set(gca, 'YScale', 'log');

% Plot 2: Purity decay over time
subplot(2, 2, 2);
hold on;
for g_idx = 1:length(gamma_range)
    plot(t, purity(g_idx, :), 'LineWidth', 2, ...
        'Color', colors(g_idx, :), ...
        'DisplayName', sprintf('γ = %.2f', gamma_range(g_idx)));
end
xlabel('Time');
ylabel('Purity Tr(ρ²)');
title('Decoherence Process: Loss of Quantum Coherence');
legend('Location', 'southwest');
grid on;

% Plot 3: Exponential suppression vs environment size
subplot(2, 2, 3);
N_range = 1:100;  % Range of environmental degrees of freedom
gamma_fixed = 0.1;
time_fixed = 10;

suppression = exp(-gamma_fixed * time_fixed * N_range);
semilogy(N_range, suppression, 'b-', 'LineWidth', 2);
hold on;
semilogy(N_range, beta_retro * suppression, 'r--', 'LineWidth', 2);

xlabel('Environmental Degrees of Freedom (N)');
ylabel('Suppression Factor');
title('Exponential Suppression with System Size');
legend('Decoherence Factor', 'β_{eff} = β × Decoherence', 'Location', 'southwest');
grid on;

% Add markers for specific system sizes
macroscopic_N = 1e23;  % Avogadro-scale system
marker_N = [1, 10, 100, macroscopic_N];
marker_suppression = exp(-gamma_fixed * time_fixed * marker_N);
loglog(marker_N, marker_suppression, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k');
text(marker_N(3)*1.1, marker_suppression(3)*0.5, 'Macroscopic', 'FontSize', 8);

% Plot 4: Quantum vs Classical regime
subplot(2, 2, 4);
S_range = linspace(0, 10, 100);  % Entropy range
beta_suppressed = beta_retro * exp(-S_range);

plot(S_range, beta_suppressed, 'm-', 'LineWidth', 2);
xlabel('Entropy S (or log(N))');
ylabel('Effective Retrocausal Coupling β_{eff}');
title('Suppression vs Entropy/Complexity');
grid on;
set(gca, 'YScale', 'log');

% Add regime labels
quantum_cutoff = 1;
classical_start = find(S_range > quantum_cutoff, 1);
if ~isempty(classical_start)
    hold on;
    area(S_range(1:classical_start), beta_suppressed(1:classical_start), ...
        'FaceColor', 'c', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    area(S_range(classical_start:end), beta_suppressed(classical_start:end), ...
        'FaceColor', 'y', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    text(mean(S_range(1:classical_start)), beta_retro/2, 'Quantum Regime', ...
        'HorizontalAlignment', 'center');
    text(mean(S_range(classical_start:end)), beta_retro/100, 'Classical Regime', ...
        'HorizontalAlignment', 'center');
end

% Add overall title
sgtitle('Retrocausal Suppression via Environmental Decoherence', ...
        'FontSize', 14, 'FontWeight', 'bold');

% Display key results
fprintf('Retrocausal Suppension via Decoherence - Simulation Results:\n');
fprintf('Initial retrocausal coupling: β = %.3f\n', beta_retro);
fprintf('\nSuppression at different scales:\n');
fprintf('Microscopic (N=1): β_eff ≈ %.3f\n', beta_retro * exp(-gamma_fixed * time_fixed * 1));
fprintf('Mesoscopic (N=1e3): β_eff ≈ %.3e\n', beta_retro * exp(-gamma_fixed * time_fixed * 1e3));
fprintf('Macroscopic (N=1e23): β_eff ≈ %.3e\n', beta_retro * exp(-gamma_fixed * time_fixed * 1e23));
fprintf('Suppression factor: %.1e\n', exp(-gamma_fixed * time_fixed * 1e23));
