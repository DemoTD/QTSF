% Retrocausal Suppression - Scaling Figures
% Figures 3 & 4: System size and entropy scaling

clear; clc; close all;

% Parameters
beta_retro = 0.1;      % Initial retrocausal coupling
gamma_fixed = 0.1;     % Fixed decoherence rate
time_fixed = 10;       % Fixed time

% Create figure for scaling behavior
figure('Position', [100, 100, 1000, 800]);

% Plot 1: Exponential suppression vs environment size
subplot(2, 1, 1);
% Corrected N_range: ending at 10^4 for better visualization
N_range = logspace(0, 4, 100);  % <--- Corrected range from 1 to 10^4
suppression = exp(-gamma_fixed * time_fixed * N_range);

loglog(N_range, suppression, 'b-', 'LineWidth', 2.5);
hold on;
loglog(N_range, beta_retro * suppression, 'r--', 'LineWidth', 2.5);

xlabel('Environmental Degrees of Freedom (N)');
ylabel('Magnitude');
title('(c) Exponential Suppression with System Size');
legend('Decoherence Factor', 'β_{eff} = β × Decoherence', 'Location', 'southwest');
grid on;

% Add markers and labels for specific system scales
micro_N = 1;
meso_N = 1e3;
% The macroscopic scale (1e23) is now outside this plot for clarity
% To show its effect, the second plot with entropy is more suitable

plot(micro_N, exp(-gamma_fixed * time_fixed * micro_N), 'ko', ...
    'MarkerSize', 10, 'MarkerFaceColor', 'k');
plot(meso_N, exp(-gamma_fixed * time_fixed * meso_N), 'ko', ...
    'MarkerSize', 10, 'MarkerFaceColor', 'k');

text(micro_N*1.5, 1e-1, 'Microscopic', 'FontSize', 10, 'FontWeight', 'bold');
text(meso_N*0.8, 1e-10, 'Mesoscopic', 'FontSize', 10, 'FontWeight', 'bold');
% Macroscopic label removed as it's not in the visible range

set(gca, 'FontSize', 12);

% Plot 2: Quantum vs Classical regime
subplot(2, 1, 2);
S_range = linspace(0, 15, 100);  % Entropy range
beta_suppressed = beta_retro * exp(-S_range);

semilogy(S_range, beta_suppressed, 'm-', 'LineWidth', 2.5);
xlabel('Entropy S (or log(N))');
ylabel('Effective Retrocausal Coupling β_{eff}');
title('(d) Suppression vs Entropy/Complexity');
grid on;
set(gca, 'FontSize', 12);

% Add regime labels and shading
quantum_cutoff = 3;
classical_start = find(S_range > quantum_cutoff, 1);
if ~isempty(classical_start)
    hold on;
    % Quantum regime shading
    fill([S_range(1), S_range(classical_start), S_range(classical_start), S_range(1)], ...
          [1e-20, 1e-20, 1, 1], [0.8, 0.9, 1], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    % Classical regime shading
    fill([S_range(classical_start), S_range(end), S_range(end), S_range(classical_start)], ...
          [1e-20, 1e-20, 1, 1], [1, 0.9, 0.8], 'FaceAlpha', 0.3, 'EdgeColor', 'none');

    text(mean(S_range(1:classical_start)), beta_retro/5, 'Quantum Regime', ...
         'FontSize', 11, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    text(mean(S_range(classical_start:end)), beta_retro/500, 'Classical Regime', ...
         'FontSize', 11, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');

    % Add vertical line at transition
    plot([quantum_cutoff, quantum_cutoff], [1e-20, 1], 'k--', 'LineWidth', 1.5);
    text(quantum_cutoff+0.2, beta_retro/50, 'Quantum-Classical Transition', ...
         'FontSize', 9, 'Rotation', 90);
end

% Save as high-resolution PNG
exportgraphics(gcf, 'retrocausal_scaling.png', 'Resolution', 300);

% Display key results
fprintf('Retrocausal Suppression Scaling Results:\n');
fprintf('Initial retrocausal coupling: β = %.3f\n', beta_retro);
fprintf('\nSuppression at different scales:\n');
fprintf('Microscopic (N=1): β_eff ≈ %.3f\n', beta_retro * exp(-gamma_fixed * time_fixed * 1));
fprintf('Mesoscopic (N=1e3): β_eff ≈ %.3e\n', beta_retro * exp(-gamma_fixed * time_fixed * 1e3));
fprintf('Macroscopic (N=1e23): β_eff ≈ %.3e\n', beta_retro * exp(-gamma_fixed * time_fixed * 1e23));
fprintf('Suppression factor: %.1e\n', exp(-gamma_fixed * time_fixed * 1e23));
