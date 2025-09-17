% K-essence Dynamics and Effective Equation of State Simulation
% For Lagrangian: P(X) = α(√(2X) - M)²
% Shows w = p/ρ = (Y - M)/(Y + M) where Y = √(2X)

clear; clc; close all;

% Parameters
alpha = 1;      % Arbitrary scaling parameter (doesn't affect w)
M = 1;          % Mass scale parameter
Y_min = 0.1;    % Minimum Y value
Y_max = 5;      % Maximum Y value
n_points = 1000; % Number of points

% Create Y array
Y = linspace(Y_min, Y_max, n_points);

% Calculate energy density ρ and pressure p
rho = alpha * (Y - M).^2 .* (3 - (Y - M)./Y);
p = alpha * (Y - M).^2;

% Calculate equation of state parameter w
w = p ./ rho;

% Alternative direct calculation: w = (Y - M)/(Y + M)
w_direct = (Y - M) ./ (Y + M);

% Create a large figure for a single plot
figure('Position', [100, 100, 1000, 700]);

% Plot a single, large graph
plot(Y, w, 'b-', 'LineWidth', 2);
hold on;
plot(Y, w_direct, 'r--', 'LineWidth', 1.5);
xlabel('Y = √(2X)', 'FontSize', 12);
ylabel('Equation of State Parameter w', 'FontSize', 12);
title('K-essence Equation of State: w = p/ρ', 'FontSize', 14, 'FontWeight', 'bold');
legend('w = p/ρ', 'w = (Y-M)/(Y+M)', 'Location', 'southeast', 'FontSize', 10);
grid on;

% Add reference lines
yline(0, 'k--', 'w = 0 (Dust/Matter)');
yline(1, 'k--', 'w = 1 (Stiff Fluid)');
xline(M, 'k--', 'Y = M');

% Display key results
fprintf('K-essence Dynamics Simulation Results:\n');
fprintf('Lagrangian: P(X) = α(√(2X) - M)²\n');
fprintf('Equation of State: w = (Y - M)/(Y + M)\n\n');
fprintf('At Y = M: w = %.4f (Dark Matter regime)\n', (M - M)/(M + M));
fprintf('At Y = 5M: w = %.4f (approaching stiff fluid)\n', (5*M - M)/(5*M + M));
fprintf('At Y = 0.1M: w = %.4f\n', (0.1*M - M)/(0.1*M + M));
