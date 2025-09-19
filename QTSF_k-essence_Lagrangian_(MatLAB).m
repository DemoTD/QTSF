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

% Create figure
figure('Position', [100, 100, 1200, 800]);

% Main plot: Equation of state w vs Y
subplot(2, 2, 1);
plot(Y, w, 'b-', 'LineWidth', 2);
hold on;
plot(Y, w_direct, 'r--', 'LineWidth', 1.5);
xlabel('Y = √(2X)');
ylabel('Equation of State Parameter w');
title('K-essence Equation of State: w = p/ρ');
legend('w = p/ρ', 'w = (Y-M)/(Y+M)', 'Location', 'southeast');
grid on;

% Add reference lines
yline(0, 'k--', 'w = 0 (Dust/Matter)', 'LabelVerticalAlignment', 'middle');
yline(1, 'k--', 'w = 1 (Stiff Fluid)', 'LabelVerticalAlignment', 'middle');
xline(M, 'k--', 'Y = M', 'LabelVerticalAlignment', 'bottom');

% Energy density ρ vs Y
subplot(2, 2, 2);
semilogy(Y, rho, 'r-', 'LineWidth', 2);
xlabel('Y = √(2X)');
ylabel('Energy Density ρ');
title('Energy Density vs Y');
grid on;
xline(M, 'k--', 'Y = M');

% Pressure p vs Y
subplot(2, 2, 3);
plot(Y, p, 'g-', 'LineWidth', 2);
xlabel('Y = √(2X)');
ylabel('Pressure p');
title('Pressure vs Y');
grid on;
xline(M, 'k--', 'Y = M');

% Zoomed view near Y = M
subplot(2, 2, 4);
Y_zoom = linspace(M*0.8, M*1.2, 200);
w_zoom = (Y_zoom - M) ./ (Y_zoom + M);
plot(Y_zoom, w_zoom, 'm-', 'LineWidth', 2);
xlabel('Y = √(2X)');
ylabel('Equation of State Parameter w');
title('Zoomed View Near Y = M');
grid on;
xline(M, 'k--', 'Y = M');
yline(0, 'k--', 'w = 0');

% Add overall title and annotation
sgtitle('K-essence Dynamics: P(X) = α(√(2X) - M)²', 'FontSize', 14, 'FontWeight', 'bold');

% Display key results
fprintf('K-essence Dynamics Simulation Results:\n');
fprintf('Lagrangian: P(X) = α(√(2X) - M)²\n');
fprintf('Equation of State: w = (Y - M)/(Y + M)\n\n');
fprintf('At Y = M: w = %.4f (Dark Matter regime)\n', (M - M)/(M + M));
fprintf('At Y = 5M: w = %.4f (approaching stiff fluid)\n', (5*M - M)/(5*M + M));
fprintf('At Y = 0.1M: w = %.4f\n', (0.1*M - M)/(0.1*M + M));
