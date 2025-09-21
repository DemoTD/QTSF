% Numerical Solution for Structure Growth with Different Equations of State
% Author: Tiziano Demaria
% Date: 11 September 2025

% Define Hubble parameter: H(a)^2 / H0^2 = Om/a^3 + Ol
H2_over_H02 = @(a, Om) Om ./ a.^3 + (1 - Om);

% Define the growth equation for a component with equation of state w
function dydx = growth_equation(a, y, w, Om)
    delta = y(1);
    ddelta_da = y(2);

    H2 = Om / a^3 + (1 - Om);
    dlnH_da = (-3 * Om / a^4) / (2 * H2);
    coefficient = (3/2) * Om * (1 - w) * (1 + 3*w) / (a^5 * H2);

    d2delta_da2 = - (dlnH_da + 3 / a) * ddelta_da + coefficient * delta;
    dydx = [ddelta_da; d2delta_da2];
end

% Integration parameters
a_initial = 0.01;
a_final = 1.0;
a_values = linspace(a_initial, a_final, 1000);

% Initial conditions: delta = 1, ddelta/da = 0
y_initial = [1; 0];

% Integrate for w = 0 (CDM)
w_cdm = 0.0;
Om = 0.3;
solution_cdm = zeros(length(a_values), 2);
solution_cdm(1, :) = y_initial';

for i = 2:length(a_values)
    a_span = [a_values(i-1), a_values(i)];
    [~, y_out] = ode45(@(a, y) growth_equation(a, y, w_cdm, Om), a_span, solution_cdm(i-1, :)');
    solution_cdm(i, :) = y_out(end, :);
end

delta_cdm = solution_cdm(:, 1);
D_cdm = delta_cdm / delta_cdm(1);

% Integrate for w = -1/3 (T-field gradient energy)
w_tfield = -1/3;
solution_w = zeros(length(a_values), 2);
solution_w(1, :) = y_initial';

for i = 2:length(a_values)
    a_span = [a_values(i-1), a_values(i)];
    [~, y_out] = ode45(@(a, y) growth_equation(a, y, w_tfield, Om), a_span, solution_w(i-1, :)');
    solution_w(i, :) = y_out(end, :);
end

delta_w = solution_w(:, 1);
D_w = delta_w / delta_w(1);

% Plot
figure;
plot(a_values, D_cdm, 'r--', 'LineWidth', 2);
hold on;
plot(a_values, D_w, 'b-', 'LineWidth', 2);
xlabel('Scale Factor (a)');
ylabel('Growth Factor D(a) = \delta(a)/\delta(a_{initial})');
title('Growth of Structure: Comparison of Equations of State');
legend('CDM (w = 0)', 'T-field Gradients (w = -1/3)');
grid on;
print('growth_comparison','-dpng');

