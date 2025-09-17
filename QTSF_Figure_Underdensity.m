% Generate Figure A1: Required Underdensity to Resolve Hubble Tension
% Author: Tiziano Demaria
% Date: 11 September 2025

% Parameters from the paper
H0_DISCREPANCY = 0.08;     % Observed ΔH₀/H₀
CONSTANT_TERM = 0.19;      % Derived constant

% Create array of f values
f = linspace(0.1, 1.0, 100);

% Calculate required underdensity: δ ≈ 0.19 / f
required_underdensity = CONSTANT_TERM ./ f;

% Create figure
figure('Position', [100, 100, 800, 500]);

% Plot the main relationship
plot(f, required_underdensity, 'Color', [65/255, 105/255, 225/255], ...
    'LineWidth', 3, 'DisplayName', 'Model: |δ∇T| = 0.19/f');
hold on;

% Highlight specific values
scatter([1.0, 0.5], [0.19, 0.38], 100, 'crimson', 'filled', ...
    'DisplayName', 'Specific solutions');
text(1.02, 0.17, 'f=1.0, δ=-19%', 'FontSize', 11, 'HorizontalAlignment', 'left');
text(0.52, 0.36, 'f=0.5, δ=-38%', 'FontSize', 11, 'HorizontalAlignment', 'left');

% Add horizontal line for observed Hubble tension
yline(H0_DISCREPANCY, '--', 'Observed ΔH₀/H₀ = 0.08', ...
    'Color', 'green', 'Alpha', 0.7, 'FontSize', 10, 'DisplayName', '');

% Customize plot
xlabel('Fraction of Dark Matter from T-field Gradients (f)', 'FontSize', 12);
ylabel('Required Underdensity (|δ∇T|)', 'FontSize', 12);
title('Required Underdensity to Resolve Hubble Tension', 'FontSize', 14);
grid on;
legend('Location', 'northeast');

% Axis limits and ticks
xlim([0.1, 1.0]);
ylim([0.0, 0.5]);

xticks(0.1:0.1:1.0);

% Add explanatory text box
textstr = sprintf(['Model prediction: |δ∇T| = 0.19/f\n\n' ...
    'For the observed Hubble tension\n' ...
    'ΔH₀/H₀ ≈ 0.08, the required underdensity\n' ...
    'depends on the fraction f of dark matter\n' ...
    'constituted by T-field gradients.']);
annotation('textbox', [0.15, 0.55, 0.3, 0.3], 'String', textstr, ...
    'FontSize', 11, 'BackgroundColor', [1, 0.9, 0.6], 'EdgeColor', 'none', ...
    'VerticalAlignment', 'top');

% Save the figure (optional)
% print('required_underdensity_plot','-dpng','-r300');

% Print derivation
fprintf('Derivation of the relationship:\n');
fprintf('From Appendix A1: ΔH₀/H₀ ≈ 0.42 · f · δ∇T\n');
fprintf('With observed ΔH₀/H₀ = %.2f:\n', H0_DISCREPANCY);
fprintf('%.2f ≈ 0.42 · f · δ∇T\n', H0_DISCREPANCY);
fprintf('Therefore: δ∇T ≈ 0.08 / (0.42 · f) ≈ 0.19 / f\n');

