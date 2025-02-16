% Implement a convergence analysis for mimetic difference approximation
% to the Darcy problem with unit permeability, using a trigonometric
% manufactured solution

ncells = [16; 32; 64; 128];

errors_p_cell_centers = zeros(length(ncells), 1);
errors_p_edge_centers = zeros(length(ncells), 1);
errors_p_nodes = zeros(length(ncells), 1);
errors_q_edge_centers = zeros(length(ncells), 1);
errors_p_q1 = zeros(length(ncells), 1);

for ii=1:length(ncells)
    errors = darcy_unit_perm_with_q1(ncells(ii), 2);
    errors_p_cell_centers(ii) = errors.p_cell_centers;
    errors_p_edge_centers(ii) = errors.p_edge_centers;
    errors_p_nodes(ii) = errors.p_nodes;
    errors_q_edge_centers(ii) = errors.q_edge_centers;
    errors_p_q1(ii) = errors.p_q1_L2;
end

% Compute log values for plotting
log_ncells = log2(ncells);
log_errors_p_cell_centers = log2(errors_p_cell_centers);
log_errors_p_edge_centers = log2(errors_p_edge_centers);
log_errors_p_nodes = log2(errors_p_nodes);
log_errors_q_edge_centers = log2(errors_q_edge_centers);
log_errors_p_q1 = log2(errors_p_q1);

% Plot errors
figure();
plot(log_ncells, log_errors_p_cell_centers, '-or', 'LineWidth', 1.5);
hold on;
plot(log_ncells, log_errors_p_edge_centers, '-ob', 'LineWidth', 1.5);
plot(log_ncells, log_errors_p_nodes, '-og', 'LineWidth', 1.5);
plot(log_ncells, log_errors_q_edge_centers, '-om', 'LineWidth', 1.5);
plot(log_ncells, log_errors_p_q1, '-ok', 'LineWidth', 1.5);

% -----> Choose your custom first y-point for the second-order reference slope
y_start = -10;  % Change this to any desired starting y-point

% Define the second-order slope (-2) while keeping the full x-range
x_ref = [log_ncells(1), log_ncells(end)];
y_ref = [y_start, y_start - 2 * (log_ncells(end) - log_ncells(1))];

plot(x_ref, y_ref, '--k', 'LineWidth', 1.5); % Dashed black line for order 2 slope

% Labels and legend
xlabel('log2(Number of cells)', 'FontSize', 15);
ylabel('log2(Discrete relative L2-error)', 'FontSize', 15);

legend({'Pressure (cell-center)', ...
        'Pressure (edge-center)', ...
        'Pressure (node)', ...
        'Flux (edge-center)', ...
        'Pressures (Q1)', ...
        'Reference slope 2'}, 'Location', 'Best', 'FontSize', 15);

grid on;
hold off;
print('convergence.png', '-depsc');
