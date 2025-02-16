function plot_q1_vs_mimetic(X_staggered, Y_staggered, U_num, X_nodes, Y_nodes, U_nodes)
    % Compute Q1 reconstruction coefficients
    A = q1_reconstruction(X_nodes, Y_nodes, U_nodes);

    Nx = size(X_nodes,1) - 1; % Number of cells in x direction
    Ny = size(Y_nodes,2) - 1; % Number of cells in y direction

    % Allocate reconstructed field
    U_q1 = zeros(size(U_num));

    % Evaluate Q1 reconstruction at staggered grid centers
    cell_idx = 1;
    for i = 1:Nx
        for j = 1:Ny
            % Get staggered grid cell center
            xc = X_staggered(i+1, j+1);
            yc = Y_staggered(i+1, j+1);

            % Reconstructed Q1 coefficients
            a = A(cell_idx, :);

            % Evaluate Q1 function at staggered cell center
            U_q1(i+1, j+1) = a(1) + a(2)*xc + a(3)*yc + a(4)*xc*yc;

            cell_idx = cell_idx + 1;
        end
    end

    % Create figure
    figure;

    % Plot mimetic solution
    subplot(1,2,1);
    surf(X_staggered, Y_staggered, U_num, 'EdgeColor', 'none');
    title('Mimetic Solution (U_{num})');
    xlabel('x');
    ylabel('y');
    zlabel('Potential');
    colorbar;
    view(3);

    % Plot reconstructed Q1 solution
    subplot(1,2,2);
    surf(X_staggered, Y_staggered, U_q1, 'EdgeColor', 'none');
    title('Reconstructed Q1 Solution (U_{Q_1})');
    xlabel('x');
    ylabel('y');
    zlabel('Potential');
    colorbar;
    view(3);

    % Display max error
    max_error = max(abs(U_q1(:) - U_num(:)));
    fprintf('Max reconstruction error: %e\n', max_error);
end

