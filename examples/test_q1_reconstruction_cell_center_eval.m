function test_q1_reconstruction_cell_center_eval(X_nodes, Y_nodes, U_nodes, X_staggered, Y_staggered, U_num)
    % Verify if the reconstructed Q1 function correctly interpolates p_num at the staggered grid cell centers.
    %
    % INPUT:
    % X_nodes  - (Nx+1, Ny+1) matrix of x-coordinates for corner nodes
    % Y_nodes  - (Nx+1, Ny+1) matrix of y-coordinates for corner nodes
    % U_nodes  - (Nx+1, Ny+1) matrix of function values at corner nodes (interpolated)
    % X_staggered - (Nx+2, Ny+2) matrix of staggered grid x-coordinates
    % Y_staggered - (Nx+2, Ny+2) matrix of staggered grid y-coordinates
    % U_num    - (Nx+2, Ny+2) staggered grid pressure values (includes boundary info)

    % Compute Q1 reconstruction coefficients
    A = q1_reconstruction(X_nodes, Y_nodes, U_nodes);

    Nx = size(X_nodes,1) - 1; % Number of cells in x direction
    Ny = size(Y_nodes,2) - 1; % Number of cells in y direction

    cell_idx = 1; % Counter for grid cell index
    max_error = 0; % Track max error
    total_error = 0; % Sum error for averaging
    count = 0; % Counter for averaging

    for i = 1:Nx
        for j = 1:Ny
            % Staggered grid cell center
            xc_staggered = X_staggered(i+1, j+1);
            yc_staggered = Y_staggered(i+1, j+1);

            % Computed element center using Cartesian nodes
            xc_cartesian = 0.25 * (X_nodes(i,j) + X_nodes(i+1,j) + X_nodes(i,j+1) + X_nodes(i+1,j+1));
            yc_cartesian = 0.25 * (Y_nodes(i,j) + Y_nodes(i+1,j) + Y_nodes(i,j+1) + Y_nodes(i+1,j+1));

            % Compute absolute difference
            dx = abs(xc_cartesian - xc_staggered);
            dy = abs(yc_cartesian - yc_staggered);

            fprintf('Cell (%d, %d): dx = %e, dy = %e\n', i, j, dx, dy);
        end
    end


    for i = 1:Nx
        for j = 1:Ny
            % Get the staggered grid cell center coordinates
            xc = X_staggered(i+1, j+1);
            yc = Y_staggered(i+1, j+1);

            % Reconstructed Q1 coefficients
            a = A(cell_idx, :);

            % Evaluate Q1 function at the staggered grid cell center
            u_rec = a(1) + a(2)*xc + a(3)*yc + a(4)*xc*yc;

            % Reference mimetic solution at the staggered grid cell center
            u_exact = U_num(i+1, j+1);

            % Compute absolute error
            error = abs(u_rec - u_exact);

            % Track max and total error
            max_error = max(max_error, error);
            total_error = total_error + error;
            count = count + 1;

            % Display warnings for large errors
            if error > 1e-10
                fprintf('Warning: Large interpolation error at cell (%d, %d): %e\n', i, j, error);
            end

            cell_idx = cell_idx + 1;
        end
    end

    % Compute average error
    avg_error = total_error / count;

    % Report results
    fprintf('Max interpolation error at staggered cell centers: %e\n', max_error);
    fprintf('Average interpolation error at staggered cell centers: %e\n', avg_error);

    if max_error < 1e-12
        fprintf('Test PASSED: Q1 reconstruction is highly accurate at staggered grid cell centers.\n');
    else
        fprintf('Test WARNING: Significant errors detected in Q1 interpolation at staggered cell centers.\n');
    end
end

