function test_q1_reconstruction_nodal_values(X_nodes, Y_nodes, U_nodes)
    % Compute Q1 reconstruction
    A = q1_reconstruction(X_nodes, Y_nodes, U_nodes);

    Nx = size(X_nodes,1) - 1; % Number of cells in x direction
    Ny = size(Y_nodes,2) - 1; % Number of cells in y direction

    cell_idx = 1; % Counter for grid cell index
    max_error = 0; % Track maximum error

    for i = 1:Nx
        for j = 1:Ny
            % Extract nodal coordinates (corners)
            x1 = X_nodes(i,j);   y1 = Y_nodes(i,j);
            x3 = X_nodes(i+1,j); y3 = Y_nodes(i+1,j);
            x7 = X_nodes(i,j+1); y7 = Y_nodes(i,j+1);
            x9 = X_nodes(i+1,j+1); y9 = Y_nodes(i+1,j+1);

            % Extract function values (corners)
            u1 = U_nodes(i,j);
            u3 = U_nodes(i+1,j);
            u7 = U_nodes(i,j+1);
            u9 = U_nodes(i+1,j+1);

            % Reconstructed coefficients for this cell
            a = A(cell_idx, :);

            % Evaluate reconstructed function at nodal points
            u1_rec = a(1) + a(2)*x1 + a(3)*y1 + a(4)*x1*y1;
            u3_rec = a(1) + a(2)*x3 + a(3)*y3 + a(4)*x3*y3;
            u7_rec = a(1) + a(2)*x7 + a(3)*y7 + a(4)*x7*y7;
            u9_rec = a(1) + a(2)*x9 + a(3)*y9 + a(4)*x9*y9;

            % Compute absolute errors
            err1 = abs(u1_rec - u1);
            err3 = abs(u3_rec - u3);
            err7 = abs(u7_rec - u7);
            err9 = abs(u9_rec - u9);

            % Update max error for reporting
            max_error = max([max_error, err1, err3, err7, err9]);

            % Display warnings if errors are too large
            if max([err1, err3, err7, err9]) > 1e-10
                fprintf('Warning: Large reconstruction error at cell (%d, %d)\n', i, j);
                fprintf('Errors: [%e, %e, %e, %e]\n', err1, err3, err7, err9);
            end

            cell_idx = cell_idx + 1;
        end
    end

    fprintf('Max reconstruction error: %e\n', max_error);

    if max_error < 1e-12
        fprintf('Test PASSED: Q1 reconstruction is exact at nodal points.\n');
    else
        fprintf('Test FAILED: Significant errors detected in nodal recovery.\n');
    end
end

