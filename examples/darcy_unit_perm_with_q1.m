function [errors] = darcy_unit_perm_with_q1(num_cells, k)
    % Model implementation of Darcy problem with unit permeability
    % Includes Q1 reconstruction and numerical quadrature-based error computation

    addpath('/Users/jvmini/Git/mole-master/src/mole_MATLAB')

    % Domain
    domain = [0, 1, 0, 1];

    % Construct grids
    nx = num_cells; ny = num_cells;
    [dx, dy] = step_size_2d([nx, ny], domain);
    [StagX, StagY] = staggered_grid_2d([dx, dy], domain);
    [HoriEdgesX, HoriEdgesY] = horizontal_edges_grid_2d([dx, dy], domain);
    [VertEdgesX, VertEdgesY] = vertical_edges_grid_2d([dx, dy], domain);
    [NodesX, NodesY] = nodes_grid_2d([dx, dy], domain);

    % Exact solution
    p  = @(x, y) sin(pi*x) .* sin(pi*y);
    qx = @(x, y) -pi * cos(pi*x) .* sin(pi*y);
    qy = @(x, y) -pi * sin(pi*x) .* cos(pi*y);
    f  = @(x, y) 2*pi*pi*sin(pi*x) .* sin(pi*y);

    % 2D Mimetic Laplacian operator
    L = lap2D(k, nx, dx, ny, dy);

    % Impose Robin BC on 2D Laplacian operator
    a = 1; b = 0;
    L = L + robinBC2D(k, nx, dx, ny, dy, a, b);

    % Exact source term at each cell-center
    source = f(StagX, StagY);

    % Assemble RHS and enforce zero BC
    RHS = source;
    RHS([1, end], :) = 0;
    RHS(:, [1, end]) = 0;
    RHS = reshape(RHS', [], 1);

    % Solve linear system to obtain mimetic pressure
    p_num = -L\RHS;

    % Retrieve flux solution
    G = grad2D(k, nx, dx, ny, dy);
    q_num = - G * p_num;
    qx_num = q_num(1:((nx+1)*ny));
    qy_num = q_num((nx+1)*ny+1:end);

    % Compute edge pressures
    [p_edges_y, p_edges_x] = cell_centered_pressure_to_edges(k, [nx, ny], p_num);
    p_edges = [p_edges_x; p_edges_y];

    % Compute nodal pressures
    p_nodes = cell_centered_pressure_to_nodes(k, [nx, ny], p_num);

    % ---------------- Q1 Reconstruction ---------------- %
    % Compute Q1 coefficients from nodal pressures
    A_q1 = q1_reconstruction(NodesX, NodesY, reshape(p_nodes, size(NodesX)));

    % ---------------- Numerical Quadrature for Error ---------------- %
    % Quadrature points (2x2 Gauss-Legendre)
    gp = [-1/sqrt(3), 1/sqrt(3)];  % Gauss points in reference element [-1,1]
    gw = [1, 1];  % Weights for 2-point quadrature

    % Compute L2 error using quadrature
    error_sum = 0;
    cell_idx = 1;

    for i = 1:nx
        for j = 1:ny
            % Get local nodal coordinates
            x1 = NodesX(i,j);   y1 = NodesY(i,j);
            x3 = NodesX(i+1,j); y3 = NodesY(i+1,j);
            x7 = NodesX(i,j+1); y7 = NodesY(i,j+1);
            x9 = NodesX(i+1,j+1); y9 = NodesY(i+1,j+1);

            % Compute Q1 coefficients for this cell
            a = A_q1(cell_idx, :);

            % Map Gauss points from reference to physical element
            for g1 = 1:2
                for g2 = 1:2
                    % Convert reference Gauss points to physical coordinates
                    xi = gp(g1); eta = gp(g2);
                    xc = 0.25 * ((1 - xi) * (1 - eta) * x1 + ...
                                 (1 + xi) * (1 - eta) * x3 + ...
                                 (1 - xi) * (1 + eta) * x7 + ...
                                 (1 + xi) * (1 + eta) * x9);
                    yc = 0.25 * ((1 - xi) * (1 - eta) * y1 + ...
                                 (1 + xi) * (1 - eta) * y3 + ...
                                 (1 - xi) * (1 + eta) * y7 + ...
                                 (1 + xi) * (1 + eta) * y9);

                    % Evaluate Q1 function at quadrature points
                    p_q1 = a(1) + a(2) * xc + a(3) * yc + a(4) * xc * yc;

                    % Exact solution at quadrature points
                    p_exact = p(xc, yc);

                    % Compute squared error
                    error_sum = error_sum + gw(g1) * gw(g2) * (p_q1 - p_exact)^2 * dx * dy;
                end
            end
            cell_idx = cell_idx + 1;
        end
    end

    % Compute L2 norm of the error
    error_p_q1 = sqrt(error_sum);

    % ---------------- Error Computation ---------------- %
    % Compute errors at cell centers (for comparison)
    p_ex = p(StagX, StagY)';
    p_ex = p_ex(:);
    error_p_cc = norm(p_num - p_ex) ./ norm(p_ex);

    % Compute errors at edges
    p_edges_ex_x = p(VertEdgesX, VertEdgesY)';
    p_edges_ex_x = p_edges_ex_x(:);
    p_edges_ex_y = p(HoriEdgesX, HoriEdgesY)';
    p_edges_ex_y = p_edges_ex_y(:);
    p_edges_ex = [p_edges_ex_x; p_edges_ex_y];
    error_p_edges = norm(p_edges_ex - p_edges) ./ norm(p_edges_ex);

    % Compute errors at nodes
    p_nodes_ex = p(NodesX, NodesY)';
    p_nodes_ex = p_nodes_ex(:);
    error_p_nodes = norm(p_nodes - p_nodes_ex) ./ norm(p_nodes_ex);

    % Compute flux errors
    q_ex_x = qx(VertEdgesX, VertEdgesY)';
    q_ex_x = q_ex_x(:);
    q_ex_y = qy(HoriEdgesX, HoriEdgesY)';
    q_ex_y = q_ex_y(:);
    q_ex = [q_ex_x; q_ex_y];
    error_q = norm(q_num - q_ex) ./ norm(q_ex);

    % ---------------- Store Errors in Struct ---------------- %
    errors = struct(...
        'p_cell_centers', error_p_cc, ...
        'p_q1_L2', error_p_q1, ...
        'p_nodes', error_p_nodes, ...
        'p_edge_centers', error_p_edges, ...
        'q_edge_centers', error_q);

end

