function [errors] = darcy_unit_perm_model(num_cells, k)
    % Model implementation of Darcy problem with unit permeability

    addpath('/Users/jvmini/Git/mole-master/src/mole_MATLAB')
    % addpath('/Users/jvpro/Documents/GitHub/mole/src/MATLAB')

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

    % 2D Mimetic laplacian operator
    L = lap2D(k, nx, dx, ny, dy);

    % Impose Robin BC on 2D Laplacian operator
    a = 1; b = 0;
    L = L + robinBC2D(k, nx, dx, ny, dy, a, b);

    % Exact source term at each cell-center
    source = f(StagX, StagY);

    % Assemble RHS and enforce zero BC.
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

    % Compute errors
    p_ex = p(StagX, StagY)';
    p_ex = p_ex(:);
    error_p_cc = norm(p_num - p_ex) ./ norm(p_ex);

    q_ex_x = qx(VertEdgesX, VertEdgesY)';
    q_ex_x = q_ex_x(:);
    q_ex_y = qy(HoriEdgesX, HoriEdgesY)';
    q_ex_y = q_ex_y(:);
    q_ex = [q_ex_x; q_ex_y];
    error_q = norm(q_num - q_ex) ./ norm(q_ex);

    p_nodes_ex = p(NodesX, NodesY)';
    p_nodes_ex = p_nodes_ex(:);
    error_p_nodes = norm(p_nodes - p_nodes_ex) ./ norm(p_nodes_ex);

    p_edges_ex_x = p(VertEdgesX, VertEdgesY)';
    p_edges_ex_x = p_edges_ex_x(:);
    p_edges_ex_y = p(HoriEdgesX, HoriEdgesY)';
    p_edges_ex_y = p_edges_ex_y(:);
    p_edges_ex = [p_edges_ex_x; p_edges_ex_y];
    error_p_edges = norm(p_edges_ex - p_edges) ./ norm(p_edges_ex);

    % Create struct
    errors = struct(...
        'p_cell_centers', error_p_cc, ...
        'p_nodes', error_p_nodes, ...
        'p_edge_centers', error_p_edges, ...
        'q_edge_centers', error_q);

end








