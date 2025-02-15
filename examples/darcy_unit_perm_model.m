function [error_p, error_p_nodes] = darcy_unit_perm_model(ncells, k)
    % Model implementation of Darcy problem with unit permeability

    % addpath('/Users/jvmini/Git/mole-master/src/mole_MATLAB')
    addpath('/Users/jvpro/Documents/GitHub/mole/src/MATLAB')

    % Input parameters
    % ncells = 40;  % number of cells in x and y directions
    % k = 2;  % degree of the mimetic operator

    % True solution
    p  = @(x, y) sin(pi*x) .* sin(pi*y);
    qx = @(x, y) -pi * cos(pi*x) .* sin(pi*y);
    qy = @(x, y) -pi * sin(pi*x) .* cos(pi*y);
    f  = @(x, y) 2*pi*pi*sin(pi*x) .* sin(pi*y);

    % Domain's limits
    west = 0; east = 1; south = 0; north = 1;

    % Number of cells and step size
    nx = ncells; dx = (east-west)/nx;
    ny = ncells; dy = (north-south)/ny;

    % 2D Mimetic laplacian operator
    L = lap2D(k, nx, dx, ny, dy);

    % Impose Robin BC on 2D Laplacian operator
    a = 1; b = 0;
    L = L + robinBC2D(k, nx, dx, ny, dy, a, b);

    % Construct the staggered grid
    x_grid = [west west+dx/2 : dx : east-dx/2 east];
    y_grid = [south south+dy/2 : dy : north-dy/2 north];
    [X, Y] = meshgrid(x_grid, y_grid);

    % Exact source term in each cell-center
    source = f(X, Y);

    % Assemble RHS and Enforce zero BC.
    RHS = source;
    RHS(1, :) = 0;
    RHS(end, :) = 0;
    RHS(:, 1) = 0;
    RHS(:, end) = 0;
    RHS = reshape(RHS, [], 1);

    % Solve linear system to obtain mimetic pressure
    p_mimetic = -L\RHS;

    % Retrieve mimetic flux solution
    G = grad2D(k, nx, dx, ny, dy);
    q_mimetic = - G * p_mimetic;

    % Now, we separate between horizontal and vertical edges
    %x_hori_edges = west + dx/2 : dx : east - dx/2;
    %y_hori_edges = south : dy : north;
    %[X_hori_edges, Y_hori_edges] = meshgrid(x_hori_edges, y_hori_edges);
    %q_hori_edges = q_mimetic(1:end/2);

    %x_vert_edges = west : dx : east;
    %y_vert_edges = south + dy/2 : dy : north - dy/2;
    %[X_vert_edges, Y_vert_edges] = meshgrid(x_vert_edges, y_vert_edges);
    %q_vert_edges = q_mimetic(end/2+1:end);

    % Exact pressure and flux values
    p_exact = p(X, Y);
    p_exact = p_exact(:);

    %q_vert_edges_exact = qx(X_vert_edges, Y_vert_edges);  % vert edges, x-component
    %q_vert_edges_exact = q_vert_edges_exact(:);
    %q_hori_edges_exact = qy(X_hori_edges, X_hori_edges);  % hori edges, y-component
    %q_hori_edges_exact = q_hori_edges_exact(:);
    %q_exact = [q_vert_edges_exact; q_hori_edges_exact];   % [qx, qy] --> makes sense

    % Error computation
    error_p = norm(p_mimetic - p_exact) / norm(p_exact);
    %error_q = norm(q_mimetic - q_exact) / norm(q_exact);
    %disp('Error Pressure'); disp(error_p);
    %disp('Error Flux'); disp(error_q);

    % Postprocessing

    # --> Center to faces interpolators
    %C2F = interpolCentersToFacesD2D(k, ncells, ncells);
    %p_faces = C2F * [p_mimetic; p_mimetic];
    %p_vert_edges = p_faces(1:end/2);
    %p_hori_edges = p_faces(end/2+1:end);

    # --> Center to nodes interpolators
    x_nodes = west:dx:east;
    y_nodes = south:dy:north;
    [X_nodes, Y_nodes] = meshgrid(x_nodes, y_nodes);
    C2N = interpolCentersToNodes2D(k, ncells, ncells);
    p_nodes_mimetic = C2N * p_mimetic;
    p_nodes_true = p(X_nodes, Y_nodes);
    error_p_nodes = norm(p_nodes_mimetic - p_nodes_true) / norm(p_nodes_true);



end








