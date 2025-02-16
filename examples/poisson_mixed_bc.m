% This code implements the Poisson Problem with mixedBC
% Left bc   -> u = 1
% Right bc  -> u = 0
% Bottom bc -> u' = 0
% Top bc    -> u' = 0

% Clear and add path
clear all; clc; close all;
%addpath('/Users/jvmini/Git/mole-master/src/mole_MATLAB')
addpath('/Users/jvpro/Documents/GitHub/mole/src/matlab')

% Set domain's and create relevant grids
domain = [0, 1, 0, 1]; % xmin, xmax, ymin, ymax
nx = 8; ny = 10;
[dx, dy] = step_size_2d([nx, ny], domain);
[StagX, StagY] = staggered_grid_2d([dx, dy], domain);
[HoriEdgesX, HoriEdgesY] = horizontal_edges_grid_2d([dx, dy], domain);
[VertEdgesX, VertEdgesY] = vertical_edges_grid_2d([dx, dy], domain);
[NodesX, NodesY] = nodes_grid_2d([dx, dy], domain);

%  Construct Laplacian operator
k = 2;
L = lap2D(k, nx, dx, ny, dy);

% Set up mixed_bc boundary
mixed_bc = mixedBC2D(k,  ...
                     nx,  dx, ...
                     ny,  dy, ...
                     "Dirichlet", 1, ...  % Dirichlet on the left side
                     "Dirichlet", 1, ...  % Dirichlet on the right side
                     "Neumann", 1, ...  % Neumann on the bottom side
                     "Neumann", 1);     % Neumann on the top side

% Update Laplacian to account for boundary conditions
L = L + mixed_bc;

% Assemble RHS and solve

b = zeros(nx+2, ny+2)';
b(2:end-1, 1) = 1;  % known value at the left boundary
% Note that we do not include the corner points here, otherwise we get
% the wrong values.

b = reshape(b', [], 1);

p = L \ b;

% Plot Numerical Pressure

% Plot
figure();
surf(StagX, StagY, reshape(p, nx+2, ny+2)');
title(['Numerical pressure, k=', num2str(k)]);
xlabel('x'); ylabel('y'); colorbar;

% Retrieve flux solution
G = grad2D(k, nx, dx, ny, dy);
q = - G * p;
qx = q(1:((nx+1)*ny));
qy = q((nx+1)*ny+1:end);

% Plot horizontal and vertical fluxes
figure();
surf(HoriEdgesX, HoriEdgesY, reshape(qy, nx, ny+1)');
title('Vertical Flux'); colorbar;

figure();
surf(VertEdgesX, VertEdgesY, reshape(qx, nx+1, ny)');
title('Horizontal Flux'); colorbar;

% Nodal interpolation
p_nodes = cell_centered_pressure_to_nodes(k, [nx, ny], p);
figure();
surf(NodesX, NodesY, reshape(p_nodes, nx+1, ny+1)');
title('Nodal pressures')

% Face pressures
[p_hori_edges, p_vert_edges] = cell_centered_pressure_to_edges(k, [nx, ny], p);

figure();
surf(HoriEdgesX, HoriEdgesY, reshape(p_hori_edges, nx, ny+1)');
title('Horizontal Edges Pressure');

figure();
surf(VertEdgesX, VertEdgesY, reshape(p_vert_edges, nx+1, ny)');
title('Vertical Edges Pressure');






