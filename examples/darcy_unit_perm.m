graphics_toolkit("qt");

% This code implements the Darcy problem with unit permeability using
% a trigonometric manufactured solution satisfying zero-boundary conditions

% Clear and add path
clear all; clc; close all;
addpath('/Users/jvmini/Git/mole-master/src/mole_MATLAB')
addpath('/Users/jvmini/Git/mimetic_a_posteriori/src')

% addpath('/Users/jvpro/Documents/GitHub/mole/src/MATLAB')

% Input parameters
ncells = 20;  % number of cells in x and y directions
k = 2;  % degree of the mimetic operator

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

% Solve linear system
p_mimetic = -L\RHS;

% Retrieve flux solution
G = grad2D(k, nx, dx, ny, dy);
q_mimetic = - G * p_mimetic;

% Grid for vertical and horizontal fluxes
x_dual_grid_h = west + dx/2 : dx : east - dx/2;
y_dual_grid_h = south : dy : north;
[Xdualh, Ydualh] = meshgrid(x_dual_grid_h, y_dual_grid_h);
qy_mimetic = q_mimetic(1:length(q_mimetic)/2);

x_dual_grid_v = west : dx : east;
y_dual_grid_v = south + dy/2 : dy : north - dy/2;
[Xdualv, Ydualv] = meshgrid(x_dual_grid_v, y_dual_grid_v);
qx_mimetic = q_mimetic(length(q_mimetic)/2+1:end);

% Exact values pressure and flux values
p_exact = p(X, Y);
p_exact = p_exact(:);

qx_exact = qx(Xdualv, Ydualv);
qx_exact = qx_exact(:);
qy_exact = qy(Xdualh, Ydualh);
qy_exact = qy_exact(:);
q_exact = [qy_exact; qx_exact];

% Error computation
error_p = norm(p_mimetic - p_exact) / norm(p_exact);
error_q = norm(q_mimetic - q_exact) / norm(q_exact);
disp('Error Pressure'); disp(error_p);
disp('Error Flux'); disp(error_q);

#% Plot Numerical Pressure
#figure();
#surf(X, Y, reshape(p_mimetic, nx+2, ny+2));
#title(['Numerical pressure, k=', num2str(k)]);
#xlabel('x'); ylabel('y'); colorbar;
#
#% Plot Exact Pressure
#figure();
#surf(X, Y, p(X, Y));
#title('Exact pressure');
#xlabel('x'); ylabel('y'); colorbar;
#
#% Plot Exact Vertical Fluxes
#figure();
#surf(Xdualh, Ydualh, qy(Xdualh, Ydualh));
#title('Exact vertical fluxes');
#xlabel('x'); ylabel('y'); colorbar;
#
#% Plot Numerical Vertical Fluxes
#figure();
#surf(Xdualh, Ydualh, reshape(qy_mimetic, ny+1, nx));
#title('Numerical vertical fluxes');
#xlabel('x'); ylabel('y'); colorbar;
#
#% Plot Exact Horizontal Fluxes
#figure();
#surf(Xdualv, Ydualv, qx(Xdualv, Ydualv));
#title('Exact horizontal fluxes');
#xlabel('x'); ylabel('y'); colorbar;
#
#% Plot Numerical Horizontal Fluxes
#figure();
#surf(Xdualv, Ydualv, reshape(qx_mimetic, ny, nx+1));
#title('Numerical horizontal fluxes');
#xlabel('x'); ylabel('y'); colorbar;


# Testing interpolators

# --> Center to faces interpolators

C2F = interpolCentersToFacesD2D(k, ncells, ncells);
p_faces = C2F * [p_mimetic; p_mimetic];
p_vertical_edges = p_faces(end/2+1:end);
p_horizontal_edges = p_faces(1:end/2);

# Pressure at the vertical edges
#figure();
#surf(Xdualv, Ydualv, reshape(p_vertical_edges, ny, nx+1));
#title('Pressure at vertical edges')

#figure();
#surf(Xdualh, Ydualh, reshape(p_horizontal_edges, ny+1, nx));
#title('Pressure at horizontal edges')

# --> Center to nodes interpolators
x_nodes = west:dx:east;
y_nodes = south:dy:north;
[X_nodes, Y_nodes] = meshgrid(x_nodes, y_nodes);
C2N = interpolCentersToNodes2D(k, ncells, ncells);
p_nodes = C2N * p_mimetic;
figure();
surf(X_nodes, Y_nodes, reshape(p_nodes, nx+1, ny+1));






