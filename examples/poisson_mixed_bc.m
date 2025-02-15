% This code implements the Poisson Problem with mixedBC
% Left bc   -> u = 1
% Right bc  -> u = 0
% Bottom bc -> u' = 0
% Top bc    -> u' = 0

% Clear and add path
clear all; clc; close all;
addpath('/Users/jvmini/Git/mole-master/src/mole_MATLAB')

% Set domain's limits
left   = 0;
right  = 1;
bottom = 0;
top    = 1;

m = 8;
n = 10;
dx = (right - left) / m;
dy = (top - bottom) / n;

%  Construct Laplacian operator
k = 2;
L = lap2D(k, m, dx, n, dy);

% Set up mixed_bc boundary
mixed_bc = mixedBC2D(k,  ...
                     m,  dx, ...
                     n,  dy, ...
                     "Robin", [1, 0], ...  % Dirichlet on the left side
                     "Robin", [1, 0], ...  % Dirichlet on the right side
                     "Robin", [0, 1], ...  % Neumann on the bottom side
                     "Robin", [0, 1]);     % Neumann on the top side

% Update Laplacian to account for boundary conditions
L = L + mixed_bc;

% Assemble RHS and solve

b = zeros(m+2, n+2)';
b(:, 1) = 100;  % known value at the left boundary
b = reshape(b', [], 1);

p = L \ b;

% Plot Numerical Pressure

% Construct staggered grid
x_grid = [left left+dx/2 : dx : right-dx/2 right];
y_grid = [bottom bottom+dy/2 : dy : top-dy/2 top];
[X, Y] = meshgrid(x_grid, y_grid);

% Plot
figure();
surf(X, Y, reshape(p, m+2, n+2)');
title(['Numerical pressure, k=', num2str(k)]);
xlabel('x'); ylabel('y'); colorbar;

% Retrieve flux solution
G = grad2D(k, m, dx, n, dy);
q = - G * p;

% Grid for vertical and horizontal fluxes
x_dual_grid_h = left + dx/2 : dx : right - dx/2;
y_dual_grid_h = bottom : dy : top;
[Xdualh, Ydualh] = meshgrid(x_dual_grid_h, y_dual_grid_h);
qy = q(1:length(q)/2);

x_dual_grid_v = left : dx : right;
y_dual_grid_v = bottom + dy/2 : dy : top - dy/2;
[Xdualv, Ydualv] = meshgrid(x_dual_grid_v, y_dual_grid_v);
qx = q(length(q)/2+1:end);








