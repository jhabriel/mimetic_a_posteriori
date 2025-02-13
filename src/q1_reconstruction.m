function A = Q1_reconstruction(X_nodes, Y_nodes, U_nodes)
  % Q1 reconstruction for a Cartesian grid using nodal values
  %
  % INPUT:
  % X_nodes  - (Nx+1, Ny+1) matrix of x-coordinates for corner nodes
  % Y_nodes  - (Nx+1, Ny+1) matrix of y-coordinates for corner nodes
  % U_nodes  - (Nx+1, Ny+1) matrix of function values at corner nodes
  %
  % OUTPUT:
  % A - (Nx * Ny, 4) matrix where each row corresponds to a grid cell
  %     and contains the Q1 polynomial coefficients [a0, a1, a2, a3]
  %
  % Each Q1 polynomial is of the form:
  %   Q_1^K(x, y) = a0 + a1*x + a2*y + a3*x*y
  %
  % The coefficients are computed by solving a 4x4 linear system for each
  % grid cell.
  %

  Nx = size(X_nodes,1) - 1; % Number of cells in x direction
  Ny = size(Y_nodes,2) - 1; % Number of cells in y direction

  A = zeros(Nx * Ny, 4); % Store coefficients for each cell

  cell_idx = 1; % Counter for grid cell index

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

          % Construct the system matrix
          M = [1, x1, y1, x1*y1;
               1, x3, y3, x3*y3;
               1, x7, y7, x7*y7;
               1, x9, y9, x9*y9];

          % Right-hand side (nodal values)
          b = [u1; u3; u7; u9];

          % Solve for Q1 coefficients
          a = M \ b;

          % Store in output matrix
          A(cell_idx, :) = a';
          cell_idx = cell_idx + 1;
      end
  end
end

