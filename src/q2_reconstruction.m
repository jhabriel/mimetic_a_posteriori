function A = Q2_reconstruction(X_nodes, Y_nodes, U_nodes, ...
                               X_mid_x, Y_mid_x, U_mid_x, ...
                               X_mid_y, Y_mid_y, U_mid_y, ...
                               X_center, Y_center, U_center)
  % Q2 reconstruction for a Cartesian grid using nodal, mid-edge, and center values
  %
  % INPUT:
  % X_nodes  - (Nx+1, Ny+1) matrix of x-coordinates for corner nodes
  % Y_nodes  - (Nx+1, Ny+1) matrix of y-coordinates for corner nodes
  % U_nodes  - (Nx+1, Ny+1) matrix of function values at corner nodes
  %
  % X_mid_x  - (Nx, Ny+1) matrix of x-coordinates for midpoints along x-direction
  % Y_mid_x  - (Nx, Ny+1) matrix of y-coordinates for midpoints along x-direction
  % U_mid_x  - (Nx, Ny+1) matrix of function values at midpoints along x-direction
  %
  % X_mid_y  - (Nx+1, Ny) matrix of x-coordinates for midpoints along y-direction
  % Y_mid_y  - (Nx+1, Ny) matrix of y-coordinates for midpoints along y-direction
  % U_mid_y  - (Nx+1, Ny) matrix of function values at midpoints along y-direction
  %
  % X_center - (Nx, Ny) matrix of x-coordinates for cell centers
  % Y_center - (Nx, Ny) matrix of y-coordinates for cell centers
  % U_center - (Nx, Ny) matrix of function values at cell centers
  %
  % OUTPUT:
  % A - (Nx * Ny, 9) matrix where each row corresponds to a grid cell
  %     and contains the Q2 polynomial coefficients [a0, a1, ..., a8]
  %
  % Each Q2 polynomial is of the form:
  %   Q_2^K(x, y) = a0 + a1*x + a2*y + a3*x*y + a4*x^2 + a5*y^2
  %                 + a6*x^2*y + a7*x*y^2 + a8*x^2*y^2
  %
  % The coefficients are computed by solving a 9x9 linear system for each grid cell.

  Nx = size(X_nodes,1) - 1; % Number of cells in x direction
  Ny = size(Y_nodes,2) - 1; % Number of cells in y direction

  A = zeros(Nx * Ny, 9); % Store coefficients for each cell

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

          % Extract mid-edge coordinates & values
          x2 = X_mid_x(i,j);   y2 = Y_mid_x(i,j);   u2 = U_mid_x(i,j);
          x4 = X_mid_y(i,j);   y4 = Y_mid_y(i,j);   u4 = U_mid_y(i,j);
          x5 = X_mid_y(i+1,j); y5 = Y_mid_y(i+1,j); u5 = U_mid_y(i+1,j);
          x8 = X_mid_x(i,j+1); y8 = Y_mid_x(i,j+1); u8 = U_mid_x(i,j+1);

          % Extract center coordinates & values
          x6 = X_center(i,j); y6 = Y_center(i,j); u6 = U_center(i,j);

          % Construct the system matrix
          M = [1, x1, y1, x1*y1, x1^2, y1^2, x1^2*y1, x1*y1^2, x1^2*y1^2;
               1, x2, y2, x2*y2, x2^2, y2^2, x2^2*y2, x2*y2^2, x2^2*y2^2;
               1, x3, y3, x3*y3, x3^2, y3^2, x3^2*y3, x3*y3^2, x3^2*y3^2;
               1, x4, y4, x4*y4, x4^2, y4^2, x4^2*y4, x4*y4^2, x4^2*y4^2;
               1, x5, y5, x5*y5, x5^2, y5^2, x5^2*y5, x5*y5^2, x5^2*y5^2;
               1, x6, y6, x6*y6, x6^2, y6^2, x6^2*y6, x6*y6^2, x6^2*y6^2;
               1, x7, y7, x7*y7, x7^2, y7^2, x7^2*y7, x7*y7^2, x7^2*y7^2;
               1, x8, y8, x8*y8, x8^2, y8^2, x8^2*y8, x8*y8^2, x8^2*y8^2;
               1, x9, y9, x9*y9, x9^2, y9^2, x9^2*y9, x9*y9^2, x9^2*y9^2];

          % Right-hand side (nodal values)
          b = [u1; u2; u3; u4; u5; u6; u7; u8; u9];

          % Solve for Q2 coefficients
          a = M \ b;

          % Store in output matrix
          A(cell_idx, :) = a';
          cell_idx = cell_idx + 1;
      end
  end
end

