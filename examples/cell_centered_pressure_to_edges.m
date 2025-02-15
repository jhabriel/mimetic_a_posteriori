function [p_hori, p_vert] = cell_centered_pressure_to_edges(k, ncells, p)

    nx = ncells(1);
    ny = ncells(2);

    C2F = interpolCentersToFacesD2D(k, nx, ny);
    p_edges = C2F * [p; p];  % stack pressure solution
    p_vert = p_edges(1:((nx+1)*ny));    % pressure at the vertical edges
    p_hori = p_edges((nx+1)*ny+1:end);  % pressure at the horizontal edges

    % The pressure at the vertical edges correspond to the x-component of a
    % vector quanitity, whereas the pressure at the horizontal edges correspond
    % to the y-component of a vector quantity. That's why the first portion
    % of p_edges correspond to the vertical edges, whereas the second portion of
    % p_edges correspond to the horizontal edges.

end
