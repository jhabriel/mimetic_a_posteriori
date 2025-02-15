function [dx, dy] = step_size_2d(ncells, domain)

    % cells
    nx = ncells(1);
    ny = ncells(2);

    % domain
    west = domain(1);
    east = domain(2);
    south = domain(3);
    north = domain(4);

    % Determine step sizes
    dx = (east - west) / nx;
    dy = (north - south) / ny;

end


