function [X, Y] = staggered_grid_2d(step_sizes, domain)

    % domain
    west = domain(1);
    east = domain(2);
    south = domain(3);
    north = domain(4);

    % step sizes
    dx = step_sizes(1);
    dy = step_sizes(2);

    % Create staggered grid
    x_staggered = [west west+dx/2 : dx : east-dx/2 east];
    y_staggered = [south south+dy/2 : dy : north-dy/2 north];
    [X, Y] = meshgrid(x_staggered, y_staggered);

end
