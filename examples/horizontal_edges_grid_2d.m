function [X, Y] = horizontal_edges_grid_2d(step_sizes, domain)

    % domain
    west = domain(1);
    east = domain(2);
    south = domain(3);
    north = domain(4);

    % step sizes
    dx = step_sizes(1);
    dy = step_sizes(2);

    % Create staggered grid
    x_edges_grid = [west+dx/2 : dx : east-dx/2];
    y_edges_grid = [south : dy : north];
    [X, Y] = meshgrid(x_edges_grid, y_edges_grid);

end
