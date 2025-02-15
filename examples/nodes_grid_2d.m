function [X, Y] = nodes_grid_2d(step_sizes, domain)

    % domain
    west = domain(1);
    east = domain(2);
    south = domain(3);
    north = domain(4);

    % step sizes
    dx = step_sizes(1);
    dy = step_sizes(2);

    % Create staggered grid
    x_nodes_grid = [west : dx : east];
    y_nodes_grid = [south : dy : north];
    [X, Y] = meshgrid(x_nodes_grid, y_nodes_grid);

end
