function [p_nodes] = cell_centered_pressure_to_nodes(k, ncells, p_cc)

    C2N = interpolCentersToNodes2D(k, ncells(1), ncells(2));
    p_nodes = C2N * p_cc;

end
