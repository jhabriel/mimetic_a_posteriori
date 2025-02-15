% Implement a convergence analysis for mimetic difference approximation
% to the Darcy problem with unit permeability, using a trigonometric
% manufactured solution

ncells = [8, 16, 32, 64];
k = [2];

errors_k2_pcenter = zeros(length(ncells), 1);
#errors_k4_p = zeros(length(ncells), 1);
errors_k2_pnodes = zeros(length(ncells), 1);
#errors_k4_q = zeros(length(ncells), 1);

for ii=1:length(k)
    for jj=1:length(ncells)

        [error_pcenter, error_pnodes] = darcy_unit_perm_model(ncells(jj), k(ii));

        if k(ii) == 2
            errors_k2_pcenter(jj) = error_pcenter;
            errors_k2_pnodes(jj) = error_pnodes;
        else
            errors_k4_p(jj) = error_p;
            errors_k4_q(jj) = error_q;
        end
    end
end

