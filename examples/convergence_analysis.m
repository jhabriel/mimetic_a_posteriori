% Implement a convergence analysis for mimetic difference approximation
% to the Darcy problem with unit permeability, using a trigonometric
% manufactured solution

ncells = [10, 20, 40, 80];
k = [2, 4];

errors_k2_p = zeros(length(ncells), 1);
errors_k4_p = zeros(length(ncells), 1);
errors_k2_q = zeros(length(ncells), 1);
errors_k4_q = zeros(length(ncells), 1);

for ii=1:length(k)
    for jj=1:length(ncells)

        [error_p, error_q] = darcy_unit_perm_model(ncells(jj), k(ii));

        if k(ii) == 2
            errors_k2_p(jj) = error_p;
            errors_k2_p(jj) = error_q;
        else
            errors_k4_p(jj) = error_p;
            errors_k4_q(jj) = error_q;
        end
    end
end

