%solve over the boundary
function [mu, nu_y] = bd_solver(u_exa, theta, src, N, gamma_dt, gamma_ddt)

f = u_exa(src(:, 1), src(:, 2));

nu_y = [gamma_dt{2}(theta) -gamma_dt{1}(theta)];

A = [src(:,1) - src(:,1)' src(:,2) - src(:,2)'];
A = (A(:, 1 : N) .* nu_y(:,1)' + A(:, N + 1 : end) .* nu_y(:,2)') ./ (A(:, 1 : N) .^ 2 + A(:, N + 1 : end) .^ 2);
A(1 : N + 1 : end) = 1/2 * (gamma_ddt{1}(theta) .* nu_y(:,1) + gamma_ddt{2}(theta) .* nu_y(:,2)) ./ (sum(nu_y .^2, 2));
A = 1/N * A - 1/2 * eye(N);

mu = A\f;
end