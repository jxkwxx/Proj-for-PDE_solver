%Integral equation based solver
function [u] = IEB_solver(u_exa, theta, src, N, tar, M, gamma_dt, gamma_ddt)

%solve over the boundary
[mu, nu_y] = bd_solver(u_exa, theta, src, N, gamma_dt, gamma_ddt);

%solve in the interior
u = in_solver(src, N, tar, M, nu_y, mu);