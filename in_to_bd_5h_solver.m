%solve over the boundary 5h-principle
function [u_near] = in_to_bd_5h_solver(beta, p, N, mu, z0, tar_near, gamma, gamma_dt)

N_taylor = floor(beta * N) + 1;
theta_taylor = ((0 : N_taylor - 1) / N_taylor * 2 * pi)';
zn = gamma{1}(theta_taylor) + sqrt(-1) * gamma{2}(theta_taylor);

v = gamma_dt{1}(theta_taylor) + sqrt(-1) * gamma_dt{2}(theta_taylor);
mu_spline = interpft(mu, N_taylor);
v = v .* mu_spline;

coe_mat = zeros(p, length(z0));
for i = 1 : p
    coe_mat0 = (zn - z0.') .^ (-i) .* v;
    coe_mat(i, :) = sum(coe_mat0) * sqrt(-1) / N_taylor;
end
coe_mat = coe_mat.';

tar_comp = tar_near(:, 1) + sqrt(-1) * tar_near(:, 2);
[~, index] = min(abs(tar_comp - z0.'), [], 2);
tar_near_t_z0 = tar_comp - z0(index);

u_near = coe_mat(index, end);
for i = p - 1 : -1 : 1
    u_near = u_near .* tar_near_t_z0 + coe_mat(index, i);
end
u_near = real(u_near);