%Solve the k-leaf--5h principle
function [tar, u, cost_time] = k_leaf_function_cor(N, m, k, u_exa)

%Initial the boundary geometry k_leaf
tic
gamma_coe = @(t) 1 + cos(k * t) / (k + 1);
gamma = {@(t) gamma_coe(t) .* cos(t);
         @(t) gamma_coe(t) .* sin(t)};
gamma_dt = {@(t) (-k / (k + 1) * sin(k * t)) .* cos(t) - gamma_coe(t) .* sin(t);
            @(t) (-k / (k + 1) * sin(k * t)) .* sin(t) + gamma_coe(t) .* cos(t)};
gamma_ddt = {@(t) (-k ^ 2 / (k + 1) * cos(k * t)) .* cos(t) - 2 * (-k / (k + 1) * sin(k * t)) .* sin(t) - gamma_coe(t) .* cos(t);
             @(t) (-k ^ 2 / (k + 1) * cos(k * t)) .* sin(t) + 2 * (-k / (k + 1) * sin(k * t)) .* cos(t) - gamma_coe(t) .* sin(t)};



%Determine the source and target
theta = (0 : N-1)' * 2 * pi / N;
h = 2 * pi / N;
src = [gamma{1}(theta) gamma{2}(theta)];

rho_cor = 0.1;
rho_cor = ((1 : m) / (m + 1)) .^ (1 / (1 + rho_cor));
tar = gamma_coe(theta) * rho_cor;
tar = [tar .* cos(theta) tar .* sin(theta)];
tar = reshape(tar, [], 2);



%solve the equation
%solve over the boundary
[mu, nu_y] = bd_solver(u_exa, theta, src, N, gamma_dt, gamma_ddt);

%solve in the interior
%partition the points near or away the boundary
chebyshev = zeros(k + 1, 2);
chebyshev(1, 1) = 1;
chebyshev(2, 2) = 1;
for i = 3 : k + 1
    if mod(i, 2) == 1
        chebyshev(1 : i, 1) = -chebyshev(1 : i, 1) + 2 * [0; chebyshev(1 : i - 1, 2)];
    else 
        chebyshev(1 : i, 2) = -chebyshev(1 : i, 2) + 2 * [0; chebyshev(1 : i - 1, 1)];
    end
end
if mod(i, 2) == 1
    chebyshev = chebyshev(:, 1);
else
    chebyshev = chebyshev(:, 2);
end

bd = @(x, y) sqrt(x .^ 2 + y .^ 2);
for i = 0 : k
    bd = @(x, y) bd(x, y) - chebyshev(i + 1) * (x ./ sqrt(x .^ 2 + y .^ 2)) .^ i / (k + 1);
end

rho_red = (1 - 5 * h) ^ 2;
tar_par = (bd(tar(:, 1), tar(:, 2)) <= rho_red);
tar_near = [tar(~tar_par, 1) tar(~tar_par, 2)];
tar_near = reshape(tar_near, [], 2);
tar_away = [tar(tar_par, 1) tar(tar_par, 2)];
tar_away = reshape(tar_away, [], 2);
[M_tar_away, ~] = size(tar_away);


%solve the points not close to the boundary
u_away = in_solver(src, N, tar_away, M_tar_away, nu_y, mu);


%solve the points close to the boundary
beta = 4;
p = 15;

N_dom = floor(N / 4) + 1;
theta_z0 = (0 : N_dom - 1)' / N_dom * 2 * pi;
z0 = (1 - 5 * h) * (gamma{1}(theta_z0) + sqrt(-1) * gamma{2}(theta_z0));

u_near = in_to_bd_5h_solver(beta, p, N, mu, z0, tar_near, gamma, gamma_dt);


u = [u_away; u_near];
tar = [tar_away; tar_near];

cost_time = toc;
end