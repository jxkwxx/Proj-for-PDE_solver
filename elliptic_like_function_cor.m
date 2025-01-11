%Solve the elliptic-like--5h principle
function [tar, u, cost_time] = elliptic_like_function_cor(N, m, a, b, alpha, u_exa)

%Initial the boundary geometry elliptic-like
tic
gamma = {@(t, r) a * r .* cos(t);
         @(t, r) b * r .* sin(t) + alpha * r .* sin(t) .^3};
gamma_fixedr = {@(t) gamma{1}(t, 1);
                @(t) gamma{2}(t, 1)};
gamma_dt = {@(t, r) -a * r .* sin(t);
            @(t, r) b * r .* cos(t) + 3 * alpha * r .* sin(t) .^ 2 .* cos(t)};
gamma_ddt = {@(t, r) -1 .* a * cos(t);
             @(t, r) -1 .* b * sin(t) + 3 * alpha * (2 * sin(t) .* cos(t) .^ 2 - sin(t) .^ 3)};
gamma_dt_fixedr = {@(t) gamma_dt{1}(t, 1);
                   @(t) gamma_dt{2}(t, 1)};


%Determine the source and target
theta = (0 : N - 1)' * 2 * pi / N;
h = 2 * pi / N;
src = [gamma{1}(theta, 1) gamma{2}(theta, 1)];

rho_cor = 0.1;
rho_cor = ((1 : m) / (m + 1)) .^ (1 / (1 + rho_cor));
tar = [gamma{1}(theta, rho_cor) gamma{2}(theta, rho_cor)];
tar = reshape(tar, [], 2);


%Solve the equation
%solve over the boundary
[mu, nu_y] = bd_solver(u_exa, theta, src, N, gamma_dt_fixedr, gamma_ddt);

%solve in the interior
%partition the points near or away the boundary
bd = @(x, y) x .^ 2 / a^2 + y .^ 2 / b^2 - 2 * alpha / b * (1 - x .^ 2 / a ^ 2) .^ 2 - alpha ^ 2 / b ^ 2 * (1 - x .^ 2 / a ^ 2) .^ 3;

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
z0 = gamma{1}(theta_z0, 1 - 5 * h) + sqrt(-1) * gamma{2}(theta_z0, 1 - 5 * h);

u_near = in_to_bd_5h_solver(beta, p, N, mu, z0, tar_near, gamma_fixedr, gamma_dt_fixedr);


u = [u_away; u_near];
tar = [tar_away; tar_near];

cost_time = toc;
end