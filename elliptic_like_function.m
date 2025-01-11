%Solve the elliptic-like
function [tar, u, cost_time] = elliptic_like_function(N, m, a, b, alpha, u_exa)

%Initial the boundary geometry elliptic-like
tic
gamma = {@(t, r) a * r .* cos(t);
         @(t, r) b * r .* sin(t) + alpha * r .* sin(t) .^3};
gamma_dt = {@(t, r) -a * r .* sin(t);
            @(t, r) b * r .* cos(t) + 3 * alpha * r .* sin(t) .^ 2 .* cos(t)};
gamma_ddt = {@(t, r) -1 .* a * cos(t);
             @(t, r) -1 .* b * sin(t) + 3 * alpha * (2 * sin(t) .* cos(t) .^ 2 - sin(t) .^ 3)};
gamma_dt_fixedr = {@(t) gamma_dt{1}(t, 1);
                   @(t) gamma_dt{2}(t, 1)};


%Determine the source and target
theta = (0 : N - 1)' * 2 * pi / N;
src = [gamma{1}(theta, 1) gamma{2}(theta, 1)];

rho_cor = 0.1;
rho_cor = ((1 : m) / (m + 1)) .^ (1 / (1 + rho_cor));
tar = [gamma{1}(theta, rho_cor) gamma{2}(theta, rho_cor)];
tar = reshape(tar, [], 2);
[M, ~] = size(tar);


%solve the equation
u = IEB_solver(u_exa, theta, src, N, tar, M, gamma_dt_fixedr, gamma_ddt);
cost_time = toc;
end

