%Solve the k-leaf
function [tar, u, cost_time] = k_leaf_function(N, m, k, u_exa)

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
theta = (0 : N - 1)' * 2 * pi / N;
h = 2 * pi / N;
src = [gamma{1}(theta) gamma{2}(theta)];

rho_cor = 0.1;
rho_cor = ((1 : m) / (m + 1)) .^ (1 / (1 + rho_cor));
tar = gamma_coe(theta) * rho_cor;
tar = [tar .* cos(theta) tar .* sin(theta)];
tar = reshape(tar, [], 2);
[M, ~] = size(tar);


%solve the equation
u = IEB_solver(u_exa, theta, src, N, tar, M, gamma_dt, gamma_ddt);
cost_time = toc;
end