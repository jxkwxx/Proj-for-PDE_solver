%Solve the elliptic

%Initial the boundary geometry
a = 7;
b = 5;
alpha = 1.3;


%Initial the exact solution
u_exa = rand(3, 1);
u_exa = u_exa / sum(u_exa);
u_exa = @(x, y) 10 + u_exa(1) / (a * b) * x .* y + u_exa(2) / exp(b) * exp(-y) .* sin(x) + u_exa(3) / exp(a) * exp(x) .* cos(y);


%Determine the number of sources and targets
N = 100;
m = 400;


%Solve the equation
%[tar, u, cost_time] = elliptic_like_function(N, m, a, b, alpha, u_exa);
[tar, u, cost_time] = elliptic_like_function_cor(N, m, a, b, alpha, u_exa);


%Plot the graph
bd = @(x, y) x .^ 2 / a^2 + y .^ 2 / b^2 - 2 * alpha / b * (1 - x .^ 2 / a ^ 2) .^ 2 - alpha ^ 2 / b ^ 2 * (1 - x .^ 2 / a ^ 2) .^ 3;
axis_range = [-a a; -b - alpha b + alpha];
N_num = [500, 500];
colorbar_range = -15 : 2 : 0;
N_contourf = length(colorbar_range);
[u_graph] = graph_plot(tar, log10(abs(u - u_exa(tar(:, 1), tar(:, 2))) ./ u_exa(tar(:, 1), tar(:, 2)) + 1e-15), bd, axis_range, N_num, N_contourf, colorbar_range);



