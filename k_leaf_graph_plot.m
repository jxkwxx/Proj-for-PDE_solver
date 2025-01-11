%Solve the k-leaf

%set the parameters
k = 9;
N = 100;
m = 400;

%determine the exact solution
a = 1 + 1 / (k + 1);
b = 1 + 1 / (k + 1);
u_exa = rand(3, 1);
u_exa = u_exa / sum(u_exa);
u_exa = @(x, y) 10 + u_exa(1) / (a * b) * x .* y + u_exa(2) / exp(b) * exp(-y) .* sin(x) + u_exa(3) / exp(a) * exp(x) .* cos(y);


%solve the equation
%[tar, u, cost_time] = k_leaf_function(100, 400, k, u_exa);
[tar, u, cost_time] = k_leaf_function_cor(100, 400, k, u_exa);


%Plot the graph
chebshev = zeros(k+1,2);
chebshev(1, 1) = 1;
chebshev(2, 2) = 1;
for i = 3 : k + 1
    if mod(i, 2) == 1
        chebshev(1 : i, 1) = -chebshev(1 : i, 1) + 2 * [0; chebshev(1 : i - 1, 2)];
    else
        chebshev(1 : i, 2) = -chebshev(1 : i, 2) + 2 * [0; chebshev(1 : i - 1, 1)];
    end
end
if mod(i, 2) == 1
    chebshev = chebshev(:, 1);
else
    chebshev = chebshev(:, 2);
end
bd = @(x, y) sqrt(x .^ 2 + y .^ 2) ;
for i = 0 : k
    bd = @(x, y) bd(x, y) - chebshev(i + 1) * (x ./ sqrt(x .^ 2 + y .^ 2)) .^ i / (k + 1);
end

axis_range = [-a a; -b b];
N_num = [500, 500];
colorbar_range = -15 : 2 : 0;
N_contourf = length(colorbar_range);
[u_plot] = graph_plot(tar, log10(abs(u - u_exa(tar(:, 1), tar(:, 2))) ./ u_exa(tar(:, 1), tar(:, 2)) + 1e-15), bd, axis_range, N_num, N_contourf, colorbar_range);