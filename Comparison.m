%Comparison

%Set the number
N = 20 : 5 : 100;
m = 400;

%Case for the elliptic-like
a = 3;
b = 2;
alpha = 0.7;

u_exa = rand(3, 1);
u_exa = u_exa / sum(u_exa);
u_exa = @(x, y) 10 + u_exa(1) / (a * b) * x .* y + u_exa(2) / exp(b) * exp(-y) .* sin(x) + u_exa(3) / exp(a) * exp(x) .* cos(y);

bd = @(x, y) x .^ 2 / a^2 + y .^ 2 / b^2 - 2 * alpha / b * (1 - x .^ 2 / a ^ 2) .^ 2 - alpha ^ 2 / b ^ 2 * (1 - x .^ 2 / a ^ 2) .^ 3;

f = @(location, state) u_exa(location.x, location.y);


cost_times_1 = zeros(length(N), 2);
errors_1 = zeros(length(N), 2);

for i = 1 : length(N)
    [tar, u, cost_times_1(i, 1)] = elliptic_like_function(N(i), m, a, b, alpha, u_exa);
    errors_1(i, 1) = norm(abs(u_exa(tar(:, 1), tar(:, 2)) - u) ./ u_exa(tar(:, 1), tar(:, 2))) / (length(tar)) * 2;
    %[res, cost_times_1(i, 2)] = FEM_solver(@elliptic, f, N(i));
    %u = interpolateSolution(res, tar(:, 1), tar(:, 2));
    %errors_1(i, 2) = norm(abs(u_exa(tar(:, 1), tar(:, 2)) - u) ./ u_exa(tar(:, 1), tar(:, 2))) / (length(tar)) * 2;
    [u, cost_times_1(i, 2), errors_1(i, 2)] = FDM_solver(a, b + alpha, N(i), N(i), u_exa, bd, 1e-5);
end



%Case for the k-leaf
k = 3;
a = 1 + 1 / (k + 1);
b = 1 + 1 / (k + 1);
u_exa = @(x, y) 10 + 1 / (a * b) * x .* y + 1 / exp(b) * exp(-y) .* sin(x) + 1 / exp(a) * exp(x) .* cos(y);

f = @(location, state) u_exa(location.x, location.y);

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

cost_times_2 = zeros(length(N), 2);
errors_2 = zeros(length(N), 2);

for i = 1 : length(N)
    [tar, u, cost_times_2(i, 1)] = k_leaf_function(N(i), m, k, u_exa);
    errors_2(i, 1) = norm(abs(u_exa(tar(:, 1), tar(:, 2)) - u) ./ u_exa(tar(:, 1), tar(:, 2))) / (length(tar)) * 2;
    %[res, cost_times_2(i, 2)] = FEM_solver(@leaf, f, N(i));
    %u = interpolateSolution(res, tar(:, 1), tar(:, 2));
    %errors_2(i, 2) = norm(abs(u_exa(tar(:, 1), tar(:, 2)) - u) ./ u_exa(tar(:, 1), tar(:, 2))) / (length(tar)) * 2;
    [u, cost_times_2(i, 2), errors_2(i, 2)] = FDM_solver(a, b, N(i), N(i), u_exa, bd, 1e-5);
end

function [x, y] = elliptic(bt, t)
switch nargin
    case 0
        x = 5;
    return
    case 1
        A = [0, 2 * pi / 5, 4 * pi / 5, 6 * pi / 5, 8 * pi / 5; 
             2 * pi / 5, 4 * pi / 5, 6 * pi / 5, 8 * pi / 5, 2 * pi; 
             1, 1, 1, 1, 1; 
             0, 0, 0, 0, 0]; 
        x = A(:,bt); 
    return
    case 2
        x = 3 .* cos(t);
        y = 2 .* sin(t) + 0.7 .* sin(t) .^ 3;
end
end

function [x, y] = leaf(bt, t)
switch nargin
    case 0
        x = 5;
    return
    case 1
        A = [0, 2 * pi / 5, 4 * pi / 5, 6 * pi / 5, 8 * pi / 5; 
             2 * pi / 5, 4 * pi / 5, 6 * pi / 5, 8 * pi / 5, 2 * pi; 
             1, 1, 1, 1, 1; 
             0, 0, 0, 0, 0]; 
        x = A(:,bt); 
    return
    case 2
        x = cos(t) .* (1 + cos(3 * t) / 4);
        y = sin(t) .* (1 + cos(3 * t) / 4);
end
end


%Plot
figure(1)
plot(N, log10(errors_1(:, 1)), 'LineWidth', 2)
hold on 
plot(N, log10(errors_1(:, 2)), 'LineWidth', 2)
xlabel('N')
ylabel('log_{10}(error)')
legend('IEB','FDM')


figure(2)
plot(N, log10(cost_times_1(:, 1)), 'LineWidth', 2)
hold on 
plot(N, log10(cost_times_1(:, 2)), 'LineWidth', 2)
xlabel('N')
ylabel('log_{10}(time)')
legend('IEB','FDM')


figure(11)
plot(N, log10(errors_2(:, 1)), 'LineWidth', 2)
hold on 
plot(N, log10(errors_2(:, 2)), 'LineWidth', 2)
xlabel('N')
ylabel('log_{10}(error)')
legend('IEB','FDM')


figure(21)
plot(N, log10(cost_times_2(:, 1)), 'LineWidth', 2)
hold on 
plot(N, log10(cost_times_2(:, 2)), 'LineWidth', 2)
xlabel('N')
ylabel('log_{10}(time)')
legend('IEB','FDM')
