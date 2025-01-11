%FDM method
function [U, cost_time, err] = FDM_solver(a, b, nx, ny, u_exa, bd, tol)

tic

x = linspace(-a, a, nx);
y = linspace(-b, b, ny);
[X, Y] = meshgrid(x, y);

U = zeros(ny, nx);

for j = 1 : nx
    for i = 1 : ny
        if bd(X(i, j), Y(i, j)) >= 1
            U(i, j) = u_exa(X(i, j), Y(i, j));
        end
    end
end

err = 1;
num_iter = 0;

while err > tol
    U_pre = U;

    for i = 2 : ny - 1
        for j = 2 : nx - 1
            if bd(X(i, j), Y(i, j)) <1
                U(i, j) = 0.25 * (U_pre(i + 1, j) + U_pre(i - 1, j) + U_pre(i, j + 1) + U_pre(i, j - 1));
            end
        end
    end
 
    err = max(max(abs(U - U_pre)));
    num_iter = num_iter + 1;
end

cost_time = toc;

err = abs(u_exa(X, Y) - U) ./ abs(u_exa(X, Y));
err = norm(err(:)) / nx * 2;
end

