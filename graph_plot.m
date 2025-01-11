%Plot the solution
function [u_graph] = graph_plot(tar, u, bd, axis_range, N_num, N_contourf, colorbar_range)
%tar and u: target's coordinate and value;
%bd: distinguish the interior points;
%axis_range and N_num: range and mesh for plotting;

%Plot the mesh
[X, Y] = meshgrid(linspace(axis_range(1,1), axis_range(1,2), N_num(1)), ... 
                  linspace(axis_range(2,1), axis_range(2,2), N_num(2)));
Z = griddata(tar(:,1), tar(:,2), u, X, Y, 'cubic');

%Distinguish the interior points
int_points = bd(X,Y) <= 1;
Z(~int_points) = nan;

%Plot the graph
u_graph = figure;
contourf(X, Y, Z, N_contourf);
hold on;
contour(X, Y, Z, N_contourf, 'k');
clim([colorbar_range(1) colorbar_range(end)]);
cb = colorbar;
cb.Ticks = colorbar_range;
axis equal;
axis off;
hold off;
end