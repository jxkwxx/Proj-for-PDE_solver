%FEM method
function [u, cost_time] = FEM_solver(Geo, f, N)

tic

model = createpde();
geometryFromEdges(model, Geo);
applyBoundaryCondition(model, 'dirichlet', 'Edge', 1:model.Geometry.NumEdges, 'u', f);
specifyCoefficients(model, 'm', 0, 'd', 0, 'c', 1, 'a', 0, 'f', 0);
mesh_default = generateMesh(model, 'Hmax', 1 / N);
u = solvepde(model);

cost_time = toc;

end