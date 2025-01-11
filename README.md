This is my projection for the PDE fast solver.

We implement the integral equation based solver through:
    *bd_solver, which solves over the boundary;
    *in_solver, which solves in the interior;
    *IEB_solver, which includes the above two;
    *in_to_bd_5h_solver, which uses the 5h-principle to solve for the points near the boundary;

We implement the FEM and FDM through:
    *FEM_solver;
    *FDM_solver;

We mainly implement our experiments for elliptic-like and k-leaf boundary through:
    *graph_plot, polt the solution;
    *elliptic_like_function, solve use IEB;
    *elliptic_like_function_cor, solve with 5h-principle;
    *elliptic_graph_plot, plot;
    *k_leaf_function, solve use IEB;
    *k_leaf_function, solve with 5h-principle;
    *k_leaf_graph_plot, plot;

