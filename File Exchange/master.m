%{
  Title: 1D Finite Element Method (FEM) Tutorial
  Author: Cheuk Lau
  Date: 2/8/2014
  Description: This is the master file of the 1D FEM tutorial. It controls
  the flow of the program. Please read the README.pdf for more information
  on each of the required files!
  Files required: (1) master.m
                  (2) a_def.m
                  (3) act_quad.m
                  (4) compute_dof_pos.m
                  (5) compute_RHS_global.m
                  (6) compute_RHS_local.m
                  (7) compute_stiff_global.m
                  (8) compute_stiff_local.m
                  (9) eval_der_shape.m
                  (10) eval_shape.m
                  (11) f_def.m
                  (12) loc_to_glob.m
                  (13) ref_quad.m
  Reference: (1) http://mathworld.wolfram.com/FiniteElementMethod.html
  Input parameters: Run program directly from master.m
  Output: Plot of analytical and computed solutions.  
%}


%% User input
h = 1 / 32; % Element width

num_edges = (h ^ -1) * 2 + 1; % Number of element edges

%% Two-point Gauss quadrature over reference element
[ref_quad_pos, quad_weights] = ref_quad();

%% Degrees of freedom locations i.e., finite element edges
dof_pos = compute_dof_pos(num_edges);

%% Reference finite element shape function values at quadrature nodes
ev = eval_shape(ref_quad_pos);

%% Reference finite elelment shape function derivative values at quadrature nodes
evder = eval_der_shape();

%% Reference element quadrature nodes mapped to each physical element
act_quad_pos = act_quad(dof_pos, ref_quad_pos, num_edges);

%% Local to global degree of freedom data structure
loc_glob = loc_to_glob(num_edges);

%% Local right-hand side integration
RHS_local = compute_RHS_local(quad_weights, dof_pos, act_quad_pos, ev, num_edges);

%% Global right-hand side integration
RHS_global = compute_RHS_global(RHS_local, loc_glob, num_edges);

%% Local stiffness matrices
stiff_local = compute_stiff_local(dof_pos, quad_weights, act_quad_pos, ev, evder, num_edges);

%% Global stiffness matrix
stiff_global = compute_stiff_global(stiff_local, num_edges);
stiff_global_sparse = sparse(stiff_global);

%% Remove the first RHS global element since first DOF is known (BC)
RHS_use = RHS_global(2:end);

%% Use GMRES iterative solver to calculate the solution
solution = gmres(stiff_global_sparse, RHS_use, 1e6);

%% Shift solution up by one to meet left BC
solution = solution + 1;

%% Insert left BC back into the solution for plotting
solution = [1; solution(1:end)];

%% Plot the solution and reference solution
mesh = 0 : h : 2; %% Spatial mesh
computed = plot(mesh, solution, 's');
set(computed, 'MarkerSize', 10);
hold on;
analytical = 1 + 3 .* mesh + mesh .* (2 - mesh) .^ 2;
analytical_plot = plot(mesh, analytical);
set(analytical_plot, 'LineWidth', 1.25, 'color', 'r');
grid on;
xlabel('x', 'FontSize', 18);
ylabel('u(x)', 'FontSize', 18);
legend('FEM', 'Analytical', 'Location', 'Southeast');