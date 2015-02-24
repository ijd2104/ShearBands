% Generate the positions of each degree of freedom positions in the problem domain

function dof_pos = compute_dof_pos(num_edges)

% Degrees of freedom positions over problem domain
dof_pos = zeros(2, 1);

% Leftmost degree of freedom at left boundary
dof_pos(1) = 0;

% Width of each finite element
h = 2 / (num_edges - 1);

% Remaining degrees of freedom locations
for i = 2 : 1 : num_edges
    
    dof_pos(i) = dof_pos(i - 1) + h;

end