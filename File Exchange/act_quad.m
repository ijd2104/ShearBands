% Map reference element quadrature nodes to each physical element

function act_quad_pos = act_quad(dof_pos, ref_quad_pos, num_edges)

% Number of elements
num_elements = num_edges - 1;

% Initialize storage of quadrature positions on physical elements
act_quad_pos = zeros(2, num_elements);

% Go through each degree of freedom
for i = 1 : 1 : 2
    
    % Go through each element
    for j = 1 : 1 : num_elements
        
        act_quad_pos(i,j) = dof_pos(j) + (dof_pos(j + 1) - dof_pos(j)) * ref_quad_pos(i);
        
    end
    
end

end
