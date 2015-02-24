% Local element integration of source function

function RHS_local = compute_RHS_local(quad_weights, dof_pos, act_quad_pos, ev, num_edges)

% Number of elements
num_elements = num_edges - 1;

% Initialize storage of local source integrations
RHS_local = zeros(num_elements, 2);

% Go over each element
for i = 1 : 1 : num_elements
    
    % Go over each shape function
    for j = 1 : 1 : 2
        
        % Go over each Gauss quadrature
        for k = 1 : 1 : 2
            
            % Special case for last degree of freedom because of natural
            % boundary condition i.e., the a(2) term in the RHS of eqn 6
            if i == num_elements && j == 2
                
                RHS_local(i,j) = RHS_local(i,j) ...
                    + (dof_pos(i + 1) - dof_pos(i)) * quad_weights(k) ...
                    * (F_def(act_quad_pos(k,i)) * ev(j,k) - ev(j,k)) ...
                    + quad_weights(k) * 3 * a_def(2);
                
            else
                
                RHS_local(i,j) = RHS_local(i,j) ...
                    + (dof_pos(i + 1) - dof_pos(i)) * quad_weights(k) ...
                    * (F_def(act_quad_pos(k,i)) * ev(j,k) - ev(j,k));
                
            end            
            
        end
        
    end
    
end

end