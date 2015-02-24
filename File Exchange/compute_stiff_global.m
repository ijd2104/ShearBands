% Assemble the global stiffness matrix using the local stiffness matrices

function stiff_global = compute_stiff_global(stiff_local, num_edges)

num_elements = num_edges - 1;

stiff_global = zeros(num_elements);

for i = 2 : 1 : num_elements
    
    stiff_global(i - 1,i - 1) = stiff_global(i - 1,i - 1) + stiff_local(i,1);
    
    stiff_global(i - 1,i) = stiff_global(i - 1,i) + stiff_local(i,2);
    
    stiff_global(i,i - 1) = stiff_global(i,i - 1) + stiff_local(i,3);
    
    stiff_global(i,i) = stiff_global(i,i) + stiff_local(i,4);
    
end

stiff_global(1,1) = stiff_global(1,1) + stiff_local(1,4);

end