% Assemble the global right-hand side approximation using local
% contributions

function RHS_global = compute_RHS_global(RHS_local, loc_glob, num_edges)

num_elements = num_edges - 1;

RHS_global = zeros(num_edges, 1);

index = zeros(num_edges, 1);

index(1) = -1 * (num_edges);

index(2) = 1;

for i = 3 : 1 : num_edges
    
    index(i) = index(i - 1) + 1;
    
end

for i = 1 : 1 : num_elements
    
    for j = 1 : 1 : 2
        
        for k = 1 : 1 : num_edges
            
            if loc_glob(i,j) == index(k)
                
                RHS_global(k) = RHS_global(k) + RHS_local(i,j);
                
            end
            
        end
        
    end
    
end

end