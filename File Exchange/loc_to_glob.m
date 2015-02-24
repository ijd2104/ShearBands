% Assemble the local-to-global degrees of freedom mapping matrix

function loc_glob = loc_to_glob(num_edges)

num_elements = num_edges - 1;

loc_glob = zeros(num_elements, 2);

loc_glob(1,1) = -1 * num_edges;

loc_glob(1,2) = 1;

for i = 2 : 1 : num_elements
    
    loc_glob(i,1) = loc_glob(i - 1, 2);
    
    loc_glob(i,2) = loc_glob(i,1) + 1;
    
end

end