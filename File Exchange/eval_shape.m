% Reference element shape function values at quadrature nodes

function ev = eval_shape(ref_quad_pos)
    
ev = zeros(2, 2);

for i = 1 : 1 : 2
                
    ev(1,i) = 1 - ref_quad_pos(i);
            
    ev(2,i) = ref_quad_pos(i);
    
end
    
end