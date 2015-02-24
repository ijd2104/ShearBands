% Reference element shape function derivatives at quadrature nodes

function evder = eval_der_shape()

evder = zeros(2, 2);

for i = 1 : 1 : 2
    
    evder(1,i) = -1;
    
    evder(2,i) = 1;
    
end

end