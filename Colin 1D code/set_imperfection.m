function prop = set_imperfection(prop,x_psi,imper)
fluc = 0;
for i = 1:prop(11)
    
    fluc =  fluc + (1-(1-prop(10))* (sech( 1/prop(12)*(x_psi - imper(i)+.0005)) )^2) ; % sig0 : yield stress
end
fluc = fluc - prop(11)+1;
prop(15) = fluc*prop(15);
prop(16) = fluc*prop(16);