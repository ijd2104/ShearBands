function [g_psi] = get_g_psi(sige,gampe,Te,xe,imper,prop,psi)

x_psi = (xe(1)*(1-psi)+xe(2)*(1+psi))/2
T_psi = (Te(1)*(1-psi)+Te(2)*(1+psi))/2;

%set imperfection
prop = set_imperfection(prop,x_psi,imper);



proplist;
g_psi = eps0dot*(sige/(sig0*(1+gampe/eps0)^nn*(1-del*(exp((T_psi-t0)/k)-1))))^m;

%{
if g_psi < eps
    g_psi = 0;
end
%}

end