function gsig = get_gsig_psi(sige,gampe,Te,xe,imper,prop,psi)

x_psi = (xe(1)*(1-psi)+xe(2)*(1+psi))/2;
T_psi = (Te(1)*(1-psi)+Te(2)*(1+psi))/2;
%set imperfection

prop = set_imperfection(prop,x_psi,imper);
       


proplist;

gsig = eps0dot*m*(sige^(1-1/m)/(sig0*(1+gampe/eps0)^nn*(1-del*(exp((T_psi-t0)/k)-1))))^m;
