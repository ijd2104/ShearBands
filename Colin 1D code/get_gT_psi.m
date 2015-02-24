function gT = get_gT_psi(stress,plastic_strain,T,xe,imper,prop,psi)

x_psi = (xe(1)*(1-psi)+xe(2)*(1+psi))/2;
       
prop = set_imperfection(prop,x_psi,imper);

proplist;

T_psi = (T(1)*(1-psi)+T(2)*(1+psi))/2;

c1 = m*del/k;
c2 = exp(T_psi/k)/(del*exp(T_psi/k)-(1+del)*exp(t0/k));
gT = c1*c2*get_g_psi(stress,plastic_strain,T,xe,imper,prop,psi);
%{
if gT < eps
    gT = 0;
end 
%}