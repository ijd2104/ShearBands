function [Residual,fv_new,fsig_new,fT_new,fgamp_new,normres] = evaluate_residual(X,nnode,nel,dt,IEN,Mv,Kv,Msig,Ksig,MT,KT,Mgamp,V0,stress0,T0,plastic_strain0,fv0,fsig0,fT0,fgamp0,prop,Velbc,x,imper)
proplist;

V_new = X(1:nnode);
stress_new = X(nnode+1:nnode+nel);
T_new = X(nnode+nel+1:2*nnode+nel);
plastic_strain_new = X(nel+2*nnode+1:2*nel+2*nnode);

fv_new = zeros(nnode,1);
fsig_new = zeros(nel,1);
fT_new = zeros(nnode,1);
fgamp_new = zeros(nel,1);

for e = 1:nel
    [ sctr ]= getsctr( e, IEN );
    [fve,fsige,fTe,fgampe] = get_forcings(e,stress_new(e),plastic_strain_new(e),T_new(sctr(:)),prop,IEN,x,imper);
    fv_new(sctr(:)) = fv_new(sctr(:)) + fve;
    fsig_new(e) = fsig_new(e) + fsige;
    fT_new(sctr(:)) = fT_new(sctr(:)) + fTe;
    fgamp_new(e) = fgamp_new(e) + fgampe;
end


%momentum
res_mo = Mv*(V_new - V0) + Kv*dt*((1-alphaparam)*stress0+alphaparam*stress_new);
res_mo(1) = V_new(1) - Velbc(1);
res_mo(nnode) = V_new(nnode) - Velbc(2);

%elasticity
res_el = Msig*(stress_new - stress0) - Ksig*dt*((1-alphaparam)*V0+alphaparam*V_new)+dt*((1-alphaparam)*fsig0 + alphaparam*fsig_new);
%energy
res_en = MT*(T_new - T0) + KT*dt*((1-alphaparam)*T0 + alphaparam*T_new)-dt*((1-alphaparam)*fT0 + alphaparam*fT_new);
%flowlaw
res_fl = Mgamp*(plastic_strain_new - plastic_strain0) - dt*((1-alphaparam)*fgamp0 + alphaparam*fgamp_new);
%Residual = Ainv*Residual;
%

%}

Residual = vertcat(res_mo,res_el,res_en,res_fl);

res_mo = res_mo./abs(Velbc(2));
res_el = res_el./sig0;
res_en = res_en./t0;
res_fl = res_fl./eps0;

normres = norm(vertcat(res_mo,res_el,res_en,res_fl));
