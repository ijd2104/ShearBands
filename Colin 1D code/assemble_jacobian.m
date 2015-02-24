function [ JAC, Picard ] = assemble_jacobian( X, nel, nnode, Mv, Kv, Ksig, Msig, MT, KT, Mgamp, prop, dt, IEN, x, imper, Velbc )
proplist;
globvarssb;

stress = X(nnode+1:nnode+nel);
T = X(nnode+nel+1:2*nnode+nel);
plastic_strain = X(nel+2*nnode+1:2*nel+2*nnode);

Gsigsig = zeros(nel,nel);
GsigT = zeros(nel,nnode);
Gsiggamp = zeros(nel,nel);

GTsig = zeros(nnode,nel);
GTT = zeros(nnode,nnode);
GTgamp = zeros(nnode,nel);

Ggampsig = zeros(nel,nel);
GgampT = zeros(nel,nnode);
Ggampgamp = zeros(nel,nel);



for e = 1:nel
    [ sctr ]= getsctr( e, IEN );
    [Gsigsige,GsigTe,Gsiggampe,GTsige,GTTe,GTgampe,Ggampsige,GgampTe,Ggampgampe] = get_nonlinear_stiffnesses(e,stress(e),plastic_strain(e),T(sctr(:)),prop,IEN,x,imper);
    Gsigsig(e,e) = Gsigsig(e,e) + Gsigsige;
    GsigT(e,sctr(:)) = GsigT(e,sctr(:)) + GsigTe;
    Gsiggamp(e,e) = Gsiggamp(e,e) + Gsiggampe;
    
    GTsig(sctr(:),e) = GTsig(sctr(:),e) + GTsige;
    GTT(sctr(:),sctr(:)) = GTT(sctr(:),sctr(:)) + GTTe;
    GTgamp(sctr(:),e) = GTgamp(sctr(:),e) + GTgampe;
    
    Ggampsig(e,e) = Ggampsig(e,e) + Ggampsige;
    GgampT(e,sctr(:)) = GgampT(e,sctr(:)) + GgampTe;
    Ggampgamp(e,e) = Ggampgamp(e,e) + Ggampgampe;    
end


%{
JAC = [Mv./Velbc alphaparam*dt*Kv./Velbc zeros(nnode,nnode) zeros(nnode,nel);...
    -alphaparam*dt*Ksig./sig0 (Msig+alphaparam*dt*Gsigsig)./sig0 -alphaparam*dt*GsigT./sig0 -alphaparam*dt*Gsiggamp./sig0;...
    zeros(nnode,nnode) -alphaparam*dt*GTsig./t0 (MT+alphaparam*dt*(KT+GTT))./t0 alphaparam*dt*GTgamp./t0;...
    zeros(nel,nnode) -alphaparam*dt*Ggampsig./eps0 alphaparam*dt*GgampT./eps0 (Mgamp+alphaparam*dt*Ggampgamp)./eps0];
JAC(1,:) = zeros(1,2*nnode+2*nel);
JAC(nnode,:) = zeros(1,2*nnode+2*nel);

JAC(1,1) = 1/Velbc;
JAC(nnode,nnode) = 1/Velbc;
Picard = [Mv./Velbc alphaparam*dt*Kv./Velbc zeros(nnode,nnode) zeros(nnode,nel);...
       -alphaparam*dt*Ksig./sig0 Msig./sig0 zeros(nel,nnode) zeros(nel,nel);...
       zeros(nnode,nnode) zeros(nnode,nel) (MT+alphaparam*dt*(KT))./t0 zeros(nnode,nel);...
       zeros(nel,nnode) zeros(nel,nel) zeros(nel,nnode) Mgamp./eps0];
Picard(1,:) = zeros(1,2*nnode+2*nel);
Picard(nnode,:) = zeros(1,2*nnode+2*nel);

Picard(1,1) = 1/Velbc;
Picard(nnode,nnode) = 1/Velbc;

%}
%
JAC = [Mv alphaparam*dt*Kv zeros(nnode,nnode) zeros(nnode,nel);...
    -alphaparam*dt*Ksig (Msig+alphaparam*dt*Gsigsig) -alphaparam*dt*GsigT -alphaparam*dt*Gsiggamp;...
    zeros(nnode,nnode) -alphaparam*dt*GTsig MT+alphaparam*dt*(KT+GTT) alphaparam*dt*GTgamp;...
    zeros(nel,nnode) -alphaparam*dt*Ggampsig alphaparam*dt*GgampT Mgamp+alphaparam*dt*Ggampgamp];
JAC(1,:) = zeros(1,2*nnode+2*nel);
JAC(nnode,:) = zeros(1,2*nnode+2*nel);
JAC(1,1) = 1;
JAC(nnode,nnode) = 1;


Picard = [Mv alphaparam*dt*Kv zeros(nnode,nnode) zeros(nnode,nel);...
       -alphaparam*dt*Ksig Msig zeros(nel,nnode) zeros(nel,nel);...
       zeros(nnode,nnode) zeros(nnode,nel) (MT+alphaparam*dt*(KT)) zeros(nnode,nel);...
       zeros(nel,nnode) zeros(nel,nel) zeros(nel,nnode) Mgamp];
Picard(1,:) = zeros(1,2*nnode+2*nel);
Picard(nnode,:) = zeros(1,2*nnode+2*nel);

Picard(1,1) = 1;
Picard(nnode,nnode) = 1;
%}
%{
Pre = [JAC(1:nnode,1:nnode) zeros(nnode,nel) zeros(nnode,nnode) zeros(nnode,nel);...
       zeros(nel,nnode) JAC(nnode+1:nel+nnode,nnode+1:nel+nnode) zeros(nel,nnode) zeros(nel,nel);...
       zeros(nnode,nnode) zeros(nnode,nel) JAC(nnode+nel+1:nel+2*nnode,nel+nnode+1:nel+2*nnode) zeros(nnode,nel);...
       zeros(nel,nnode) zeros(nel,nel) zeros(nel,nnode) JAC(2*nnode+nel+1:2*nel+2*nnode,nel+2*nnode+1:2*nel+2*nnode)];

%}
%
 

%}
JAC = sparse(JAC);
Picard = sparse(Picard);