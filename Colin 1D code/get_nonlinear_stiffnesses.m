function [Gsigsige,GsigTe,Gsiggampe,GTsige,GTTe,GTgampe,Ggampsige,GgampTe,Ggampgampe] = get_nonlinear_stiffnesses(e,stress,plastic_strain,T,prop,IEN,x,imper)
globvarssb;
proplist;
IENe      = IEN(e,:);               % nodes of current element
xe        = x(IENe);                % x coordinates of current element
he        = abs(xe(2) - xe(1));   
J         = he/2;                 % Jacobian of the transformation (1D)

Gsigsige = 0;
GsigTe = zeros(1,2);
Gsiggampe = 0;

GTsige = zeros(2,1);
GTTe = zeros(2,2);
GTgampe = zeros(2,1);

Ggampsige = 0;
GgampTe = zeros(1,2);
Ggampgampe = 0;

ngp = 2;
[wgp,pgp] = gauss(ngp);

for i = 1:ngp
    psi = pgp(i);
    g = get_g_psi(stress,plastic_strain,T,xe,imper,prop,psi);
    gsig = get_gsig_psi(stress,plastic_strain,T,xe,imper,prop,psi);
    gT = get_gT_psi(stress,plastic_strain,T,xe,imper,prop,psi);
    N = Nmatrix1D(psi);
    
    Gsigsige = Gsigsige + E*wgp(i)*gsig*J;
    GsigTe = GsigTe + wgp(i)*E*gT*N*J;
    Gsiggampe = Gsiggampe + wgp(i)*E*m*nn/(eps0+plastic_strain)*g*J;
    
    GTsige = GTsige + wgp(i)*N'*(1+m)*kai/(rho*cp)*g*J;
    GTTe = GTTe + wgp(i)*(N'*N)*kai/(rho*cp)*stress*gT*J;
    GTgampe = GTgampe + wgp(i)*N'*kai/(rho*cp)*stress*m*nn/(eps0+plastic_strain)*g*J;
    
    Ggampsige = Ggampsige + wgp(i)*gsig*J;
    GgampTe = GgampTe + wgp(i)*gT*N*J;
    Ggampgampe = Ggampgampe + wgp(i)*m*nn/(eps0+plastic_strain)*g*J;
    
end