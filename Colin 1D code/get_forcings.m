% generate the element stiffness, mass and axial force matrix
function [fve,fsige,fTe,fgampe] = get_forcings(e,sige,gampe,Te,prop,IEN,x,imper)
globvarssb;
proplist;
IENe      = IEN(e,:);               % nodes of current element
xe        = x(IENe);                % x coordinates of current element
he        = abs(xe(2) - xe(1));   
J         = he/2;                 % Jacobian of the transformation (1D) 

% initialize element matrices
                  
fve = zeros(2,1);                  
                
fsige = 0;
                 
fTe = zeros(2,1);

fgampe = 0;

ngauss_f = 1;

[wgp_f,pgp_f] = gauss(ngauss_f);

for i = 1:ngauss_f
    
    psi = pgp_f(i);
    N = Nmatrix1D(psi);
    g_psi = get_g_psi(sige,gampe,Te,xe,imper,prop,psi);
    fsige = fsige + E*wgp_f(i)*g_psi*J;
    fTe = fTe + N'*sige*kai/(rho*cp)*wgp_f(i)*g_psi*J;
    fgampe = fgampe + wgp_f(i)*g_psi*J;
end
