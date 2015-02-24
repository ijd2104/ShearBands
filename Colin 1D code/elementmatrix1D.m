% generate the element stiffness, mass and axial force matrix
function [mve,kve,msige,ksige,mTe,kTe,mgampe,mvle] = elementmatrix1D(e,prop,IEN,x)
globvarssb;
proplist;
IENe      = IEN(e,:);               % nodes of current element
xe        = x(IENe);                % x coordinates of current element
he        = abs(xe(2) - xe(1));   
J         = he/2;                   % Jacobian of the transformation (1D) 

% initialize element matrices
me_CG1 = zeros(2,2);
mvle = zeros(2,2);
kve = zeros(2,1);                  

me_DG0 = 0;
ksige = zeros(1,2);                  

kTe = zeros(2,2);                  

ngauss_m = 2;

[wgp_m,pgp_m] = gauss(ngauss_m);

for i = 1:2  
    psi  =  pgp_m(i);
    N_CG1 = Nmatrix1D(psi);
    B_CG1 = Bmatrix1D(he);
    
    me_CG1 = me_CG1 + wgp_m(i)*(N_CG1'*N_CG1)*J;
    kve = kve + 1/rho*wgp_m(i)*(B_CG1'*1)*J;  
    kTe = kTe + tcond/(rho*cp)*wgp_m(i)*(B_CG1'*B_CG1)*J;
    me_DG0 = me_DG0 + wgp_m(i)*J;
    ksige = ksige + E*wgp_m(i)*(1*B_CG1)*J;
    if e == 11
        N_CG1
    end
end

mve = me_CG1;
mTe = me_CG1;
msige = me_DG0;
mgampe = me_DG0;

mvle(1,1) = sum(mve(1,:));
mvle(2,2) = sum(mve(2,:));
    