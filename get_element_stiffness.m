function [r,m,k] = get_element_stiffness( xe,flg )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Linear terms
    mv = zeros(2,2);
    mT = zeros(2,2);
    ms = 0;
    mg = 0;
    kvs = zeros(2,1);
    ksv = zeros(1,2);
    kTT = zeros(2,2); 
    
    ngp = 2; %number of Gauss points
    [W,xi] = gaussQuad(ngp);
    for i = 1:ngp
        [N,B] = get_shape_functions(xi(i),h);
        
        %Mass matrix
        mv = mv+W(i)*rho*(N'*N)*J;
        mT = mT+W(i)*rho*cp*(N'*N)*J;
        ms = ms+W(i)*J;
        mg = mg+W(i)*J;
        
        %Linear Stiffness matrix
        kvs = kvs-W(i)*B'*J;
        ksv = ksv+W(i)*E*B*J;
        kTT = kTT-W(i)*lambda*(B'*B)*J;
        
        %Non linear stiffness matrix
        [dgdg,dgds,dgdT] = get_plastic_strain_rate();
    end

%% Non linear terms
    Gss = 0;
    GsT = zeros(1,2);
    Gsg = 0;
    
    GTs = zeros(2,1);
    GTT = zeros(2,2);
    GTg = zeros(2,1);
    
    Ggs = 0;
    GgT = zeros(1,2);
    Ggg = 0;
    
    ngp = 2; %number of Gauss points
    [W,xi] = gaussQuad(ngp);
    for i = 1:ngp
        N = (1/2)*[1-xi(i) 1+xi(i)];
        N_T = N';
        
        g = get_g;
        dgdp = get_dgdp(g);
        dgds = get_dgds(g);
        dgdT = get_dgdT(g);
        
        Gss = Gss-W(i)*E*dgds*J;
        GsT = GsT-W(i)*E*dgdT*N*J;
        Gsg = Gsg-W(i)*E*dgdp*J;
        
        GTs = GTs+W(i)*N';
        GTT = GTT+W(i)*(N'*N)*dgdT;
        GTg = GTg+W(i)*N'*dgdp;
        
        Ggs = Ggs+W(i)*dgds*J;
        GgT = GgT+W(i)*dgdT*N*J;
        Ggg = Ggg+W(i)*dgdp*J;
    end

%% Other
    
end

function [B,N] = get_shape_functions(xi,h)
    %For velocity and temperature
    N = (1/2)*[1-xi 1+xi];
    B = (1/h)*[-1, 1];
end

function [dgdg,dgds,dgdT] = get_plastic_strain_rate(model)
    %Input variables
    %eps0dot            reference strain rate
    %sig0               yield stress
    %eps0               yield strain
    %n                  strain hardening exponent
    %T0                 reference temperature
    %m                  thermal softening exponent
    
    switch model
        case 1
            %Litonski-Batra Model
        case 2
            %Power Law Model
        case 3       
            %Johnson-Cook Model
            x_psi = (xe(1)*(1-psi)+xe(2)*(1+psi))/2;
            T_psi = (Te(1)*(1-psi)+Te(2)*(1+psi))/2;

            g = 
        case 4
            %Bodner-Parton Model
        otherwise
            error('Incorrect constitutive model')
    end
end