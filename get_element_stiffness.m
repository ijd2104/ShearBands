function [r,m,k] = get_element_stiffness( xe,flg )
%GET_ELEMENT_STIFFNESS Computes the terms for FEM matrices of an element
%   Detailed explanation goes here

    %Initialize matrices
    %Mass terms
    m.v = zeros(2,2);
    m.T = zeros(2,2);
    m.s = 0;
    m.g = 0;
    
    %Linear terms
    k.vs = zeros(2,1);
    k.sv = zeros(1,2);
    k.TT = zeros(2,2);
    
    %Non linear terms
    k.ss = 0;
    k.sT = zeros(1,2);
    k.sg = 0;
    
    k.Ts = zeros(2,1);
    k.TT = zeros(2,2);
    k.Tg = zeros(2,1);
    
    k.gs = 0;
    k.gT = zeros(1,2);
    k.gg = 0;
    
    %Gauss integration
    ngp = 2;
    [W,xi] = gaussQuad(ngp);
    get_properties      %retrieve material properties
    for i = 1:ngp
        [N,B] = get_shape_functions(xi(i),h);
        T_xi = dot(T,N);
        
        %Mass matrix
        m.v = m.v+W(i)*rho*(N'*N)*J;
        m.T = m.T+W(i)*rho*cp*(N'*N)*J;
        m.s = m.s+W(i)*J;
        m.g = m.g+W(i)*J;
        
        %Linear Stiffness matrix
        k.vs = k.vs-W(i)*B'*J;
        k.sv = k.sv+W(i)*E*B*J;
        k.TT = k.TT-W(i)*lambda*(B'*B)*J;
        
        %Non linear stiffness matrix
        [g,dgdp,dgds,dgdT] = get_plastic_strain_rate(s,T_xi,p);
        k.ss = k.ss-W(i)*E*dgds*J;
        k.sT = k.sT-W(i)*E*dgdT*N*J;
        k.sg = k.sg-W(i)*E*dgdp*J;
        
        [F1,F2] = get_F(s,g,dgds);
        k.Ts = k.Ts+W(i)*N'*F1*s*J;
        k.TT = k.TT+W(i)*(N'*N)*F2*s*dgdT*J;
        k.Tg = k.Tg+W(i)*N'*F2*s*dgdp*J;
        
        k.gs = k.gs+W(i)*dgds*J;
        k.gT = k.gT+W(i)*dgdT*N*J;
        k.gg = k.gg+W(i)*dgdp*J;
    end
end

function [B,N] = get_shape_functions(xi,h)
    %For velocity and temperature
    N = (1/2)*[1-xi 1+xi];
    B = (1/h)*[-1, 1];
end

function [g,dgdp,dgds,dgdT] = get_plastic_strain_rate(s,T,p)        
    %Johnson-Cook Model
    get_parameters          %get parameters for the constitutive model
    
    P = 1-((T-To)/(Tm-To))^m;
    Q = A+B*p^N;
    
    g = pdotr*e^((s/(P*Q)-1)/c);
    dgds = g/(c*P*Q);
    dgdp = -dgds*B*N*p^(N-1)*s/Q;
    dgdT = dgds*m*s*((T-To)/(Tm-To))^(m-1)/(P*(Tm-To));
end

function [F1,F2] = get_F(s,g,dgds)
    %For constant Taylor Quinney Coefficient
    n = 1;
    dndg = 0;
    get_parameters          %get parameters for the constitutive model
    
    F1 = chi*(dndg*dgds*g+n*g/s+n*dgds);
    F2 = chi*(dndg*g+n);
end