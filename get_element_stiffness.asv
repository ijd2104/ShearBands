function [r,m,k] = get_element_stiffness( xe,flg )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Linear terms
    B = (1/h)*[-1, 1];
    B_T = B';
    
    m_CG = zeros(2,2);
    m_DG = 0;
    kvs = zeros(2,1);
    ksv = zeros(1,2);
    kTT = zeros(2,2); 
    
    ngp = 2; %number of Gauss points
    [W,xi] = gaussQuad(ngp);
    for i = 1:ngp
        N = (1/2)*[1-xi(i) 1+xi(i)];
        N_T = N';
        
        
        m_CG = m_CG+W(i)*(N_T*N);
        m_DG = m_DG+W(i);
        
        kvs = kvs+W(i)*B_T;
        ksv = ksv+W(i)*V;
        kTT = kTT+W(i)*(B_T*B);
        
        
    end
    mv = rho*m_CG*J;
    mT = rho*c*m_CG*J;
    ms = m_DG*J;
    mg = m_DG*J;
    kvs = -kvs*J;
    ksv = E*ksv*J;
    kTT = -lambda*kTT*J;
    
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
    B = (1/h)*[-1, 1];
    B_T = B';
    for i = 1:ngp
        N = (1/2)*[1-xi(i) 1+xi(i)];
        N_T = N';
        
        gg = get_gg
        gs = get_gs
        gT = get_gT
        
        Gss = Gss+W(i)*gs;
        GsT = GsT+W(i)*gT*N;
        Gsg = Gsg+W(i)*gg;
        
        GTs = GTs+W(i)*N_T;
        GTT = GTT+W(i)*(N_T*N)*gT;
        GTg = GTg+W(i)*N_T*gg;
        
        Ggs = Ggs+W(i)*gs;
        GgT = GgT+W(i)*gT*N;
        Ggg = Ggg+W(i)*gg;
    end

%% Other
    
end

