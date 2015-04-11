function [r,j] = get_element_stiffness(x0,x,h) 
    %% Initialization
    global t TimeIntPar matProp
    
    v0 = x0(1:2);
    s0 = x0(3);
    T0 = x0(4:5);
    p0 = x0(6);

    v = x(1:2);
    s = x(3);
    T = x(4:5);
    p = x(6);

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
    %k.TT = zeros(2,2);
    k.Tg = zeros(2,1);
    
    k.gs = 0;
    k.gT = zeros(1,2);
    k.gg = 0;
    
    %Load material properties
    rho = matProp.rho;
    cp = matProp.cp;
    lambda = matProp.lambda;
    %E = matProp.E;
    G = matProp.G;
    
    %% Gauss integration
    ngp = 2;
    [W,xi] = gaussian_quadrature(ngp);
    J = h/2;
    
    for i = 1:ngp
        [Nv,Ns,NT,Ng,Bv,BT] = get_shape_functions(xi(i),h);
        T_xi = dot(T,NT);
        
        %Mass matrix
        m.v = m.v+W(i)*rho*(Nv'*Nv)*J;
        m.s = m.s+W(i)*(Ns'*Ns)*J;
        m.T = m.T+W(i)*rho*cp*(NT'*NT)*J;
        m.g = m.g+W(i)*(Ng'*Ng)*J;
        
        %Linear Stiffness matrix
        k.vs = k.vs-W(i)*(Bv'*Ns)*J;
        k.sv = k.sv+W(i)*G*(Ns'*Bv)*J;
        k.TT = k.TT-W(i)*lambda*(BT'*BT)*J;
        
        %Non linear stiffness matrix
        [g,dgdp,dgds,dgdT] = get_plastic_strain_rate(s,T_xi,p);
        k.ss = k.ss-W(i)*G*dgds*(Ns'*Ns)*J;
        k.sT = k.sT-W(i)*G*dgdT*(Ns'*NT)*J;
        k.sg = k.sg-W(i)*G*dgdp*(Ns'*Ng)*J;
        
        [F1,F2] = get_F(s,g,dgds);
        k.Ts = k.Ts+W(i)*F1*s*(NT'*Ns)*J;
        k.TT = k.TT+W(i)*F2*s*dgdT*(NT'*NT)*J;
        k.Tg = k.Tg+W(i)*F2*s*dgdp*(NT'*Ng)*J;
        
        k.gs = k.gs+W(i)*dgds*(Ng'*Ns)*J;
        k.gT = k.gT+W(i)*dgdT*(Ng'*NT)*J;
        k.gg = k.gg+W(i)*dgdp*(Ng'*Ng)*J;
    end
    
    %% Compute residual
    dt = t.dt;
    a = TimeIntPar.alpha;
    
    vdot = (v-v0)/dt;
    sdot = (s-s0)/dt;
    Tdot = (T-T0)/dt;
    pdot = (p-p0)/dt;
    
    fs = 0;
    fT = zeros(2,1);
    fg =0;
    
    ngp = 1;
    [W,xi] = gaussian_quadrature(ngp);
    for i = 1:ngp
        [Nv,Ns,NT,Ng,Bv,BT] = get_shape_functions(xi(i),h);
        fs = fs+W(i)*E*g*J;
        fT = fT+W(i)*N'*s*chi*g*J;
        fg = fg+W(i)*g*J;
    end
    
    r.v = m.v*vdot-(1-a)*(k.vs*s0)-a*(k.vs*s);
    r.T = m.s*sdot-(1-a)*(k.sv*v0+k.ss*s0+k.sT*T0+k.sg*p0)-a*(k.sv*v+k.ss*s+k.sT*T+k.sg*p);
    r.s = m.T*Tdot-(1-a)*(k.Ts*s0+k.TT*T0+k.Tg*p0)-a*(k.Ts*s+k.TT*T+k.Tg*p);
    r.g = m.g*pdot-(1-a)*(k.gs*s0+k.gT*T0+k.gg*p0)-a*(k.gs*s+k.gT*T+k.gg*p);
    
    %% Compute analytical Jacobian
    j.vv = m.v/dt;
    j.vs = -a*k.vs;
    %j.vT = 0;
    %j.vg = 0;
    
    j.sv = -a*k.sv;
    j.ss = m.s/dt-a*k.ss;
    j.sT = -a*k.sT;
    j.sg = -a*k.sg;
    
    %j.Tv = 0;
    j.Ts = -a*k.Ts;
    j.TT = m.T/dt-a*k.TT;
    j.Tg = -a*k.Tg;
    
    %j.gv = 0;
    j.gs = -a*k.gs;
    j.gT = -a*k.gT;
    j.gg =  m.g/dt-a*k.gg;
end

function [Nv,Ns,NT,Ng,Bv,BT] = get_shape_functions(xi,h)
    %For velocity
    Nv = (1/2)*[1-xi 1+xi];
    Bv = (1/h)*[-1, 1];
    
    %For stress
    Ns = 1;
    %Bs = 0;
    
    %For temperature
    NT = (1/2)*[1-xi 1+xi];
    BT = (1/h)*[-1, 1];
    
    %For plastic strain
    Ng = 1;
    %Bg = 0;
end

function [g,dgdp,dgds,dgdT] = get_plastic_strain_rate(s,T,p)
    %Load constitutive model parameters
    global modelPar
    A = modelPar.A;
    B = modelPar.B;
    N = modelPar.N;
    To = modelPar.To;
    Tm = modelPar.Tm;
    m = modelPar.m;
    c = modelPar.c;
    pdotr = modelPar.pdotr;
    
    P = 1-((T-To)/(Tm-To))^m;
    Q = A+B*p^N;
    
    g = pdotr*e^((s/(P*Q)-1)/c);
    dgds = g/(c*P*Q);
    dgdp = -dgds*B*N*p^(N-1)*s/Q;
    dgdT = dgds*m*s*((T-To)/(Tm-To))^(m-1)/(P*(Tm-To));
end

function [F1,F2] = get_F(s,g,dgds)
    global matProp
    chi = matProp.chi;
    %For constant Taylor Quinney Coefficient
    n = 1;
    dndg = 0;
    
    F1 = chi*(dndg*dgds*g+n*g/s+n*dgds);
    F2 = chi*(dndg*g+n);
end