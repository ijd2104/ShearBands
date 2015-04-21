function [r,j,f] = get_element_stiffness(x0,x,h,f0,ue,flg) 
% Inputs
% x0  -  x at the previous time step
% x   -  x at the current Newton iteration
% h   -  element size
% f0  -  forcings from the previous Newton iteration
% flg -  flag
%        3 for residual only
%        5 for jacobian only
%        6 for jacobian + residual

% Outputs
% r   -  element residual
% j   -  element jacobian
% f   -  element forcings at the current Newton iteration

    %% Initialization
    global t TimeIntPar matProp v_BC
    
    v0 = x0(1:2);
    s0 = x0(3);
    T0 = x0(4:5);
    p0 = x0(6);
    
    fs0 = f0(3);
    fT0 = f0(4:5);
    fg0 = f0(6);

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
    
    %Forcings
    f.s = 0;
    f.T = zeros(2,1);
    f.g =0;
    
    %Non linear terms
    G.ss = 0;
    G.sT = zeros(1,2);
    G.sg = 0;
    
    G.Ts = zeros(2,1);
    G.TT = zeros(2,2);
    G.Tg = zeros(2,1);
    
    G.gs = 0;
    G.gT = zeros(1,2);
    G.gg = 0;
    
    %Load material properties
    rho = matProp.rho;
    cp = matProp.cp;
    lambda = matProp.lambda;
    chi = matProp.chi;
    E = matProp.G;
    
    %% Gauss integration
    ngp = 2;
    [W,xi] = gaussian_quadrature(ngp);
    J = h/2;
    
    for i = 1:ngp
        [Nv,Ns,NT,Ng,Bv,BT] = get_shape_functions(xi(i),h);
        T_xi = dot(T,NT);
        x_xi = dot(ue,Nv);
        nimp = set_imperfection(x_xi);
        
        [g,dgdp,dgds,dgdT] = get_plastic_strain_rate(s,T_xi,p,nimp);
        [F1,F2] = get_F(s,g,dgds);
        
        %Mass matrix
        m.v = m.v+W(i)*rho*(Nv'*Nv)*J;
        m.s = m.s+W(i)*(Ns'*Ns)*J;
        m.T = m.T+W(i)*rho*cp*(NT'*NT)*J;
        m.g = m.g+W(i)*(Ng'*Ng)*J;
        
        %Linear Stiffness matrix
        k.vs = k.vs-W(i)*(Bv'*Ns)*J;
        k.sv = k.sv+W(i)*E*(Ns'*Bv)*J;
        k.TT = k.TT-W(i)*lambda*(BT'*BT)*J;
        
        if flg==3 || flg==6 %residual+jacobian
            %Forcings
            f.s = f.s-W(i)*E*g*Ns'*J;
            f.T = f.T+W(i)*chi*s*g*NT'*J;
            f.g = f.g+W(i)*g*Ng'*J;
            if flg==3 %residual only
                continue
            end
        end
        
        %Non linear stiffness matrix
        G.ss = G.ss-W(i)*E*dgds*(Ns'*Ns)*J;
        G.sT = G.sT-W(i)*E*dgdT*(Ns'*NT)*J;
        G.sg = G.sg-W(i)*E*dgdp*(Ns'*Ng)*J;
        
        G.Ts = G.Ts+W(i)*F1*s*(NT'*Ns)*J;
        G.TT = G.TT+W(i)*F2*s*dgdT*(NT'*NT)*J;
        G.Tg = G.Tg+W(i)*F2*s*dgdp*(NT'*Ng)*J;
        
        G.gs = G.gs+W(i)*dgds*(Ng'*Ns)*J;
        G.gT = G.gT+W(i)*dgdT*(Ng'*NT)*J;
        G.gg = G.gg+W(i)*dgdp*(Ng'*Ng)*J;
    end
    
    %% Compute residual
    if flg==3 || flg==6 %residual / residual+jacobian
        dt = t.dt;
        a = TimeIntPar.alpha;
        
        vdot = (v-v0)/dt;
        sdot = (s-s0)/dt;
        Tdot = (T-T0)/dt;
        pdot = (p-p0)/dt;
        
        r.v = m.v*vdot-(1-a)*(k.vs*s0)-a*(k.vs*s);
        r.s = m.s*sdot-(1-a)*(k.sv*v0+fs0)-a*(k.sv*v+f.s);
        r.T = m.T*Tdot-(1-a)*(k.TT*T0+fT0)-a*(k.TT*T+f.T);
        r.g = m.g*pdot-(1-a)*(fg0)-a*(f.g);
    elseif flg==5 %jacobian only
        r = [];
        f = f0;
    end
    %% Compute analytical Jacobian
    if flg==5 || flg==6 %jacobian / residual+jacobian
        j.vv = m.v/dt;
        j.vs = -a*k.vs;
        %j.vT = 0;
        %j.vg = 0;
        
        j.sv = -a*k.sv;
        j.ss = m.s/dt-a*G.ss;
        j.sT = -a*G.sT;
        j.sg = -a*G.sg;
        
        %j.Tv = 0;
        j.Ts = -a*G.Ts;
        j.TT = m.T/dt-a*(k.TT+G.TT);
        j.Tg = -a*G.Tg;
        
        %j.gv = 0;
        j.gs = -a*G.gs;
        j.gT = -a*G.gT;
        j.gg =  m.g/dt-a*G.gg;
    elseif flg==3 %residual only
        j = [];
    end
    
    %% Normalization
    normv = v_BC; %maximum vel (equivalent to 15)
    norms = 457.3E6; %yield shear stress
    normT = 298; %reference temperature
    normg = 1E-8;%457.3E6/200E9; %yield strain
    
    r.v = r.v/normv;
    r.s = r.s/norms;
    r.T = r.T/normT;
    r.g = r.g/normg;
    
    j.vv = j.vv/normv;
    j.vs = j.vs/normv;
    
    j.sv = j.sv/norms;
    j.ss = j.ss/norms;
    j.sT = j.sT/norms;
    j.sg = j.sg/norms;
    
    j.Ts = j.Ts/normT;
    j.TT = j.TT/normT;
    j.Tg = j.Tg/normT;
    
    j.gs = j.gs/normg;
    j.gT = j.gT/normg;
    j.gg = j.gg/normg;
    
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

function [g,dgdp,dgds,dgdT] = get_plastic_strain_rate(s,T,p,nimp)
    %Load constitutive model parameters
    global modelPar
    A = modelPar.A*nimp;
    B = modelPar.B*nimp;
    N = modelPar.N;
    To = modelPar.To;
    Tm = modelPar.Tm;
    m = modelPar.m;
    c = modelPar.c;
    pdotr = modelPar.pdotr;
    
    P = 1-((T-To)/(Tm-To))^m;
    Q = A+B*p^N;
    
    g = pdotr*exp((s/(P*Q)-1)/c);
    dgds = g/(c*P*Q);
    dgdp = -dgds*B*N*p^(N-1)*s/Q;
    dgdT = dgds*m*s*((T-To)/(Tm-To))^(m-1)/(P*(Tm-To));
    
    if g <= 1E-16
        %True for the first iteration
        g = 0;
        dgds = 0;
        dgdp = 0;
        dgdT = 0;
    elseif p == 0
        dgdp = 0;
    elseif isnan(g)
        %Stress s is NaN, P is NaN, Q is huge
        error('g is NaN')
    end
    if isinf(dgds)
        error('dgds is infinite')
    end
    
end

function [F1,F2] = get_F(s,g,dgds)
    global matProp
    chi = matProp.chi;
    %For constant Taylor Quinney Coefficient
    n = 1;
    dndg = 0;
    
    F1 = chi*(dndg*dgds*g+n*g/s+n*dgds);
    %F1 becomes NaN because stress is zero, g is zero
    if s == 0
        F1 = 0;
    end
    F2 = chi*(dndg*g+n);
end

function nimp = set_imperfection(x)
    e = 1E1;
    ared = 0.01;
    r0 = 10E-6;
    x0 = 1E-3/2;
    rn = abs(x-x0)/r0;
    nimp = 1-ared*(2/(e^rn+e^(-rn)));
end