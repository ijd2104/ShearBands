function [r,j,f] = get_element_stiffness(x0,x,h,f0) 
    %% Initialization
    global t TimeIntPar matProp
    
    v0 = x0(1:2);
    s0 = x0(3);
    T0 = x0(4:5);
    p0 = x0(6);
    
    fs0 = f0(3);
    fT0 = f0(4:5);
    fg0 = f0(6);

    v = x(1:2);
    s = x(3);
    if s > 1E6
        %pause
    end
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
    chi = matProp.chi;
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
    
    f.s = 0;
    f.T = zeros(2,1);
    f.g =0;
    ngp = 2;
    [W,xi] = gaussian_quadrature(ngp);
    for i = 1:ngp
        [g,dgdp,dgds,dgdT] = get_plastic_strain_rate(s,T_xi,p);
        [Nv,Ns,NT,Ng,Bv,BT] = get_shape_functions(xi(i),h);
        f.s = f.s-W(i)*G*g*(Ns'*Ns)*J;
        f.T = f.T+W(i)*chi*s*g*(NT'*Ns)*J;
        f.g = f.g+W(i)*g*(Ng'*Ng)*J;
    end
    
    r.v = m.v*vdot-(1-a)*(k.vs*s0)-a*(k.vs*s);
    r.s = m.s*sdot-(1-a)*(k.sv*v0+fs0)-a*(k.sv*v+f.s);
    r.T = m.T*Tdot-(1-a)*(k.Ts*s0+fT0)-a*(k.Ts*s+f.T);
    r.g = m.g*pdot-(1-a)*(fg0)-a*(f.g);
    
    %% Compute analytical Jacobian
    j.vv = m.v/dt;
    if isnan(sum(sum(j.vv))) || isinf(sum(sum(j.vv)))
        %j.vv
    end
    j.vs = -a*k.vs;
    if isnan(sum(sum(j.vs))) || isinf(sum(sum(j.vs)))
        %j.vs
    end
    %j.vT = 0;
    %j.vg = 0;
    
    j.sv = -a*k.sv;
    if isnan(sum(sum(j.sv))) || isinf(sum(sum(j.sv)))
        %j.sv
    end
    j.ss = m.s/dt-a*k.ss;
    if isnan(sum(sum(j.ss))) || isinf(sum(sum(j.ss)))
        %j.ss
    end
    
    j.sT = -a*k.sT;
    if isnan(sum(sum(j.sT))) || isinf(sum(sum(j.sT)))
        %j.sT
    end
    j.sg = -a*k.sg;
    if isnan(sum(sum(j.sg))) || isinf(sum(sum(j.sg)))
        %j.sg
    end
    
    %j.Tv = 0;
    j.Ts = -a*k.Ts;
    if isnan(sum(sum(j.Ts))) || isinf(sum(sum(j.Ts)))
        %j.Ts
    end
    j.TT = m.T/dt-a*k.TT;
    if isnan(sum(sum(j.TT))) || isinf(sum(sum(j.TT)))
        %j.TT
    end
    j.Tg = -a*k.Tg;
    if isnan(sum(sum(j.Tg))) || isinf(sum(sum(j.Tg)))
        %j.Tg
    end
%     
    %j.gv = 0;
    j.gs = -a*k.gs;
    if isnan(sum(sum(j.gs))) || isinf(sum(sum(j.gs)))
        %j.gs
    end
     j.gT = -a*k.gT;
    if isnan(sum(sum(j.gT))) || isinf(sum(sum(j.gT)))
        %j.gT
    end
    j.gg =  m.g/dt-a*k.gg;
    if isnan(sum(sum(j.gg))) || isinf(sum(sum(j.gg)))
        %j.gg
    end
    
    %% Normalization
    normv = 1E4; %wild guess
    norms = 457.3E6; %yield shear stress
    normT = 298; %reference temperature
    normg = 457.3E6/200E9; %yield strain
    
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
    
    g = pdotr*exp((s/(P*Q)-1)/c);
    dgds = g/(c*P*Q);
    dgdp = -dgds*B*N*p^(N-1)*s/Q;
    %if p is zero, dgdp becomes NaN
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
        error('g is undefined')
    end
    
    if isinf(dgds)
        %dgds
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