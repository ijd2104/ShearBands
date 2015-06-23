function coupled_shear_bands()
%LINEARELASTICITY Solves linear elasticity problem using mixed FEM
%   Saves velocity and displacement at intervals
    global t N

    setupData();
    getMesh();
    
    tol = 1E-3;
    e = floor(N.elem/2);
    xhist = zeros(1,4);
    flgRL = 1;

    X0 = zeros(N.node,1);
    
    %Initialize temperature
    T1 = N.vnode+N.snode+1;
    Tn = N.vnode+N.snode+N.Tnode;
    X0(T1:Tn) = 310;
    strs = 0;

    Xt = X0;
    for n = 2:t.steps
        X0 = Xt;
        t.iter = n;
        t.curr = t.curr+t.dt;
        Xt = newtonIter(X0);
        
        [xe,he,ue] = getx(Xt,e);
        xhist(1:3) = xhist(2:4);
        xhist(4) = xe;
        Me = get_element_stiffness(xhist(3),xhist(4),he,[],ue,8);
        Ke = get_element_stiffness(xhist(3),xhist(4),he,[],ue,9);
        [w,v] = myeig(Ke,Me);
        if real(w(1)) > tol && flgRL
            compute_root_locus(xhist(1),xhist(2),he,w,v);
            flgRL = 0;
            error('Program terminated')
        end
        
        if Xt(N.vnode+1)<0.8*strs
            break
        elseif Xt(N.vnode+1)>strs
            strs = Xt(N.vnode+1);
        end
    end
end

function setupData()
    global t TimeIntPar NewtonPar matProp modelPar v_BC
    
    %Time data
    t.total = 2.5E-4;
    t.steps = 500;
    t.dt = t.total/t.steps;
    if t.dt<5E-11
        pause
    end
    t.ramp = t.total/2.5;
    t.curr = 0;
    
    %Initial velocity
    v_BC = 5E0;
    
    %Newton's method parameters
    NewtonPar.NormTol = 1E-5;
    NewtonPar.iter = 25;
    
    %Time integration parameters
    TimeIntPar.alpha = 1;
    
    %Material properties for 4340 Steel
    matProp.rho = 7830;             %kg/m^3 mass density
    matProp.cp = 477;               %J/kgK specific heat
    matProp.lambda = 38;            %W/mK thermal conductivity
    matProp.E = 200E9;              %Pa Elastic modulus
    matProp.mu = 0.29;              % Poisson ratio
    matProp.G = 77.5E9;             %Pa shear modulus
    matProp.chi = 0.9;              % Taylor-Quinney coefficient
    
    %Johnson Cook Parameters for 4340 Steel
    modelPar.A = 457.3E6;            %Pa yield shear stress
    modelPar.B = 294.4E6;            %Pa shear stress hardening parameter
    modelPar.N = 0.26;               % strain hardening parameter
    modelPar.To = 298;               %K reference temperature
    modelPar.Tm = 1793;              %K melting temperature
    modelPar.m = 1.03;               % thermal softening exponent
    modelPar.c = 0.014;              % strain-rate hardening parameter
    modelPar.pdotr = 1.0;            %1/s reference strain rate
end

function getMesh()
    global N u
    
    L = 1E-3;
  
    N.elem = 400;
    N.vnode = N.elem+1;
    N.Tnode = N.elem+1;
    N.snode = N.elem;
    N.gnode = N.elem;
    N.conn = [1:N.elem; 2:N.vnode]';
    N.node = N.vnode+N.Tnode+N.snode+N.gnode;
    
    de = 2E-7;
    B = 0.25;
    P = 5;
    nnode = round(N.vnode/2);
    S0 = de/(L/2)*nnode;
    
    u = zeros(1,nnode);
    for i = 1:nnode
        xrel = i/nnode;
        if xrel<=B
            u(i) = S0*xrel;
        elseif xrel>B
            u(i) = S0*xrel+(1-S0)*((xrel-B)/(1-B))^P;
        end
    end
    u_flp = -fliplr(u);
    if mod(nnode,2)
        u = [u_flp(1:end-1) u];
    else
        u = [u_flp u];
    end
    u = u/2E3; %position of nodes
end

function Xk = newtonIter(Xt)
    global NewtonPar N
    
    niter = NewtonPar.iter;
    ntol = NewtonPar.NormTol;
    
    Xk = Xt;
    F0 = [zeros(N.vnode,1);  %velocity
         zeros(N.snode,1);  %stress
         zeros(N.Tnode,1);  %temperature
         zeros(N.gnode,1)]; %plastic strain
    
    [J,R,F] = matrixAssembly(Xt,Xk,F0);
    [J,R,X] = applyBC(J,R,Xk);
    
    velBC = get_vbc(X);
    dX.D = -X.D+velBC;
    R.D = -dX.D;
    
    dX.N = -J.NN\(R.N+J.ND*dX.D);
    Xk = [X.D+dX.D; X.N+dX.N];
    Xk = unPartition(Xk);
    
    if ~isreal(Xk)
        error('Xk became imaginary')
    end
    
    for k = 1:niter
        [J,R,F] = matrixAssembly(Xt,Xk,F);
        [J,R,X] = applyBC(J,R,Xk);
        
        if norm(R.N) < ntol
            break
        end
        
        if k==20
            error('Exceeded 20 Newton iterations')
        end
        
        dX.N = -J.NN\R.N;
        
        Xk = [X.D; X.N+dX.N];
        
        Xk = unPartition(Xk);
        
        if ~isreal(Xk)
            error('Xk became imaginary')
        end
    end
end

function [J,R,F] = matrixAssembly(Xt,Xn,F0)
    global N u
    
    %Initialize matrices
    R.v = zeros(N.vnode,1);
    R.s = zeros(N.snode,1);
    R.T = zeros(N.Tnode,1);
    R.g = zeros(N.gnode,1);
    
    F.v = zeros(N.vnode,1);
    F.s = zeros(N.snode,1);
    F.T = zeros(N.Tnode,1);
    F.g = zeros(N.gnode,1);
    
    J.vv = zeros(N.vnode,N.vnode);
    J.vs = zeros(N.vnode,N.snode);
    J.vT = zeros(N.vnode,N.Tnode); %
    J.vg = zeros(N.vnode,N.gnode); %
    
    J.sv = zeros(N.snode,N.vnode);
    J.ss = zeros(N.snode,N.snode);
    J.sT = zeros(N.snode,N.Tnode);
    J.sg = zeros(N.snode,N.gnode);
    
    J.Tv = zeros(N.Tnode,N.vnode); %
    J.Ts = zeros(N.Tnode,N.snode);
    J.TT = zeros(N.Tnode,N.Tnode);
    J.Tg = zeros(N.Tnode,N.gnode);
    
    J.gv = zeros(N.gnode,N.vnode); %
    J.gs = zeros(N.gnode,N.snode);
    J.gT = zeros(N.gnode,N.Tnode);
    J.gg = zeros(N.gnode,N.gnode);
    
    %Assemble matrices
    for e = 1:N.elem
        
        L = N.conn(e,:);
        ue = u(L);
        h = abs(ue(2)-ue(1));
        
        Lv = L;
        Ls = e+N.vnode;
        LT = Lv+N.vnode+N.elem;
        Lg = Ls+N.vnode+N.elem;
        LL = [Lv Ls LT Lg];
        
        x0 = Xt(LL);
        x  = Xn(LL);
        f0 = F0(LL);
        
        [r,j,f] = get_element_stiffness(x0,x,h,f0,ue,6);
        
        % Residual vector
        R.v(L) = R.v(L)+r.v;
        R.s(e) = R.s(e)+r.s;
        R.T(L) = R.T(L)+r.T;
        R.g(e) = R.g(e)+r.g;
        
        %Forcing vector
        F.s(e) = F.s(e)+f.s;
        F.T(L) = F.T(L)+f.T;
        F.g(e) = F.g(e)+f.g;
        
        % Jacobian matrix
        J.vv(L,L) = J.vv(L,L)+j.vv;
        J.vs(L,e) = J.vs(L,e)+j.vs;
        %J.vT
        %J.vg
        
        J.sv(e,L) = J.sv(e,L)+j.sv;
        J.ss(e,e) = J.ss(e,e)+j.ss;
        J.sT(e,L) = J.sT(e,L)+j.sT;
        J.sg(e,e) = J.sg(e,e)+j.sg;
        
        %J.Tv
        J.Ts(L,e) = J.Ts(L,e)+j.Ts;
        J.TT(L,L) = J.TT(L,L)+j.TT;
        J.Tg(L,e) = J.Tg(L,e)+j.Tg;
        
        %J.gv
        J.gs(e,e) = J.gs(e,e)+j.gs;
        J.gT(e,L) = J.gT(e,L)+j.gT;
        J.gg(e,e) = J.gg(e,e)+j.gg;
    end
    
    R = [R.v; R.s; R.T; R.g];
    F = [F.v; F.s; F.T; F.g];
    
    JJ = [J.vv J.vs J.vT J.vg;
          J.sv J.ss J.sT J.sg;
          J.Tv J.Ts J.TT J.Tg;
          J.gv J.gs J.gT J.gg];
    
    J = sparse(JJ);
end

function [J,R,X] = applyBC(J0,R0,X0)
    global N BC eN eD
    
    BC = [1 -1; N.vnode 1];
    nBC=size(BC,1);
    eN = [1:N.node]';
    eD = BC(:,1);
    for i = 1:nBC
        eN(eN==eD(i)) = [];
    end    
    
    X.D = X0(eD);
    X.N = X0(eN);
    J.ND = J0(eN,eD);
    J.NN = J0(eN,eN);
    J.DN = 0*J0(eD,eN);
    J.DD = eye(nBC);
    R.N = R0(eN);
    R.D = R0(eD);
    
    X.BCp = zeros(size(BC,1),1);
    for i = 1:size(BC,1)
        X.BCp(i) = BC(i,2);
    end
end

function X = unPartition(X0)
    global eN eD
    e = [eD;eN];

    X = zeros(numel(e),1);
    for i = 1:numel(e)
        X(e(i)) = X0(i);
    end
end

function v = get_vbc(X)
    vt = get_vel();
    v = X.BCp*vt;
end

function v = get_vel()
    global t v_BC
    
    a = t.curr/t.ramp;
    if a < 1
        f = a;
    else
        f = 1;
    end
    v = v_BC*f;
end

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
%        8 Me
%        9 Ke

% Outputs
% r   -  element residual
% j   -  element jacobian
% f   -  element forcings at the current Newton iteration

    %% Initialization
    global t TimeIntPar matProp modelPar v_BC
    
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
        %j.vT = zeros(2,2);
        %j.vg = zeros(2,1);
        
        j.sv = -a*k.sv;
        j.ss = m.s/dt-a*G.ss;
        j.sT = -a*G.sT;
        j.sg = -a*G.sg;
        
        %j.Tv = zeros(2,2);
        j.Ts = -a*G.Ts;
        j.TT = m.T/dt-a*(k.TT+G.TT);
        j.Tg = -a*G.Tg;
        
        %j.gv = zeros(1,2);
        j.gs = -a*G.gs;
        j.gT = -a*G.gT;
        j.gg =  m.g/dt-a*G.gg;
    elseif flg==3 %residual only
        j = [];
    end
    
    %% Normalization
    normv = v_BC;           %velocity at the boundary condition
    norms = modelPar.A;     %yield shear stress
    normT = modelPar.To;    %reference temperature
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
    ared = 0.01;
    r0 = 10E-6;
    x0 = 1E-3/2;
    rn = abs(x-x0)/r0;
    nimp = 1-ared*(2/(exp(rn)+exp(-rn)));
end


function [x,h,ue] = getx(X,e)
    global u N
    x1 = 6*(e-1)+1;
    xn = x1+5;
    x = X(x1:xn);
    
    L = N.conn(e,:);
    ue = u(L);
    h = abs(ue(2)-ue(1));
end

function [w,v] = myeig(k,m)
    [v,w] = eig(k,m,'vector');
    [w,i] = sort(real(w),'ascend');
    v = v(:,i);
end

function compute_root_locus(x0,xm,h,w,v)
