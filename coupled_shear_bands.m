function X = coupled_shear_bands()
%LINEARELASTICITY Solves linear elasticity problem using mixed FEM
%   Saves velocity and displacement at intervals
    global t N
    
    setupData();
    getMesh();
    
    i = 10; %save interval
    X = num2cell(zeros(N.node,t.steps/i),1);
    X0 = zeros(N.node,1);
    
    Xt = X0;
    for n = 2:t.steps
        t.iter = n;
        t.curr = t.curr+t.dt;
        Xt = newtonIter(Xt);
        if mod(n,10)==0
            X{n/i} = Xt;
        end
    end
    X = [X0 cell2mat(X)];
end

function setupData()
    global t TimeIntPar NewtonPar matProp modelPar
    
    %Time data
    t.total = 2/5048; %Time for wave to travel 1m, 1/5000
    t.steps = 1000;
    t.dt = t.total/t.steps;
    t.ramp = t.total/10;
    t.curr = 0;
    
    %Newton's method parameters
    NewtonPar.NormTol = 1E-5;
    NewtonPar.iter = 20;
    
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
    
    L = 1;
    N.elem = 100;
    u = 0:(1/N.elem):L; %position of nodes
    N.vnode = numel(u);
    N.eqn = N.vnode;
    N.conn = [1:N.elem; 2:N.vnode]';
    N.node = N.vnode+N.elem;
end

function Xk = newtonIter(Xt)
    global NewtonPar
    
    niter = NewtonPar.iter;
    ntol = NewtonPar.NormTol;
    
    Xk = Xt;
    [J,R] = matrixAssembly(Xt,Xk);
    %[J,R,X] = matrixPartition(J,R,Xk);
    [J,R,X] = applyBC(J,R,Xk);
    
    velBC = get_vbc(X);
    dX.D = -X.D+velBC;
    R.D = -dX.D;
    dX.N = -J.NN\(R.N+J.ND*dX.D);
    Xk = [X.D+dX.D; X.N+dX.N];
    Xk = unPartition(Xk);
    
    for k = 1:niter
        %fprintf('Newton iteration %d\n',k);
        [J,R] = matrixAssembly(Xt,Xk);
        [J,R,X] = applyBC(J,R,Xk);

        %norm(R.N)
        if norm(R.N) < ntol
            break
        end
        if k==20,
            error('stop')
        end
        
        dX.N = -J.NN\R.N;
        
        Xk = [X.D; X.N+dX.N];
        
        Xk = unPartition(Xk);
    end
        figure(1)
        plot(Xk(1:101))
        axis([0 101 -8E5 8E5])
        xlabel('Nodes')
        ylabel('Velocity m/s')
        pause(.001)
end

function [J,R] = matrixAssembly(Xt,Xn)
    global N u
    
    %Initialize matrices
    R.v = zeros(N.vnode,1);
    R.s = zeros(N.elem,1);
    R.T = zeros(N.vnode,1);
    R.g = zeros(N.elem,1);
    
    J.vv = zeros(N.vnode,N.vnode);
    J.vs = zeros(N.vnode,N.elem);
    %J.vT = zeros(N.vnode,N.vnode);
    %J.vg = zeros(N.vnode,N.elem);
    
    J.sv = zeros(N.elem,N.vnode);
    J.ss = zeros(N.elem,N.elem);
    J.sT = zeros(N.elem,N.vnode);
    J.sg = zeros(N.elem,N.elem);
    
    %J.Tv = J.vv;
    J.Ts = J.vs;
    J.TT = J.vv;
    J.Tg = J.vg;
    
    %J.gv = J.sv;
    J.gs = J.ss;
    J.gT = J.sv;
    J.gg = J.ss;
    
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
        
        [r,j] = get_element_stiffness(x0,x,h);
        
        % Residual vector
        R.v(L) = R.v(L)+r.v;
        R.s(e)  = R.s(e)+r.s;
        R.T(L) = R.T(L)+r.T;
        R.g(e)  = R.g(e)+r.g;
        
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
    
    R = [R.v; R.s; R.T R.g];
    
    JJ = blkdiag([J.vv J.vs; J.sv J.ss],[J.TT J.Tg; J.gT J.gg]);
    JJ(Ls,LT) = J.sT;
    JJ(Ls,Lg) = J.sg;
    JJ(LT,Ls) = J.Ts;
    JJ(Lg,Ls) = J.gs;
    
    J = sparse(JJ);
end

function [J,R,X] = applyBC(J0,R0,X0)
    global N BC eN eD
    
    BC = [1 0; N.vnode 1];
    nBC=size(BC,1);
    eN = [1:N.node]';
    eD = BC(:,1);
    for i = 1:nBC
        eN(eN==eD(i)) = [];
    end    
    
    R0(1:N.vnode) = R0(1:N.vnode)/1E4;
    R0(N.vnode+1:end) = R0(N.vnode+1:end)/250E9;
    J0(1:N.vnode,:) = J0(1:N.vnode,:)/1E4;
    J0(N.vnode+1:end,:)=J0(N.vnode+1:end,:)/250E9;
    
    X.D = X0(eD);
    X.N = X0(eN);
    J.ND = J0(eN,eD);
    J.NN = J0(eN,eN);
    J.DN = 0.0.*J0(eD,eN);
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
    global t
    v0 = 10;
    
    a = t.curr/t.ramp;
    if a < 1
        f = (pi/t.ramp)*sin(a*2*pi);
    else
        f = 0;
    end
    v = v0*f;
end