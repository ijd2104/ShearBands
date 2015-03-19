function X = linearElasticity()
%LINEARELASTICITY Solves linear elasticity problem using mixed FEM
%   Saves velocity and displacement at intervals
    global t N
    
    setupData();
    getMesh();
    
    i = 10; %save interval
    X = num2cell(zeros(N.node,t.steps/i),1);
    X0 = zeros(N.node,1);
    
    Xt = X0;
    for n = 1:t.steps
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
    global t TimeIntPar NewtonPar
    
    %Time data
    t.total = 2E-4; %Time for wave to travel 1m, 1/5000
    t.steps = 2000;
    t.dt = t.total/t.steps;
    t.ramp = t.total/4;
    t.curr = 0;
    
    %Newton's method parameters
    NewtonPar.NormTol = 1E-8;
    NewtonPar.iter = 500;
    
    %Time integration parameters
    TimeIntPar.alpha = 1;
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
    
    X0 = Xt;
    [J,R] = matrixAssembly(Xt,X0);
    [J,R,X] = matrixPartition(J,R,X0);
    X.BC = get_vbc(X);
    dX.D = -X.D+X.BC;
    dX.N = J.NN\(-R.N+J.ND*dX.D);
    Xk = [X.BC; X.N+dX.N];
    Xk = unPartition(Xk);
    
    for k = 1:niter
        [J,R] = matrixAssembly(Xt,Xk);
        if norm(R) < ntol
            break
        end
        [J,R,X] = matrixPartition(J,R,Xk);
        dX.N = -J.NN\R.N;
        Xk = [X.BC; X.N+dX.N];
        Xk = unPartition(Xk);
    end
end

function [J,R] = matrixAssembly(Xt,Xn)
    global N
    
    %Initialize matrices
    R.v = zeros(N.vnode,1);
    R.s = zeros(N.elem,1);
    J.vv = zeros(N.vnode,N.vnode);
    J.vs = zeros(N.vnode,N.elem);
    J.ss = zeros(N.elem,N.elem);
    J.sv = zeros(N.elem,N.vnode);
    
    %Assemble matrices
    for e = 1:N.elem
        elementMatrix(e); %returns [m,k]
        
        L = N.conn(e,:);
        xv = [Xt(L) Xn(L)];
        xs = [Xt(e+N.vnode) Xn(e+N.vnode)];
        
        r = elementResidual(xv,xs);
        R.v(L) = R.v(L)+r.v;
        R.s(e) = R.s(e)+r.s;
        
        j = elementJacobian(xv,xs);
        J.vv(L,L) = J.vv(L,L)+j.vv;
        J.vs(L,e) = J.vs(L,e)+j.vs;
        J.sv(e,L) = J.sv(e,L)+j.sv;
        J.ss(e,e) = J.ss(e,e)+j.ss;
    end
    R = [R.v; R.s];
    J = [J.vv J.vs; J.sv J.ss];
    J = sparse(J);
end

function elementMatrix(e)
    global N u
    global m k
    
    E = 200E9; %Pa
    rho = 7850; %kg/m^3
    
    L = N.conn(e,:);
    ue = u(L);
    h = abs(ue(2)-ue(1));
    J = h/2;
    
    Nv = @(x) (1/2)*[1-x 1+x];
    Bv = (1/h)*[-1 1];
    Ns = [1];
    fvv = @(x) rho*Nv(x)'*Nv(x);
    fvs = @(x) Bv'*Ns;
    fsv = @(x) -E*Ns'*Bv;
    fss = @(x) Ns'*Ns;
    
    m.vv=0; m.ss=0;
    k.vs=0; k.sv=0;
    ngp = 2; %number of Gauss points
    [W,xi] = gaussQuad(ngp);
    for i = 1:ngp
        m.vv = m.vv+W(i)*fvv(xi(i));
        k.vs = k.vs+W(i)*fvs(xi(i));
        k.sv = k.sv+W(i)*fsv(xi(i));
        m.ss = m.ss+W(i)*fss(xi(i));
    end
    m.vv = m.vv*J;
    k.vs = k.vs*J;
    k.sv = k.sv*J;
    m.ss = m.ss*J;
end

function r = elementResidual(xv,xs)
    global TimeIntPar t
    global m k
    
    xv = num2cell(xv,1);
    xs = num2cell(xs,1);
    a = TimeIntPar.alpha;
    dt = t.dt;
    n = 1;
    
    r.v = (m.vv/dt)*(xv{n+1}-xv{n})-(1-a)*k.vs*xs{n}-a*k.vs*xs{n+1};
    r.s = (m.ss/dt)*(xs{n+1}-xs{n})-(1-a)*k.sv*xv{n}-a*k.sv*xv{n+1};
end

function j = elementJacobian(xv,xs)
    global TimeIntPar t
    global m k
    
    a = TimeIntPar.alpha;
    dt = t.dt;
    
    j.vv = (1/dt)*m.vv;
    j.vs = -a*k.vs;
    j.sv = -a*k.sv;
    j.ss = (1/dt)*m.ss;
end

function [J,R,X] = matrixPartition(J0,R0,X0)
    global N BC eN eD
    
    BC = [1 0; N.vnode 1];
    eN = [1:N.node]';
    eD = BC(:,1);
    for i = 1:numel(eD)
        eN(eN==eD(i)) = [];
    end
    
    X.D = X0(eD);
    X.N = X0(eN);
    J.DD = J0(eD,eD);
    J.DN = J0(eD,eN);
    J.ND = J0(eN,eD);
    J.NN = J0(eN,eN);
    R.D = R0(eD);
    R.N = R0(eN);
    
    X.BC = zeros(size(BC,1),1);
    for i = 1:size(BC,1)
        X.BC(i) = BC(i,2);
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
    v = X.BC*vt;
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
