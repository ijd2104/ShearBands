function x = shearBands()
%SHEARBANDS Solves linear elasticity

    global t X N
    
    setupData();
    getMesh();
    matrixAssembly();
    
    nnode = 2*N.elem+1;
    
    v0 = 10;
    X0 = [v0 zeros(1,N.elem-1) -v0 zeros(1,N.elem)]';
    % Advance in time
    X = [X0 zeros(nnode,t.steps)];
    X = num2cell(X,1);
    for n = 1:t.steps
        t.iter = n;
        disp(t.iter/t.steps)
        X{n+1} = newtonIter();
    end
    x = cell2mat(X);
end

function setupData()
    global t NewtonPar IntPar
    
    t.total = 2.5E-6;
    t.dt = 5E-9;
    t.steps = floor(t.total/t.dt);
    t.curr = 0;
    
    %Newton data
    NewtonPar.NormTol = 1E-8;
    NewtonPar.iter = 500;
    
    %Time Inegration
    IntPar.alpha = 1;
end

function getMesh()
%GETMESH Generates a 1D mesh
     global N y
     L = 1E-3;   %length of bar
    divs = [10 30 21 30 10]';   %number of domain subdivisions
    doms = [.0004 .000098 .000004 .000098 .0004]';    %upper bounds for domains
    
    N.elem = sum(divs);
    N.node = N.elem+1;
    N.eqn = N.node;
    N.conn = [1:N.node-1;2:N.node]';
    
    %Node positions
    y = cell(1,numel(doms));
    for i = 0:numel(doms)-1
        coords = linspace(sum(doms(1:i)),sum(doms(1:i+1)),divs(i+1)+1)-L/2;
        if i == 0
            y{i+1} = coords;
        else
            y{i+1} = coords(2:divs(i+1)+1);
        end
    end
    y = cell2mat(y);
% y = [5E-4 0 5E-4];
% N.elem = 2;
% N.node = N.elem+1;
% N.eqn = N.node;
% N.conn = [1:N.node-1;2:N.node]';
end

function [Mvv,Kvs,Mss,Ksv] = elementMatrix(e)
%Generates element matrix using two point Gauss quadrature
    global N y
    
    L = N.conn(e,:); %Scatter matrix
    ye = y(L);

    E = 200E9;
    rho = 7830;
    
    h = abs(ye(2)-ye(1));
    J = h/2;
    
    W = [1 1];
    xi = [-0.577350269189626 0.577350269189626];
    
    %Shape functions and derivatives
    Nv = @(xi) (1/2)*[1-xi 1+xi];
    Bv = (1/h)*[-1 1];
    Ns = [1];
    
    f = @(xi) rho*Nv(xi)'*Nv(xi);
    Mvv = J*(W(1)*f(xi(1))+W(2)*f(xi(2)));
    
    f = @(xi) Bv'*Ns;
    Kvs = J*(W(1)*f(xi(1))+W(2)*f(xi(2)));
    
    f = @(xi) -E*Ns'*Bv;
    Ksv = J*(W(1)*f(xi(1))+W(2)*f(xi(2)));
    
    f = @(xi) Ns'*Ns;
    Mss = J*(W(1)*f(xi(1))+W(2)*f(xi(2)));
end

function matrixAssembly()
    global N M K
    % Initializations
    M.vv = zeros(N.eqn,N.eqn);
    K.vs = zeros(N.eqn,N.elem);
    M.ss = zeros(N.elem,N.elem);
    K.sv = zeros(N.elem,N.eqn);
    
    %Assembly
    for e = 1:N.elem
        L = N.conn(e,:);    %scatter matrix
        [mvv,kvs,mss,ksv] = elementMatrix(e);
        M.vv(L(:),L(:)) = M.vv(L(:),L(:))+mvv;
        K.vs(L(:),e) = K.vs(L(:),e)+kvs;
        M.ss(e,e) = M.ss(e,e)+mss;
        K.sv(e,L(:)) = K.sv(e,L(:))+ksv;
    end
    M.vv = sparse(M.vv);
    K.vs = sparse(K.vs);
    M.ss = sparse(M.ss);
    K.sv = sparse(K.sv);
end

function x = newtonIter()
    global N NewtonPar t X
    
    niter = NewtonPar.iter;
    ntol = NewtonPar.NormTol;
    nnode = 2*N.elem+1;
    
    x0 = X{t.iter};
    x = [x0 zeros(nnode,niter)];
    x = num2cell(x,1);
    
    for k = 1:niter
        
        J = assembleJacobian();
        if condest(J) == Inf
            J = J+diag(ones(1,nnode)*1E-7);
        end
        
        R = evaluateResidual(x{k});
        if normest(R) < ntol
            break
        end
        dx = -J\R;
        x{k+1} = x{k}+dx;
    end
    x = x{k};
end

function R = evaluateResidual(x)
    global t X M K IntPar N
    a = IntPar.alpha;
    
    xt = X{t.iter};
    xv = num2cell([xt(1:N.node) x(1:N.node)],1);
    xs = num2cell([xt(N.node+1:end) x(N.node+1:end)],1);
    n = 1;
    
    dt = t.dt;
    
    rv = (M.vv/dt)*(xv{n+1}-xv{n})-(1-a)*K.vs*xs{n}-a*K.vs*xs{n+1};
    rs = (M.ss/dt)*(xs{n+1}-xs{n})-(1-a)*K.sv*xv{n}-a*K.sv*xv{n+1};
    R = [rv;rs];
end

function J = assembleJacobian()
    global t M K IntPar
    a = IntPar.alpha;
    dt = t.dt;
    
    Jvv = (1/dt)*M.vv;
    Jvs = -a*K.vs;
    Jsv = -a*K.sv;
    Jss = (1/dt)*M.ss;
    J = [Jvv Jvs; Jsv Jss];
end