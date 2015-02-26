function x = shearBands()
%SHEARBANDS Solves linear elasticity

    global t X N
    
    setupData();
    getMesh();
    
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
end

function v = get_vbct()
    global t
    v = sin(t.curr);
end

function elementMatrix(e)
%Generates element matrix using two point Gauss quadrature
    global N y
    global m k
    
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
    m.vv = J*(W(1)*f(xi(1))+W(2)*f(xi(2)));
    
    f = @(xi) Bv'*Ns;
    k.vs = J*(W(1)*f(xi(1))+W(2)*f(xi(2)));
    
    f = @(xi) -E*Ns'*Bv;
    k.sv = J*(W(1)*f(xi(1))+W(2)*f(xi(2)));
    
    f = @(xi) Ns'*Ns;
    m.ss = J*(W(1)*f(xi(1))+W(2)*f(xi(2)));
end

% function matrixAssembly()
%     global N M K
%     
%     % Initializations
%     M.vv = zeros(N.eqn,N.eqn);
%     K.vs = zeros(N.eqn,N.elem);
%     M.ss = zeros(N.elem,N.elem);
%     K.sv = zeros(N.elem,N.eqn);
%     
%     %Assembly
%     for e = 1:N.elem
%         L = N.conn(e,:);    %scatter matrix
%         [m,k] = elementMatrix(e);
%         M.vv(L(:),L(:)) = M.vv(L(:),L(:))+m.vv;
%         K.vs(L(:),e) = K.vs(L(:),e)+k.vs;
%         M.ss(e,e) = M.ss(e,e)+m.ss;
%         K.sv(e,L(:)) = K.sv(e,L(:))+k.sv;
%     end
%     M.vv = sparse(M.vv);
%     K.vs = sparse(K.vs);
%     M.ss = sparse(M.ss);
%     K.sv = sparse(K.sv);
% end

function r = elementResidual(xv,xs)
    global IntPar t
    global m k
    
    a = IntPar.alpha;
    
    n = 1;
    
    dt = t.dt;
    
    r.v = (m.vv/dt)*(xv{n+1}-xv{n})-(1-a)*k.vs*xs{n}-a*k.vs*xs{n+1};
    r.s = (m.ss/dt)*(xs{n+1}-xs{n})-(1-a)*k.sv*xv{n}-a*k.sv*xv{n+1};
end

function j = elementJacobian()
    global t IntPar
    global m k
    
    a = IntPar.alpha;
    dt = t.dt;
    
    j.vv = (1/dt)*m.vv;
    j.vs = -a*k.vs;
    j.sv = -a*k.sv;
    j.ss = (1/dt)*m.ss;
end

function [J,R] = matrixAssembly(xt,xn)
    global N
    
    R.v = zeros(N.eqn,1);
    R.s = zeros(N.elem,1);
    J.vv = zeros(N.eqn,N.eqn);
    J.vs = zeros(N.eqn,N.elem);
    J.ss = zeros(N.elem,N.elem);
    J.sv = zeros(N.elem,N.eqn);

    for e = 1:N.elem
        L = N.conn(e,:);
        xv = num2cell([xt(L) xn(L)],1);
        xs = num2cell([xt(e+N.node) xn(e+N.node)],1);
        elementMatrix(e);
        
        r = elementResidual(xv,xs);
        R.v(L) = R.v(L)+r.v;
        R.s(e) = R.s(e)+r.s;
        
        j = elementJacobian();
        J.vv(L(:),L(:)) = J.vv(L(:),L(:))+j.vv;
        J.vs(L(:),e) = J.vs(L(:),e)+j.vs;
        J.ss(e,e) = J.ss(e,e)+j.ss;
        J.sv(e,L(:)) = J.sv(e,L(:))+j.sv;
    end
    R = [R.v; R.s];
    J = [J.vv J.vs; J.sv J.ss];
    J = sparse(J);
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
        
        [J,R] = matrixAssembly(x0,x{k});

        if condest(J) == Inf
            J = J+diag(ones(1,nnode)*1E-7);
        end
        
        if normest(R) < ntol
            break
        end
        dx = -J\R;
        x{k+1} = x{k}+dx;
    end
    x = x{k};
end