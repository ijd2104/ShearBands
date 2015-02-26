function x = shearBands()
%SHEARBANDS Solves linear elasticity

    global t N
    
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
        X{n+1} = newtonIter(X{n});
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
    
    y = [-1 0 1];
    N.elem = 2;
    N.node = N.elem+1;
    N.eqn = N.node;
    N.conn = [1:N.node-1;2:N.node]';
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

function [J,R] = matrixAssembly(Xn,Xk)
    global N
    
    %Initialize matrices
    R.v = zeros(N.eqn,1);
    R.s = zeros(N.elem,1);
    J.vv = zeros(N.eqn,N.eqn);
    J.vs = zeros(N.eqn,N.elem);
    J.ss = zeros(N.elem,N.elem);
    J.sv = zeros(N.elem,N.eqn);

    for e = 1:N.elem
        L = N.conn(e,:);
        xv = num2cell([Xn(L) Xk(L)],1);
        xs = num2cell([Xn(e+N.node) Xk(e+N.node)],1);
        elementMatrix(e);
        
        r = elementResidual(xv,xs);
        R.v(L) = R.v(L)+r.v;
        R.s(e) = R.s(e)+r.s;
        
        j = elementJacobian();
        J.vv(L(:),L(:)) = J.vv(L(:),L(:))+j.vv;
        J.vs(L(:),e) = J.vs(L(:),e)+j.vs;
        J.sv(e,L(:)) = J.sv(e,L(:))+j.sv;
        J.ss(e,e) = J.ss(e,e)+j.ss;
    end
    R = [R.v; R.s];
    J = [J.vv J.vs; J.sv J.ss];
    J = sparse(J);
end

function Xk = newtonIter(Xn)
    global N NewtonPar
    
    niter = NewtonPar.iter;
    ntol = NewtonPar.NormTol;
    nnode = 2*N.elem+1;
    
    Xk = Xn;
    
    for k = 1:niter
        
        [J,R] = matrixAssembly(Xn,Xk);
        %Jacobian
        if condest(J) == Inf
            J = J+diag(ones(1,nnode)*1E-7);
        end
        %Residual
        if normest(R) < ntol
            break
        end
        dXk = -J\R;
        
        Xk = Xk+dXk;
    end

end