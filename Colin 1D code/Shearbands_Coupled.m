function Shearbands_Coupled()
globvarssb;
vo = 10;
tsteps = 2000;
totaltime = 2/5048;
timeramp = totaltime/10;%1.0e-6;
dt = totaltime/tsteps;


%input options
       
%Jacobian = 'Numerical';
Jacobian = 'Analytic';
material_model = 'Modified Litonski, Exponential Thermal Softening';
integration_rule = 'Backward Euler';
%integration_rule = 'Foreward Euler';
%integration_rule = 'Trapezoidal Rule';

%{
L = 1;
%{
pl = 5;
x = zeros(nnode,1);
for n = 1:nnode
    x(n) = 16*L*((n-nnode/2-1/2)/nel)^pl;
end
%}
% set mesh
% ---------
%{
ldomain1 = .0002;
ldomain2 = .000098;
ldomain3 = .000004;
ldomain4 = .000098;
ldomain5 = .0002;
ldomain6 = .000098;
ldomain7 = .000004;
ldomain8 = .000098;
ldomain9 = .0002;


domain=[ldomain1, ldomain2, ldomain3, ldomain4, ldomain5, ldomain6, ldomain7, ldomain8, ldomain9]; % size of domain

%}
% ldomain1 = .0004;
% ldomain2 = .000098;
% ldomain3 = .000004;
% ldomain4 = .000098;
% ldomain5 = .0004;
% 
% 
% 
% domain=[ldomain1, ldomain2, ldomain3, ldomain4, ldomain5]; % size of domain
% create 1d mesh
% ---------------
%}
% [x,IEN,nel,nnode]= get1dmesh(domain,divide,L);

L = 1;
nel = 100;
nnode = nel+1;
x = 0:(1/nnode):L;
IEN = [1:nel; 2:nnode]';

%load disc3;
neq = nnode;
%x   =  linspace(0,L,nnode);
options{1} = ['Integration Rule = ',integration_rule];
options{2} = ['Jacobian = ',Jacobian];
options{3} = ['Material Model = ',material_model];
options{4} = ['Number of Elements = ',num2str(nel)];
options{5} = ['Applied Velocity = ', num2str(vo)];

par.integration_rule = integration_rule;
par.Jacobian = Jacobian;
par.material_model = material_model;
[prop,options,imper] = getprop(par,options);
proplist;


% set velocity BC
% --------------
nvbcnode = 2; % total number of velocity BC node


vbcnode(1) = 1; % node #
vbc(1) = 0; % prescribed velocity

vbcnode(2) = nnode; 
vbc(2) = vo;
% connectivity array

date_and_time = datestr(now);


%{
if Save_options == 1
    fid_options = fopen(['Shearband1D_Implicit_Options_',date_and_time,'.txt'],'w');
    for i = 1:length(options)
        fprintf(fid_options,'%s \n',options{i}); %minimum number of digits to be printed is 4, number of digits
    end
    fclose(fid_options);
end

if Save_mesh ==1
    fid_mesh = fopen(['Shearband1D_Implicit_mesh_',date_and_time,'.txt'],'a');    
    fprintf(fid_mesh,' %e',x); %minimum number of digits to be printed is 4, number of digits    
    fclose(fid_mesh);    
end
%}

%% Initial values
V = zeros(neq,1);
plastic_strain = zeros(nel,1);
stress =  zeros(nel,1);
T = 293*ones(neq,1);
d = zeros(neq,1);

%% system matrices

% Initializations
fv0     = zeros(neq,1);
fsig0   = zeros(nel,1);
fT0     = zeros(neq,1);
fgamp0  = zeros(nel,1);

Mv      = zeros(neq,neq);
Kv      = zeros(neq,nel);
Msig    = zeros(nel,nel);
Ksig    = zeros(nel,neq);
MT      = zeros(neq,neq);
KT      = zeros(neq,neq);
Mgamp   = zeros(nel,nel);

% Assembly
for e = 1:nel
    [ sctr ]= getsctr( e, IEN );
    [ mve, kve, msige, ksige, mTe, kTe, mgampe ] = elementmatrix1D(e,prop,IEN,x);
    
    Mv(sctr(:),sctr(:)) = Mv(sctr(:),sctr(:)) + mve;
    Kv(sctr(:),e) = Kv(sctr(:),e) + kve;
    Msig(e,e) = Msig(e,e) + msige;
    Ksig(e,sctr(:)) = Ksig(e,sctr(:)) +ksige;
    MT(sctr(:),sctr(:)) = MT(sctr(:),sctr(:)) + mTe;
    KT(sctr(:),sctr(:)) = KT(sctr(:),sctr(:)) + kTe;
    Mgamp(e,e) = Mgamp(e,e) + mgampe;
    
    [fve,fsige,fTe,fgampe] = get_forcings(e,stress(e),plastic_strain(e),T(sctr(:)),prop,IEN,x,imper);
    fv0(sctr(:)) = fv0(sctr(:)) + fve;
    fsig0(e) = fsig0(e) + fsige;
    fT0(sctr(:)) = fT0(sctr(:)) + fTe;
    fgamp0(e) = fgamp0(e) + fgampe;
    
end
X0 = vertcat(V,stress,T,plastic_strain);
d0 = d;
tol = 1e-5;
time = 0;
% save init_vec_mat Mv Kv Msig Ksig MT KT Mgamp fv0 fsig0 fT0 fgamp0;
% return

%% Loop in time

luc_count = 0;

while time < totaltime

    time = time+dt;
    Velbc = zeros(nvbcnode,1);
    for inode=1:nvbcnode
        if ( time <= timeramp ) % within ramp function region
            V(vbcnode(inode)) = (pi/timeramp)*sin(time/timeramp*2*pi)*vbc(inode);
            Velbc(inode) = (pi/timeramp)*sin(time/timeramp*2*pi)*vbc(inode);
        else
            V(vbcnode(inode)) = vbc(inode);
            Velbc(inode) = vbc(inode);
        end
    end
    
    %{
    deformation_rate = Ksig*V;
    [stress,T,plastic_strain] = RTMM(deformation_rate,stress,T,plastic_strain,dt,prop);
    %}
    X = X0;
    V0 = X0(1:nnode);
    [Residual,fv_new,fsig_new,fT_new,fgamp_new,normres] = evaluate_residual(X,nnode,nel,dt,IEN,Mv,Kv,Msig,Ksig,MT,KT,Mgamp,V,stress,T,plastic_strain,fv0,fsig0,fT0,fgamp0,prop,Velbc,x,imper);
    JAC = assemble_jacobian(X,nel,nnode,Mv,Kv,Ksig,Msig,MT,KT,Mgamp,prop,dt,IEN,x,imper,Velbc(2));
    disp(['Time = ',num2str(time,'%10.4e'),', Residual = ',num2str(norm(Residual),'%10.4e')])
    Newton_count = 0;
    
    %%%%%%%%%%%%%%%%%%%%%
    % Newton Iterations %
    %%%%%%%%%%%%%%%%%%%%%
    while normres > tol 
        
        Newton_count = Newton_count + 1
        if Newton_count > 10
            dt = dt/2;
            X = X0;
            %[Residual] = evaluate_residual(X,nnode,nel,dt,IEN,Mv,Kv,Msig,Ksig,MT,KT,Mgamp,V,stress,T,plastic_strain,fv0,fsig0,fT0,fgamp0,prop,Velbc,x,imper);
            JAC = assemble_jacobian(X,nel,nnode,Mv,Kv,Ksig,Msig,MT,KT,Mgamp,prop,dt,IEN,x,imper,Velbc(2));
            Newton_count = 1;
            if dt < 5e-11
                return
            end
        end
        %Compute search direction with a newton step
        delX = JAC \ (-1*Residual);
        
        [Residual,fv_new,fsig_new,fT_new,fgamp_new,normres] = evaluate_residual(X+delX,nnode,nel,dt,IEN,Mv,Kv,Msig,Ksig,MT,KT,Mgamp,V,stress,T,plastic_strain,fv0,fsig0,fT0,fgamp0,prop,Velbc,x,imper);
        X = X + delX;
        %Simple Line search, if needed
        %{
        step_length = 1;
        if norm(Residual_trial) < norm(Residual)
            X = X+delX;
            Residual = Residual_trial;
        else
            while norm(Residual_trial) >= norm(Residual)
                step_length = .5*step_length;
                if step_length < 5*eps
                    return
                end
                [Residual_trial,fv_new,fsig_new,fT_new,fgamp_new] = evaluate_residual(X+step_length*delX,nnode,nel,dt,IEN,Mv,Kv,Msig,Ksig,MT,KT,Mgamp,V,stress,T,plastic_strain,fv0,fsig0,fT0,fgamp0,prop,Velbc,x);
            end
            X = X + step_length*delX;
            Residual = Residual_trial;
        end
        %}
        disp(['Time = ',num2str(time,'%10.4e'),', Residual = ',num2str(normres,'%10.4e'),', Time Step = ',num2str(dt,'%10.4e')])
        JAC = assemble_jacobian(X,nel,nnode,Mv,Kv,Ksig,Msig,MT,KT,Mgamp,prop,dt,IEN,x,imper,Velbc(2));
        
    end
    
    luc_count = luc_count + 1;
%     if luc_count == 5
%         save jacobian JAC Residual;
%         return;
%     end
    
    %%%%%%%%%%%%%%
    % End Newton %
    %%%%%%%%%%%%%%

  %  [Residual,JAC,fv_new,fsig_new,fT_new,fgamp_new] = evaluate_residual_FC(X,nnode,nel,dt,IEN,Mv,Kv,Msig,Ksig,MT,KT,Mgamp,V,stress,T,plastic_strain,fv0,fsig0,fT0,fgamp0,prop,Velbc);

    V = X(1:nnode);
    
    figure(1)
    plot(V)
    axis([0 101 -8E5 8E5])
    pause(.001)
            
    stress = X(nnode+1:nnode+nel);
    T = X(nnode+nel+1:2*nnode+nel);
    plastic_strain = X(nel+2*nnode+1:2*nel+2*nnode);
    d = d0 + dt*V0 + .5*dt*(V-V0);
    
    X0 = X;
    d0 = d;
    fv0 = fv_new;
    fsig0 = fsig_new;
    fT0 = fT_new;
    fgamp0 = fgamp_new;
    
end