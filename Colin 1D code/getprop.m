function [prop,options,imper] = getprop(par,options)

switch par.material_model
    case 'Modified Litonski, Exponential Thermal Softening'
        
        imper = .0005;
        prop(1)=  200.e9;               % e : young's modulus
        prop(2)= 7830.0d0;              % rho0 : mass density
        prop(3)=  448.0d0;              % cp : specific heat
        prop(4)=  11.2e-6;              % thet : coefficient of thermal expansion
        prop(5)=  0.90d0;               % kai : the fraction of plastic work converted to heat
        prop(6)= 0;                     % heat transfer coefficent for convection boundary condition
        prop(7)= 1.0d0;                 % cross section area of 1d model
        prop(8)= 803.5;                 % thermal conductivity
        switch par.integration_rule
            case 'Backward Euler'
                prop(9) = 1;            % integration parameter theta
            case 'Foreward Euler'
                prop(9) = 0;
            case 'Trapezoidal Rule'
                prop(9) = .5;
        end
        prop(10) = .99;                 % imperfection percent
        prop(11) = length(imper);       % number of imperfections
        prop(12) = 5e-5;                % imperfection radius
        prop(13)=  1.e-3;               % eps0dot : reference strain rate
        prop(14)=  70.0d0;              % m: rate senstivity parameter
        prop(15)=  2000.e6;             % sig0 : yield stress
        prop(16)=  2000.e6/200.e9;      % eps0 : yield strain
        prop(17)=  0.010d0;             % n : strain hardening exponent
        prop(18)=  293.0d0;             % t0 : reference temperature
        prop(19)=  0.80d0;              % del : thermal softening parameter
        prop(20)=  500.0d0;             % k : thermal softening parameter

        
        a = length(options);
        options{a+1}  = 'Material Properties:';
        options{a+2}  = ['Youngs Modulus, E = ',num2str(prop(1))];
        options{a+3}  = ['Mass Density, \\rho = ',num2str(prop(2))];
        options{a+4}  = ['Specific Heat, \\c_p = ',num2str(prop(3))];
        options{a+5}  = ['Coefficient of Thermal Expansion, \\alpha = ',num2str(prop(4))];
        options{a+6}  = ['Taylor-Quinney Coefficient, \\kai = ',num2str(prop(5))];
        options{a+7}  = ['Heat Transfer Coefficient, h = ',num2str(prop(6))];
        options{a+8}  = ['Cross Sectional Area, A = ',num2str(prop(7))];
        options{a+9}  = ['Thermal Conductivity,  \\kappa = ',num2str(prop(8))];
        options{a+10} = ['Integration Parameter  \\theta = ',num2str(prop(9))];
        options{a+11} = ['Imperfection Magnitude = ',num2str(prop(10))];
        options{a+12} = ['Imperfection Center= ',num2str(prop(11))];
        options{a+13} = ['Imperfection Radius,  \\kappa = ',num2str(prop(12))];              
        options{a+14} = ['Reference Strain Rate, \\dot{\epsilon_0} = ',num2str(prop(13))];
        options{a+15} = ['Rate Sensitivity Parameter, m = ',num2str(prop(14))];
        options{a+16} = ['Yield Stress, \\sigma_0 = ',num2str(prop(15))];
        options{a+17} = ['Yield Strain, \\epsilon_0 = ',num2str(prop(16))];
        options{a+18} = ['Strain Hardening Exponent, n = ',num2str(prop(17))];
        options{a+19} = ['Reference Temperature, \\T_0 = ',num2str(prop(18))];
        options{a+20} = ['Thermal Softening Parameter, \\delta = ',num2str(prop(19))];
        options{a+21} = ['Thermal Softening Parameter, k = ',num2str(prop(20))];

        
    case 'Modified Litonski, Linear Thermal Softening'
    case 'Johnson-Cook'
end