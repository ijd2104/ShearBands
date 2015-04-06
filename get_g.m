function g = get_g( input_args )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    g = gampdotr*exp((s/((1-((T-T0)/(Tm-T0))^m)*(As+Bs*gamp^N))-1)/c);
    dgds = -g/(c*(1-((T-T0)/(Tm-T0))^m)*(As+Bs*gamp^N));
    dgdp = dgds*Bs*Ns*gamp^(N-1)/(As+Bs*gamp^N);
    dgdT = dgds*m*s*((To-T)/(T0-Tm))^m/((T-T0)*((T0-T)/(T0-Tm))^m);
end

