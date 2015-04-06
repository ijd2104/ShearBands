function dgdp = get_dgdp( g )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    dgdp = Bs*Ns*gamp^(n-1)*g/(c*(1-((T-T0)/(Tm-T0))^m)*(As+Bs*gamp^n))^2;

end