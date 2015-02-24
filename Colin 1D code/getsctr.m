function [sctr]= getsctr(e,IEN)

% initialize
sctr= [];

% scatter matrix for regular dof
sctr(1)= IEN(e,1);
sctr(2)= IEN(e,2);



